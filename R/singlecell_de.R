#' Run single-cell differential expression
#' 
#' Run differential expression using traditional single-cell methods
#' Note that this is effectively a wrapper of the FindMarkers function in Seurat.
#' 
#' @param input a single-cell matrix to be converted, with features (genes) in rows
#'   and cells in columns. Alternatively, a \code{Seurat}, \code{monocle3}, or 
#'   or \code{SingleCellExperiment} object can be directly input.
#' @param meta the accompanying meta data whereby the rownames match the column
#'   names of \code{input}. 
#' @param cell_type_col the vector in \code{meta} containing the cell type 
#'   information. Defaults to \code{cell_type}.
#' @param label_col the vector in \code{meta} containing the experimental
#'   label. Defaults to \code{label}. 
#' @param de_method the mixed model type to use. Defaults to wilcox.
#' @param min_cells the minimum number of cells in a cell type to retain it.
#'   Defaults to \code{3}.
#' @param min_features the minimum number of expressing cells (or replicates) 
#'   for a gene to retain it. Defaults to \code{0}.   
#' @param binarization binarization for single-cell ATAC-seq only
#' @param normalization normalization for Seurat/Signac methods
#' @param latent_vars latent variables for Seurat/Signac methods
#' @param input_type refers to either scRNA or scATAC
#' @return a data frame containing differential expression results.
#'  
#' @importFrom magrittr extract set_rownames %<>%
#' @importFrom Matrix rowSums
#' @importFrom dplyr %>% mutate bind_rows
#' @importFrom tibble rownames_to_column 
#' @import Seurat
#' 

singlecell_de = function(
        input,
        meta = NULL,
        cell_type_col = 'cell_type',
        label_col = 'label',
        de_method = 'wilcox',
        min_cells = 3,
        min_features = 0,
        normalization = 'log_tp10k',
        binarization = FALSE,
        latent_vars = NULL,
        input_type = 'scRNA'
) {
    
    if (input_type == 'scRNA') {
        # check the arguments
        if (!de_method %in% c("wilcox", "bimod", "t", "negbinom", "poisson", "LR", "MAST"))
            stop("Please select one of: wilcox, bimod, t, negbinom, poisson, LR,
           or MAST as the de_method")  
    } else if (input_type == 'scATAC') {
        if (!de_method %in% c("wilcox", "t", "negbinom", "LR", "fisher", "binomial", "LR_peaks", "permutation"))
            stop("Please select one of: wilcox, t, negbinom, LR, fisher, binomial, LR_peaks or permutation as the de_method")
    }
    
    if ((binarization == TRUE & normalization != 'log_tp10k')) {
        stop("Please choose either binarization effect or normalization effect. Running both at the same time seems irresponsible.")
    }
    
    if ((binarization == TRUE & input_type != 'scATAC')) {
        stop("Please select binarization for scATAC-seq only.")
    }
    
    # first, make sure inputs are correct
    inputs = check_inputs(
        input, 
        meta,
        replicate_col = NULL,
        cell_type_col = cell_type_col,
        label_col = label_col
    )
    expr = inputs$expr
    meta = inputs$meta
    
    # define labels, reverse to avoid logFC conflicts
    labels = unique(meta$label)
    if (is.factor(labels)) {
        label1 = levels(labels)[2]
        label2 = levels(labels)[1]
    } else {
        label1 = labels[2]
        label2 = labels[1]
    }
    
    # get cell types
    cell_types = unique(meta$cell_type)
    
    # create a seurat object
    rownames(meta) = colnames(expr)
    sc = CreateSeuratObject(expr, meta.data = meta)
    
    # make sure idents are set right in the Seurat object
    Idents(sc) = sc$cell_type
    
    # check if integer or already normalized, normalize if needed
    mat = GetAssayData(sc, slot = 'counts')
    if ((sum(mat %% 1 == 0) == length(mat)) == T) {
        if (normalization == 'log_tp10k'){
            sc %<>% NormalizeData()
        } else if (normalization == 'tp10k'){
            sc %<>% NormalizeData(normalization.method='RC')
        } else if (normalization == 'log_tp_median'){
            median_cpc = median(colSums(mat))
            sc %<>% NormalizeData(scale.factor = median_cpc)
        } else if (normalization == 'tp_median'){
            median_cpc = median(colSums(mat))
            sc %<>% NormalizeData(normalization.method='RC', scale.factor = median_cpc)
        } else if (normalization == 'TFIDF'){
            
            # extracted from Signac https://github.com/stuart-lab/signac/blob/2ad6c3c9c0c8dd31f7e1433b2efd5050d8606f27/R/preprocessing.R#L621
            npeaks = Matrix::colSums(mat)
            tf = Matrix::tcrossprod(x = mat, y = Diagonal(x = 1 / npeaks))
            rsums = rowSums(mat)
            idf = ncol(mat)/rsums
            idf = log(1+idf) # since precomputed idf = FALSE (by default)
            norm_mat = Diagonal(length(idf), idf) %*% tf
            slot(object = norm_mat, name = "x") = log1p(slot(object = norm_mat, name = "x") * 1e4)
            colnames(norm_mat) = colnames(mat)
            rownames(norm_mat) = rownames(mat)
            vals = slot(object = norm_mat, name = "x")
            vals[is.na(x = vals)] = 0
            slot(object = norm_mat, name = "x") = vals
                
            if (packageVersion("SeuratObject") >= 5) {
                sc[['RNA']]$data = norm_mat
            } else {
                sc[['RNA']]@data = norm_mat
            }
        }
    } else {
        if (packageVersion("SeuratObject") >= 5) {
                sc[['RNA']]$data = norm_mat
            } else {
                sc[['RNA']]@data = norm_mat
            }
    }
    
    if (binarization) {
        mat = GetAssayData(sc, slot='counts')
        mat@x[mat@x > 0] = 1
            
        if (packageVersion("SeuratObject") >= 5) {
                sc[['RNA']]$data = norm_mat
            } else {
                sc[['RNA']]@data = norm_mat
            }
    }
    
    # run single cell DE using Seurat
    DE = list()
    for (cell_type_idx in seq_along(cell_types)) {
        cell_type = cell_types[cell_type_idx]
        message("[", cell_type_idx, "/", length(cell_types),
                "] working on cell type: ", cell_type, " ...")
        
        # check to make sure there are enough cells
        n_cells = table(meta$cell_type, meta$label)
        if (min(n_cells[cell_type, ]) < min_cells) {
            message(" .. not enough cells, skipping ...")
            next
        } else {
            tryCatch({
                # subset to the right cell type
                Idents(sc) = sc$cell_type
                sub = sc %>% subset(idents = cell_type)
                # drop genes below threshold
                keep = rowSums(sub) >= min_features
                sub = sub[keep,]
                # seurat-based methods
                if (!de_method %in% c("fisher", "binomial", "LR_peaks", "permutation", "snapatac")){
                    # run DE analysis
                    res = FindMarkers(sub, ident.1 = label1, ident.2 = label2,
                                      assay = 'RNA', min.pct = -Inf, 
                                      min.cells.feature = 0,
                                      min.cells.group = min_cells, 
                                      logfc.threshold = -Inf,
                                      group.by = 'label', 
                                      subset.ident = cell_type,
                                      test.use = de_method,
                                      latent.vars = latent_vars
                    ) %>%
                        rownames_to_column('gene') %>%
                        mutate(
                            cell_type = cell_type,
                            de_family = 'singlecell',
                            de_method = de_method,
                            de_type = 'singlecell'
                        )
                } else {
                    meta = sub@meta.data
                        if (packageVersion("SeuratObject") >= 5) {
                                mat = sc[['RNA']]$data
                            } else {
                                mat = sc[['RNA']]@data
                            }
                    res = da_function_wrapper(mat, meta, method, cell_type)
                }
                DE[[cell_type]] = res
            }, error = function(e) message(e))
        }
    }
    DE %<>% bind_rows()
    return(DE)
}

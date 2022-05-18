#' Calculate delta variance
#' 
#' Calculate delta variance from a single-cell matrix.
#' 
#' @param input a single-cell matrix to be converted, with features (genes) in rows
#'   and cells in columns. Alternatively, a \code{Seurat}, \code{monocle3}, or 
#'   or \code{SingleCellExperiment} object can be directly input.
#' @param meta the accompanying meta data whereby the rownames match the column
#'   names of \code{input}. If a \code{Seurat}, \code{monocle3} or
#'   \code{SingleCellExperiment} object is provided this can be null.
#' @param replicate_col the vector in \code{meta} containing the replicate 
#'   information. Defaults to \code{replicate}.
#' @param cell_type_col the vector in \code{meta} containing the cell type 
#'   information. Defaults to \code{cell_type}.
#' @param label_col the vector in \code{meta} containing the experimental
#'   label. Defaults to \code{label}. 
#' @param min_cells the minimum number of cells in a cell type to retain it.
#'   Defaults to \code{3}.
#' @param min_reps the minimum number of replicates in a cell type to retain it.
#'   Defaults to \code{2}.
#' @param min_features the minimum number of replicates expressing  a gene 
#'   to retain it. Defaults to \code{0}   
#' @return a list of pseudobulk matrices, for each cell type.
#'  
#' @importFrom magrittr %<>% extract
#' @importFrom dplyr %>% rename_ count group_by filter pull n_distinct distinct
#'   summarise
#' @importFrom purrr map map_int
#' @importFrom Matrix rowSums colSums
#' @importFrom matrixStats rowVars
#' @importFrom edgeR cpm
#' @importFrom stats setNames
#' @importFrom methods is
#' @export
#' 
calculate_delta_variance = function(input, 
                         meta = NULL,
                         replicate_col = 'replicate',
                         cell_type_col = 'cell_type',
                         label_col = 'label',
                         min_cells = 3,
                         min_reps = 2,
                         min_features = 0) {
  # first, make sure inputs are correct
  inputs = check_inputs(
    input, 
    meta,
    replicate_col = replicate_col,
    cell_type_col = cell_type_col,
    label_col = label_col)
  expr = inputs$expr
  meta = inputs$meta
  
  # keep only cell types with enough cells
  keep = meta %>%
    dplyr::count(cell_type, label) %>%
    group_by(cell_type) %>%
    filter(all(n >= min_cells)) %>%
    pull(cell_type) %>%
    unique()
  
  # process data into gene x replicate x cell_type matrices
  DV = keep %>%
    map( ~ {
      print(.)
      cell_type = .
      meta0 = meta %>% filter(cell_type == !!cell_type)
      expr0 = expr %>% extract(, rownames(meta0))
      # catch cell types without replicates or conditions
      if (n_distinct(meta0$label) < 2)
        return(NA)
      replicate_counts = distinct(meta0, label, replicate) %>%
        group_by(label) %>%
        summarise(replicates = n_distinct(replicate)) %>%
        pull(replicates)
      if (any(replicate_counts < min_reps))
        return(NA)
      
      # process data into gene X replicate X cell_type matrice
      mm = model.matrix(~ 0 + replicate:label, data = meta0)
      mat_mm = expr0 %*% mm
      keep_genes = rowSums(mat_mm > 0) > min_features
      # return NA if no genes are retained
      if (length(keep_genes) == 0)
        return(NA)
      mat_mm = mat_mm[keep_genes, ] %>% as.data.frame()
      mat_mm %<>% as.data.frame()
      colnames(mat_mm) = gsub("replicate|label", "", colnames(mat_mm))
      # drop empty columns
      keep_samples = colSums(mat_mm) > 0
      mat_mm %<>% extract(, keep_samples)
      
      # shuffle the replicates
      meta00 = meta0 %>%
        group_by(label) %>%
        mutate(replicate = sample(replicate))
      mm2 = model.matrix(~ 0 + replicate:label, data = meta00)
      mat_mm2 = expr0 %*% mm2
      keep_genes = rowSums(mat_mm2 > 0) > min_features
      # return NA if no genes are retained
      if (length(keep_genes) == 0)
        return(NA)
      mat_mm2 = mat_mm2[keep_genes, ] %>% as.data.frame()
      mat_mm2 %<>% as.data.frame()
      colnames(mat_mm2) = gsub("replicate|label", "", colnames(mat_mm2))
      # drop empty columns
      keep_samples = colSums(mat_mm2) > 0
      mat_mm2 %<>% extract(, keep_samples)
      
      # normalize each matrix
      norm1 = cpm(mat_mm)
      norm2 = cpm(mat_mm2)
      
      # calculate the variance of each matrix
      var1 = rowVars(norm1)
      var2 = rowVars(norm2)
      
      # return the delta variance with the gene information
      delta = var2 - var1
      out = data.frame(gene = rownames(mat_mm), DV = delta)
      return(out)
    }) %>%
    setNames(keep)
}

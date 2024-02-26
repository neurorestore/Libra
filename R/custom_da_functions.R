library(lmtest)
library(parallel)

fisher_exact_helper = function(temp, compar1_barcodes) {
  temp[temp>0]=1
  temp = data.frame('count'=temp, 'barcode'=names(temp)) %>%
    mutate('label'=ifelse(barcode %in% compar1_barcodes, 1, 0))
  temp = table(temp[,c('count', 'label')])
  if (nrow(temp) > 1){
    res = fisher.test(temp)
    return (c(res$p.value, res$estimate))    
  }
  else
    return(c(NA,NA))
}
  
get_fisher_exact_pval = function(mat, compar1_barcodes, compar2_barcodes){
  mat = mat[rowSums(mat) > 0,]
  # res_df = data.frame(t(apply(mat, 1, fisher_exact_helper, compar1_barcodes=compar1_barcodes)))
  res_df = mclapply(1:nrow(mat), mc.cores=16, function(x){
    fisher_exact_helper(mat[x,], compar1_barcodes)
  })
  res_df = data.frame(do.call(rbind, res_df))
  colnames(res_df) = c('p_val', 'test_statistic')
  rownames(res_df) = rownames(mat)
  res_df$p_val_adj = p.adjust(res_df$p_val, method='BH')
  return(res_df)
}

get_snapatac_exact_pval = function(mat, compar1_barcodes, compar2_barcodes){
  bcv=0.1
  mat = mat[rowSums(mat) > 0,]
  data_sum = data.frame('compar1'=rowSums(mat[,compar1_barcodes]),'compar2' = rowSums(mat[,compar2_barcodes]))
  group = factor(c(1,2))
  design = model.matrix(~group)
  d = DGEList(counts=data_sum, group = group)
  de = exactTest(d, dispersion=bcv^2)
  temp = de$table
  temp$test_statistic = -log(temp$PValue)
  return(temp)
}

get_binomial_pval = function(mat, compar1_barcodes, compar2_barcodes){
  mat = mat[rowSums(mat) > 0,]
  mat[mat>0] = 1
  mat1 = mat[,compar1_barcodes]
  mat2 = mat[,compar2_barcodes]
  
  n1 = ncol(mat1)
  n2 = ncol(mat2)
  s1 = rowSums(mat1)
  s2 = rowSums(mat2)
  m1 = s1/n1
  m2 = s2/n2

  pvals = do.call(rbind, mclapply(seq_along(m1), mc.cores=16, function(x){
    if(m1[x] >= m2[x]){
      test_stat = pbinom(q = s1[x]-1, size = n1, prob = (max(c(s2[x],1))) / n2, lower.tail=FALSE, log.p = TRUE)
    }else{
      test_stat = pbinom(q = s2[x]-1, size = n2, prob = (max(c(s1[x],1))) / n1, lower.tail=FALSE, log.p = TRUE)
    }
    p = min(c(2 * exp(test_stat), 1)) #handle 2-sided test
    return (c(p, test_stat))
  }))
  pvals = cbind(names(s1), pvals)
  colnames(pvals) = c('gene', 'p_val', 'test_statistic')
  pvals %<>%
    as_tibble() %>%
    mutate(
      p_val=as.numeric(p_val),
      test_statistic=as.numeric(test_statistic)
    )
  pvals$p_val_adj = p.adjust(pvals$p_val, method = 'BH')
  return(data.frame(pvals))
}

get_LR_peaks_pval = function(mat, compar1_barcodes, compar2_barcodes){
  mat = mat[rowSums(mat) > 0, ]
  no_of_features = nrow(mat)
  mat[mat > 0] = 1
  m1 = mat[,compar1_barcodes]
  m2 = mat[,compar2_barcodes]

  fmla = as.formula(object = "GENE ~ group")
  fmla2 = as.formula(object = "GENE ~ 1")
  pvals = do.call(rbind, mclapply(rownames(mat), mc.cores=1, function(x){
    test_df = rbind(
      data.frame('GENE'=m1[x,], 'group'=0),
      data.frame('GENE'=m2[x,], 'group'=1)
    )

    model1 = glm(formula = fmla, data = test_df, family = "binomial")
    model2 = glm(formula = fmla2, data = test_df, family = "binomial")
    lrtest = lrtest(model1, model2)
    out = c(lrtest$Pr[2], lrtest$Chisq[2])
    return(out)
  }))
  pvals = cbind(rownames(mat), pvals)
  colnames(pvals) = c('gene', 'p_val', 'test_statistic')
  pvals %<>%
    as_tibble() %>%
    mutate(
      p_val=as.numeric(p_val),
      test_statistic=as.numeric(test_statistic)
    )
  pvals$p_val_adj = p.adjust(pvals$p_val, method = 'BH')
  return(pvals)
}

get_permutation_test_pval = function(mat, compar1_barcodes, compar2_barcodes){

  n_trials=1000
  n_cores=1

  compar1_size = length(compar1_barcodes)
  compar2_size = length(compar2_barcodes)

  mat = mat[rowSums(mat) > 0,]
  m1 = mat[,compar1_barcodes]
  m2 = mat[,compar2_barcodes]

  p1_prime = sum(rowSums(m1 > 0))
  p2_prime = sum(rowSums(m2 > 0))

  combined_barcodes = c(compar1_barcodes, compar2_barcodes)

  permutations = do.call(cbind, mclapply(1:n_trials, mc.cores=ncores, function(i){

    random_compar1_barcodes = sample(combined_barcodes, compar1_size)
    random_compar2_barcodes = combined_barcodes[!combined_barcodes %in% random_compar1_barcodes]
    random_m1 = mat[,random_compar1_barcodes]
    random_m2 = mat[,random_compar2_barcodes]
    
    random_p1_prime = sum(rowSums(random_m1 > 0))
    random_p2_prime = sum(rowSums(random_m2 > 0))
    
    m1_random_means = rowMeans(random_m1 > 0)
    m2_random_means = rowMeans(random_m2 > 0)
    random_test_stats = (m1_random_means/random_p1_prime) - (m2_random_means/random_p2_prime)
    return (random_test_stats)
  }))

  colnames(permutations) = paste0('Trial_', 1:n_trials)
  rownames(permutations) = rownames(mat)
  permutations_means = rowMeans(permutations)
  permutations_sds = rowSds(permutations)
  permutations_res = data.frame('gene'=names(permutations_means), 'permutation_mean'=permutations_means, 'permutation_sd'=permutations_sds)

  p1_s = rowMeans(m1 > 0)
  p2_s = rowMeans(m2 > 0)
  test_stats = (p1_s/p1_prime) - (p2_s/p2_prime)
  test_stats_df = data.frame('gene'=names(test_stats), 'test_statistic'=test_stats)
  test_stats_df %<>% left_join(permutations_res, by='gene')
  pvals = test_stats_df %<>% mutate(
    z_score=(test_statistic-permutation_mean)/permutation_sd,
    p_val=pnorm(z_score)
    ) %>%
    select(gene, p_val, z_score)

  colnames(pvals) = c('gene', 'p_val', 'test_statistic')
  pvals %<>%
    as_tibble() %>%
    mutate(
      p_val=as.numeric(p_val),
      test_statistic=as.numeric(test_statistic)
    )
  pvals$p_val_adj = p.adjust(pvals$p_val, method = 'BH')
  return(pvals)
}

da_function_wrapper = function(mat, meta, da_method, cell_type_name) {
  comparisons = levels(meta$label)
  compar1_barcodes = meta %>%
    filter(label==comparisons[1]) %>%
    rownames()
  compar2_barcodes = meta %>%
    filter(label==comparisons[2]) %>%
    rownames()

  mean1 = log(rowMeans(mat[, compar1_barcodes]+1))
  mean2 = log(rowMeans(mat[, compar2_barcodes]+1))
  logFC_df = data.frame(
    gene = rownames(mat),
    avg_logFC = mean2 - mean1
  )

  if (da_method == 'fisher'){
    temp_res = get_fisher_exact_pval(mat, compar1_barcodes, compar2_barcodes)
    temp_res$gene = rownames(temp_res)
    rownames(temp_res) = NULL
    temp_res %<>%
      as_tibble() %>%
      mutate(
        cell_type = cell_type_name,
        de_family = 'singlecell',
        de_method = 'fisher',
        de_type = 'singlecell'
      ) %>% 
      left_join(logFC_df) %>%
      select(cell_type, gene, avg_logFC, test_statistic, p_val, p_val_adj, de_family, de_method, de_type)

  } else if (da_method == 'snapatac_findDAR'){
    temp_res = get_snapatac_exact_pval(mat, compar1_barcodes, compar2_barcodes)
    temp_res$gene = rownames(temp_res)
    temp_res$p_val_adj = p.adjust(temp_res$PValue, method='BH')
    rownames(temp_res) = NULL
    temp_res %<>%
      as_tibble() %>%
      mutate(
        cell_type = cell_type_name,
        de_family = 'snapatac_findDAR',
        de_method = da_method,
        de_type = 'snapatac_findDAR'
      ) %>% 
      dplyr::rename(p_val = PValue) %>%
      dplyr::rename(avg_logFC = logFC) %>%
      select(cell_type, gene, avg_logFC, test_statistic, p_val, p_val_adj, de_family, de_method, de_type)
  } else if (da_method == 'binomial'){
    temp_res = get_binomial_pval(mat, compar1_barcodes, compar2_barcodes)
    temp_res %<>%
      as_tibble() %>%
      mutate(
        cell_type = cell_type_name,
        de_family = 'singlecell',
        de_method = 'binomial',
        de_type = 'singlecell'
      ) %>% 
      left_join(logFC_df) %>%
      select(cell_type, gene, avg_logFC, test_statistic, p_val, p_val_adj, de_family, de_method, de_type)
  } else if (da_method == 'LR_peaks'){
      temp_res = get_LR_peaks_pval(mat, compar1_barcodes, compar2_barcodes, features_check=FALSE)
      temp_res %<>%
        as_tibble() %>%
        mutate(
          cell_type = cell_type_name,
          de_family = 'singlecell',
          de_method = 'LR_peaks',
          de_type = 'singlecell'
        ) %>% 
        left_join(logFC_df) %>%
        select(cell_type, gene, avg_logFC, test_statistic, p_val, p_val_adj, de_family, de_method, de_type)
    } else if (da_method == 'permutation'){
      temp_res = get_permutation_test_pval(mat, compar1_barcodes, compar2_barcodes)
      temp_res %<>%
        as_tibble() %>%
        mutate(
          cell_type = cell_type_name,
          de_family = 'singlecell',
          de_method = 'permuation',
          de_type = 'singlecell'
        ) %>% 
        left_join(logFC_df) %>%
        select(cell_type, gene, avg_logFC, test_statistic, p_val, p_val_adj, de_family, de_method, de_type)
    }
    return(temp_res)
}
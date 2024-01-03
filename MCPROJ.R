################################################################################
# 1. Seeding and mixture modeling
# Given a metacell footprint, represent the target metacell as a weighted combination 
# of the anchors
# https://cran.r-project.org/web/packages/CVXR/vignettes/cvxr_intro.html
################################################################################
# R implementation of MCProj
get_nn = function(query_mat,target_mat,k = 10){
  simmat = t(apply(query_mat,2,FUN = function(x) apply(target_mat,2,FUN = function(y) cor(x,y))))
  if(k>1){
    simplify = F
  }else{
    simplify = T
  }
  nn = apply(simmat,1,FUN = function(x) names(sort(x,decreasing = T))[1:k],simplify = simplify)
  return(nn)
}

mixture_model_mc = function(Y,X,eps = 1e-5,method = 'CVXR'){
  if(method == 'CVXR'){
    require('CVXR')
    # CAVE: this method takes longer!
    
    # Given the vector Y (query metacell), represent this vector as a weighted sum of X
    # returns: a list of projections
    #  Y = query_mc_fp[,q_id]
    #  X = target_mc_fp[,t_id]
    p = ncol(X)
    betaHat <- Variable(p)
    objective  <- Minimize(sum((Y - X %*% betaHat)^2))
    #problem <- Problem(objective)
    problem <- Problem(objective, constraints = list(betaHat >= eps,sum(betaHat) == 1))
    result <- solve(problem)
    betas = result$getValue(betaHat)
    if(!all(is.na(betas))){
      betas[betas<=eps] = 0
      rownames(betas) = colnames(X)
      # there are some genes unaccounted for indeed. 
      mc_proj = as.matrix(X) %*% as.matrix(betas)
      r2 = summary(lm(Y~as.matrix(X) %*% as.matrix(betas)))$adj.r.squared
    }else{
      mc_proj = NA
      betas = NA
      r2 = NA
    }
  }else if(method == 'nnls'){
    require('nnls')
    nnls.fit = nnls::nnls(X, Y)
    betas = as.matrix(nnls.fit$x,ncol = 1)
    rownames(betas) = colnames(X)
    betas[betas<=eps] = 0
    mc_proj = as.matrix(X) %*% as.matrix(betas)
    r2 = summary(lm(Y~as.matrix(X) %*% as.matrix(betas)))$adj.r.squared
    # rescale betas to sum to 1:
    #betas = betas/sum(betas)
  }
  return(list(proj = mc_proj,betas = betas, r2 = r2))
}

get_proj_matrix = function(projections){
  out = sapply(projections,FUN = function(x) x$proj)
  rownames(out) = rownames(projections[[1]]$proj)
  return(out)
}

# debug
if(F){
  # for some reason, doesn't work with e_gc values 
  k = 10 # number of anchors
  anchor = get_nn(query_mat,target_mat,k = 1)
  nn = get_nn(target_mat,target_mat,k = 10)
  # remove self-matches:
  # nn = setNames(lapply(seq_along(nn),FUN = function(i) nn[[i]][nn[[i]]!=names(nn)[i]]),names(nn))
  # extend the set of anchors
  mc_id = '315'
  mc_id = '545'
  nns = nn[[anchor[mc_id]]]
  Y = as.matrix(query_mat[,mc_id])[ids,]
  X = as.matrix(target_mat[,nns])[ids,]
  dim(X)
  dim(Y)
  p = ncol(X)
  betaHat <- Variable(p)
  objective  <- Minimize(sum((Y - X %*% betaHat)^2))
  #problem <- Problem(objective)
  eps = 1e-5
  problem <- Problem(objective, constraints = list(betaHat >= eps,sum(betaHat) == 1))
  result <- solve(problem)
  betas = result$getValue(betaHat)
  betas[betas<=eps] = 0
  sum(betas)
  sort(as.numeric(betas),decreasing = T)*100
  nnls.fit <- nnls::nnls(X, Y)$x
  sum(nnls.fit)
  
  library('nnls')
  
  if (requireNamespace("nnls", quietly = TRUE)) {
    nnls.fit <- nnls::nnls(X, Y)$x
  } else {
    nnls.fit <- rep(NA, p)
  }
  #
  
}
# somehow, it doesn't work with 10k genes 

# pipeline function 
MCPROJ_seeding_and_mixture_modeling = function(query_mat,target_mat,k = k,ncpu = ncpu,method = 'nnls'){
  require(parallel)
  # method - either "CVXR" or "nnls"
  message(sprintf('Mixture modeling of %s states using %s target states with %s nearest neighbors.',ncol(query_mat),ncol(target_mat),k))
  # input - matrices 
  # output - mixture projections
  # get anchors:
  anchor = get_nn(query_mat,target_mat,k = 1)
  # extend to the closest 10 metacells to anchor including itself
  nn = get_nn(target_mat,target_mat,k = 10)
  # Project metacells
  projections = mclapply(colnames(query_mat),FUN = function(mc_id){
    nns = nn[[anchor[mc_id]]]
    mixture_model_mc(Y=as.matrix(query_mat[,mc_id]),X = as.matrix(target_mat[,nns]),method = method)
  },mc.cores = ncpu)
  names(projections) = colnames(query_mat)
  sapply(projections,FUN = function(x) x$r2)
  return(projections)
}


################################################################################
# 2. Multiplicative gene correction
################################################################################
# If gene correction is needed (as determined by a user parameter), we compute for each gene the correlation between its query and projection expression over all metacells 
# and its overall expression fold factor 
# Genes with and a high degree of overall fold change  are corrected by setting. 
# Given corrected query data, we can reiterate the seeding and mixture phase. 
# We repeat this procedure at most three times (or less if no genes are corrected).
get_genes_to_correct = function(query_mat,proj_mat){
  # Q to P gene correlations 
  gene_corr = sapply(rownames(query_mat),FUN = function(x) cor(as.numeric(query_mat[x,]),as.numeric(proj_mat[x,])))
  summary(gene_corr)
  # expression fold factor
  ff = sapply(rownames(query_mat),FUN = function(x) sum(as.numeric(proj_mat[x,]))/sum(as.numeric(query_mat[x,])))
  # Correct genes with corr >= 0.8 and log(ff)>= log(1.15)
  genes_to_correct = names(which(gene_corr>=0.8&abs(log(ff))>=log(1.5)))
  return(list(gene_corr = gene_corr,fold_factor = ff,genes_to_correct = genes_to_correct))
}

correct_genes = function(query_mat,fold_factor,genes_to_correct){
  query_mat_corr = query_mat
  # correct the genes
  for(gene_id in genes_to_correct){
    query_mat_corr[gene_id,] = query_mat[gene_id,]*fold_factor[gene_id]  
  }
  return(query_mat_corr)
}

MCPROJ_multiplicative_gene_correction = function(query_mat,target_mat,proj_mat,k_nn = 10,max_iter=3,ncpu = 1,min_genes_to_correct = 1){
  counter = 1
  message('Getting genes to correct ...')
  genes_to_correct = get_genes_to_correct(query_mat,proj_mat)
  message(sprintf('it %s: %s genes to correct.',counter,length(genes_to_correct$genes_to_correct)))
  if(length(genes_to_correct$genes_to_correct)>=min_genes_to_correct){
    query_mat = correct_genes(query_mat,fold_factor = genes_to_correct$fold_factor,genes_to_correct$genes_to_correct)
    message(sprintf('it %s: Re-computing projections with corrected matrix.',counter))
    projections = MCPROJ_seeding_and_mixture_modeling(query_mat,target_mat,k = k,ncpu = ncpu)
    proj_mat = get_proj_matrix(projections)
    while(counter<=max_iter & length(genes_to_correct$genes_to_correct)>min_genes_to_correct){
      genes_to_correct = get_genes_to_correct(query_mat,proj_mat)
      message(sprintf('it %s: %s genes to correct.',counter,length(genes_to_correct$genes_to_correct)))
      query_mat = correct_genes(query_mat,fold_factor = genes_to_correct$fold_factor,genes_to_correct$genes_to_correct)
      # Re-compute the projections 
      message(sprintf('it %s: Re-computing projections with corrected matrix.',counter))
      projections = MCPROJ_seeding_and_mixture_modeling(query_mat,target_mat,k = k_nn,ncpu = ncpu)
      proj_mat = get_proj_matrix(projections)
      counter = counter + 1
    }
  }else{
    message('no genes to correct - skipping...')
  }
  return(list(query_mat = query_mat,proj_mat = proj_mat,projections = projections))
}



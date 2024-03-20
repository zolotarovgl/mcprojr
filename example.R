# Replicate MCProj algorithm 
library(metacell)
library(parallel)
ncpu = 5 # number of parallel threads 
k = 10 # number of target metacells selected to represent the queries
out_fn = 'results_mcproj/'
dir.create(out_fn)
source('mcprojr.R')
# TODOs:
# Solve strict non-negativity
################################################################################
# Data import 
################################################################################
if(F){
  setwd('~/Documents/projects/mlei_development/1110_joint_adult/')
  metacell::scdb_init("data_cdna/scdb",force_reinit=TRUE)
  run_name = "joint_Mlei02040506_CTfilt" 
  mc = metacell::scdb_mc(id = 'mc_filt')
  mat = metacell::scdb_mat(id = run_name)
  mlei2018_mc_fp = read.table('data_comparison/Mlei2018_v2/Mlei2018v2_metacell_gene_FC')
  colnames(target_mc_fp) = gsub('X','',colnames(target_mc_fp))
  metacell::scdb_init("../Mlei01/data_cdna/scdb",force_reinit=TRUE)
  mc_mlei01 = metacell::scdb_mc(id = '220523_Mlei01_10x_CTfilt')
  mlei01_mc = mc_mlei01
  mlei01_mat = metacell::scdb_mat(id = '220523_Mlei01_10x_CTfilt')
  mlei01_mc_anno = read.table('results_cydipid/metacell_annotation_transferred_mlei01.tsv',sep = '\t',header = T)
  mlei01_mc2type = setNames(mlei01_mc_anno$cell_type,mlei01_mc_anno$metacell)
}


###################################################################
mc_target = mlei01_mc
mc_query = mc
# This is only to filter the genes for faster executions 
fp_thr = 3
var_genes_query = names(apply(mc_query@mc_fp,1,max)[apply(mc_query@mc_fp,1,max)>=fp_thr])
var_genes_target = names(apply(mc_target@mc_fp,1,max)[apply(mc_target@mc_fp,1,max)>=fp_thr])

# Select noisy genes?
# E_gc input ###################################################################
# they are never negative 
eps = 1e-5
query_mat_input = log2(mc_query@e_gc + eps)
target_mat_input = log2(mc_target@e_gc + eps)
query_mat_input = query_mat_input[var_genes_query,]
target_mat_input = target_mat_input[var_genes_target,]
#query_mat_input = query_mat_input[apply(query_mat_input,1,var)>0.1,]
#target_mat_input = target_mat_input[apply(target_mat_input,1,var)>0.1,]
################################################################################
# 0. Compare genes
# Output the common set of genes 
################################################################################
query_genes = unique(rownames(query_mat_input))
target_genes = unique(rownames(target_mat_input))
int_genes = intersect(query_genes,target_genes)
message(sprintf('%s genes in query',length(query_genes)))
message(sprintf('%s genes in target',length(target_genes)))
message(sprintf('%s genes in common.',length(int_genes)))

genes = int_genes
message(sprintf('%s genes retained for downstream analysis.',length(genes)))
query_mat = query_mat_input[genes,]
target_mat = target_mat_input[genes,]
# Run initial set of projections:
projections = MCPROJ_seeding_and_mixture_modeling(query_mat,target_mat,k = 30,ncpu = ncpu,method = 'nnls')
# get projection matrix 
proj_mat = get_proj_matrix(projections)
# save initial round of projections 
#saveRDS(projections,sprintf('%s/projections.RData',out_fn))
#write.table(proj_mat,sprintf('%s/projection_matrix.tab',out_fn))
#projections = readRDS('results_mcproj/projections.RData')
#proj_mat = read.table('results_mcproj/projection_matrix.tab')
#colnames(proj_mat) = gsub('X','',colnames(proj_mat)) 
hist(sapply(projections,FUN = function(x) x$r2),breaks = 100,main = 'R2 distribution')
plot(sapply(projections,FUN = function(x) max(x$betas)),sapply(projections,FUN = function(x) x$r2),ylab = 'R^2', xlab = 'Max mc weight',col = mc@colors,pch = 16)

################################################################################
# 2. Multiplicative gene correction
################################################################################
old_r2 = sapply(projections,FUN = function(x) x$r2)
MGC_out = MCPROJ_multiplicative_gene_correction(query_mat,target_mat,proj_mat, k_nn = k,max_iter = 3,ncpu = 5,min_genes_to_correct = 1)
query_mat = MGC_out$query_mat
proj_mat = MGC_out$proj_mat
projections = MGC_out$projections
# save 
saveRDS(projections,'results_mcproj/projections.RData')
write.table(proj_mat,'results_mcproj/projection_matrix.tab')
projections = readRDS('results_mcproj/projections.RData')
proj_mat = read.table('results_mcproj/projection_matrix.tab')
colnames(proj_mat) = gsub('X','',colnames(proj_mat)) 

#plot(density(sapply(projections,FUN = function(x) x$r2)))
new_r2 = sapply(projections,FUN = function(x) x$r2)
plot(old_r2,new_r2,xlab = 'R^2 before gene correction',ylab = 'R^2 after gene correction')
abline(a = 0,b = 1)
hist(sapply(projections,FUN = function(x) x$r2),breaks = 100)
################################################################################
# 3. Range-based filtering
################################################################################
# For each projected gene, we compute the range of its (possibly corrected) query expressions and similarly its projected range. 
# If the shared expression range is less than 50% of the query range, then we filter out the gene as being too different between the query and the atlases
# this is not an expression range - this is gene enrichment range. 
dim(query_mat)
dim(target_mat)
dim(proj_mat)
gene_range_fraction = sapply(rownames(proj_mat),FUN = function(x){
  query_exp = as.numeric(query_mat[x,])
  proj_exp = as.numeric(proj_mat[x,])
  query_exp = query_exp[query_exp!=0]
  proj_exp = proj_exp[proj_exp!=0]
  if(!all(proj_exp == 0)){
    query_range = c(quantile(query_exp,.02),quantile(query_exp,.98))
    proj_range = c(quantile(proj_exp,.02),quantile(proj_exp,.98))
    shared_range = c(max(query_range[1],proj_range[1]),min(query_range[2],proj_range[2]))
    return(as.numeric(diff(shared_range)/diff(query_range)))
  }else{
    # gene isn't predicted to be expressed in the projection
    return(NA)
  }
})
message(sprintf('%s (%s %%) NA values ',sum(is.nan(gene_range_fraction)),round(sum(is.nan(gene_range_fraction))/length(gene_range_fraction)*100)))
gene_range_fraction[is.nan(gene_range_fraction)] = 0
gene_range_fraction[gene_range_fraction<0] = 0
hist(gene_range_fraction,breaks = 100)
# what if there are no expression in the projection? 
low_range_fraction = names(gene_range_fraction[gene_range_fraction<=.5])
message(sprintf('%s (%s%%) genes removed due to low range fraction.',length(low_range_fraction),round(length(low_range_fraction)/length(rownames(proj_mat))*100)))
# Remove  these genes from the comparisons 
# too many genes with low range fraction after multiple corrections 
################################################################################
# 4. Filtering genes within atlas types
################################################################################
eps = 1e-5
skew_mat = log2((proj_mat + eps)/(query_mat + eps))
# skewed genes per metacell
skewed_genes =  names(rowSums(abs(skew_mat)>=3)[rowSums(abs(skew_mat)>=3)>0])
message(sprintf('%s skewed genes (log2 fc >=3).',length(skewed_genes)))
################################################################################
# Identify type-skewed genes 
# T.B.A 
################################################################################
# metacell can be assigned to the type if it doesn't contain a lot of skewed genes 

################################################################################
# 5. Composite projection
# Model unassigned metacells as more flexible combinations of target states
################################################################################

################################################################################
# 6. Projection QC
################################################################################

################################################################################
# Explore
################################################################################

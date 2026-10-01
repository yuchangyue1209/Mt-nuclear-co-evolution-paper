#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(readxl);library(data.table)})
a <- commandArgs(trailingOnly=TRUE)
if(length(a)<3L) stop('Usage: Rscript test_watershed_mtcluster.R TABLE_S1_XLSX MT_CLUSTER_TSV OUTPUT_DIR [N_PERM]')
out <- a[3]; dir.create(out,recursive=TRUE,showWarnings=FALSE)
B <- if(length(a)>=4) as.integer(a[4]) else 9999L
if(is.na(B)||B<99L)stop('Use at least 99 permutations')
set.seed(20260930)
# One row per population, not one row per fish or SNP.
sheets <- excel_sheets(a[1])
meta <- NULL
for(sh in sheets) {
 z <- as.data.table(read_excel(a[1],sheet=sh))
 setnames(z,trimws(names(z)))
 if(all(c('Population','Region','Habitat','Watershed') %in% names(z))) {meta<-z;break}
}
if(is.null(meta))stop('No sheet has Population, Region, Habitat, Watershed columns')
meta <- meta[,.(pop=toupper(trimws(as.character(Population))),
 region=trimws(as.character(Region)),habitat=trimws(as.character(Habitat)),
 watershed=trimws(as.character(Watershed)))]
meta[,region:=fifelse(tolower(region)%in%c('alaska','ak'),'AK',fifelse(tolower(region)%in%c('bc','british columbia'),'BC',region))]
meta <- meta[!is.na(pop)&nzchar(pop)]
if(anyDuplicated(meta$pop))stop('Duplicate population IDs in Table S1')
cl <- fread(a[2], sep='\t', header=TRUE, fill=TRUE, blank.lines.skip=TRUE);setnames(cl,trimws(names(cl)))
if(!'pop'%in%names(cl))stop('Cluster table needs a pop column')
cluster_col <- intersect(c('mtCluster','mt_lineage','mt_cluster'),names(cl))
if(!length(cluster_col))stop('Cluster table needs mtCluster or mt_lineage')
cl <- cl[,.(pop=toupper(trimws(as.character(pop))),mt_cluster=trimws(as.character(get(cluster_col[1]))))]
cl <- unique(cl[!is.na(pop)&nzchar(pop)&!is.na(mt_cluster)&nzchar(mt_cluster)])
if(anyDuplicated(cl$pop))stop('Conflicting cluster assignments for a population')
merged <- merge(meta,cl,by='pop',all.x=TRUE,sort=FALSE)
# Primary analysis: established and recently colonized freshwater populations.
merged[,excluded_reason:=fifelse(grepl('marine',habitat,ignore.case=TRUE),'Marine reference',
 fifelse(is.na(mt_cluster)|!nzchar(mt_cluster),'Missing mt cluster',
 fifelse(is.na(watershed)|!nzchar(watershed)|tolower(watershed)%in%c('na','n/a','unknown','-'),'Missing watershed',
 fifelse(is.na(region)|!nzchar(region),'Missing region',''))))]
fwrite(merged,file.path(out,'population_input_audit.tsv'),sep='\t')
d <- merged[excluded_reason=='']
if(nrow(d)<4||uniqueN(d$mt_cluster)<2)stop('Too few complete freshwater populations or clusters')
fwrite(d,file.path(out,'populations_analyzed.tsv'),sep='\t')
counts <- d[,.(n_populations=.N,populations=paste(sort(pop),collapse=', ')),by=.(region,watershed,mt_cluster)]
fwrite(counts,file.path(out,'watershed_mtcluster_counts.tsv'),sep='\t')
wide <- dcast(counts,region+watershed~mt_cluster,value.var='n_populations',fill=0)
fwrite(wide,file.path(out,'watershed_mtcluster_crosstab.tsv'),sep='\t')
ws <- d[,.(n_populations=.N,n_clusters=uniqueN(mt_cluster),
 clusters=paste(sort(unique(mt_cluster)),collapse=', ')),by=.(region,watershed)]
fwrite(ws,file.path(out,'watershed_summary.tsv'),sep='\t')

# Test watershed-cluster association conditional on region.
region_indices <- split(seq_len(nrow(d)),d$region)
watershed_groups <- interaction(d$region,d$watershed,drop=TRUE)
association_stat <- function(y) {
 sum(vapply(region_indices,function(ii) {
  t <- table(d$watershed[ii],y[ii])
  if(nrow(t)<2||ncol(t)<2)return(0)
  expected <- outer(rowSums(t),colSums(t))/sum(t)
  sum((t-expected)^2/expected)
 },numeric(1)))
}
# LOOCV: predict only where another population shares region+watershed.
# Majority ties receive fractional credit; no arbitrary tie-breaking.
eligible <- which(ave(seq_len(nrow(d)),watershed_groups,FUN=length)>1)
majority <- function(v) {
 t<-table(v);names(t)[t==max(t)]
}
cv <- function(y,detail=FALSE) {
 if(!length(eligible))return(if(detail)data.table() else c(watershed=NA_real_,region=NA_real_,gain=NA_real_))
 rows <- lapply(eligible,function(i) {
  same_ws<-which(watershed_groups==watershed_groups[i]);same_ws<-setdiff(same_ws,i)
  same_reg<-setdiff(which(d$region==d$region[i]),i)
  w<-majority(y[same_ws]);r<-majority(y[same_reg])
  data.table(pop=d$pop[i],region=d$region[i],watershed=d$watershed[i],observed=y[i],
   predicted_watershed=paste(w,collapse=' | '),predicted_region=paste(r,collapse=' | '),
   watershed_score=as.numeric(y[i]%in%w)/length(w),region_score=as.numeric(y[i]%in%r)/length(r),
   n_training_same_watershed=length(same_ws))
 })
 z<-rbindlist(rows)
 if(detail)return(z)
 c(watershed=mean(z$watershed_score),region=mean(z$region_score),gain=mean(z$watershed_score-z$region_score))
}
y<-d$mt_cluster;obs_assoc<-association_stat(y);obs_cv<-cv(y)
pred<-cv(y,TRUE)
fwrite(pred,file.path(out,'leave_one_population_out_predictions.tsv'),sep='\t')
unassessable<-d[setdiff(seq_len(nrow(d)),eligible)]
fwrite(unassessable,file.path(out,'populations_without_watershed_training_peers.tsv'),sep='\t')
perm_assoc<-numeric(B);perm_gain<-rep(NA_real_,B)
for(b in seq_len(B)) {
 yp<-y
 for(ii in region_indices)yp[ii]<-sample(y[ii],length(ii),replace=FALSE)
 perm_assoc[b]<-association_stat(yp)
 perm_gain[b]<-cv(yp)['gain']
 if(b%%1000L==0L)cat('Permutations:',b,'/',B,'\n')
}
p_assoc<-(1+sum(perm_assoc>=obs_assoc-1e-12))/(B+1)
p_gain<-if(is.finite(obs_cv['gain']))(1+sum(perm_gain>=obs_cv['gain']-1e-12))/(B+1) else NA_real_
summary<-data.table(metric=c('freshwater_populations_analyzed','regions','region_watershed_groups',
 'watersheds_containing_multiple_clusters','cv_evaluable_populations','cv_coverage',
 'watershed_cv_accuracy_tie_adjusted','region_cv_accuracy_tie_adjusted','cv_accuracy_gain',
 'association_statistic_within_region','association_permutation_P','cv_gain_permutation_P','permutations'),
 value=c(nrow(d),uniqueN(d$region),nrow(ws),sum(ws$n_clusters>1),length(eligible),length(eligible)/nrow(d),
 obs_cv['watershed'],obs_cv['region'],obs_cv['gain'],obs_assoc,p_assoc,p_gain,B))
fwrite(summary,file.path(out,'watershed_prediction_summary.tsv'),sep='\t')
print(summary)
writeLines(c(
 'One population is one observation. Marine references are excluded; recent freshwater populations are included when cluster assignments are available.',
 'Missing metadata and cluster assignments are recorded in population_input_audit.tsv.',
 'Association and prediction-gain permutation tests shuffle mt-cluster labels within AK/BC, preserving regional cluster counts.',
 'Cross-validation holds out one population and predicts its cluster from other populations in the same region and watershed.',
 'Singleton watersheds are not scored. Region-only accuracy uses the same held-out populations.',
 'Majority ties receive fractional credit. Accuracy is therefore an expected accuracy, not an arbitrarily tie-broken classification rate.',
 'These predictions apply to additional sampled populations from represented watersheds, not to previously unsampled watersheds.',
 'Table S1 watershed labels include broad geographic areas; significance does not establish a causal drainage effect.'
),file.path(out,'analysis_notes.txt'))
cat('[OK] Results saved to:',out,'\n')

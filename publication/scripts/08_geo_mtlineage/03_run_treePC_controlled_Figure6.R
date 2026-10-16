# Historical treePC-adjusted alternative; current entry point: 03_run_Figure6_associations.R
#!/usr/bin/env Rscript
# A: geographic SNP associations. B: cluster association controlling treePC1/2.
# C: within-gene multi-SNP profiles (one point per gene).
# Separate output: one SNP per gene, population profiles across genes.
# Usage: Rscript script.R LD_PRUNED.tsv.gz METADATA.tsv NUCLEAR_TREE_OR_PCS GENES.tsv OUTDIR [N_PERM]
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(ggrepel)})
a<-commandArgs(trailingOnly=TRUE)
if(length(a)<5L)stop('Usage: Rscript 03_run_treePC_controlled_Figure6.R LD_PRUNED.tsv.gz METADATA.tsv NUCLEAR_TREE_OR_PCS GENES.tsv OUTDIR [N_PERM]')
out<-a[5];dir.create(out,recursive=TRUE,showWarnings=FALSE)
B<-if(length(a)>5L)as.integer(a[6])else 9999L
if(is.na(B)||B<99L)stop('N_PERM must be >=99')
save_tsv<-function(x,name)fwrite(x,file.path(out,name),sep='\t',na='NA')
rename_first<-function(d,target,aliases){
 if(!target%in%names(d)){
  hit<-aliases[aliases%in%names(d)]
  if(length(hit))setnames(d,hit[1],target)
 };d
}
norm_pop<-function(x)toupper(trimws(as.character(x)))
meta<-fread(a[2]);meta<-rename_first(meta,'pop',c('Population','population','Pop'))
meta<-rename_first(meta,'latitude',c('Latitude','lat'))
meta<-rename_first(meta,'longitude',c('Longitude','lon','long'))
meta<-rename_first(meta,'mt_lineage',c('mtCluster','mitochondrial_cluster','cluster'))
meta<-rename_first(meta,'region',c('Region'))
required<-c('pop','region','latitude','longitude','mt_lineage')
if(!all(required%in%names(meta)))stop('Metadata needs: ',paste(required,collapse=', '))
meta[,pop:=norm_pop(pop)];meta[,region:=norm_pop(region)]
if(anyDuplicated(meta$pop))stop('Duplicate populations in metadata')
longterm<-c('FG','LG','SL','SR','TL','WB','WK','WT','SWA','THE','JOE','BEA','MUC','PYE','AMO','ROB','BOOT','ECHO','GOS','LAW')
meta<-meta[pop%in%setdiff(longterm,'AMO')]
if(nrow(meta)!=19L)stop('Metadata must contain all 19 long-term freshwater populations excluding AMO')
meta[,latitude:=as.numeric(latitude)];meta[,longitude:=as.numeric(longitude)]
if(any(!is.finite(meta$latitude))||any(!is.finite(meta$longitude))||anyNA(meta$mt_lineage))stop('Missing geography or cluster metadata')
meta[,region:=factor(region,levels=c('AK','BC'))]
meta[,mt_lineage:=factor(mt_lineage)]
if(anyNA(meta$region))stop('Region must be AK or BC')

# Third input may be an existing PC table or the original nuclear population tree.
first<-trimws(readLines(a[3],n=1,warn=FALSE))
if(startsWith(first,'(')){
 if(!requireNamespace('ape',quietly=TRUE))stop('Package ape is needed for tree input')
 tr<-ape::read.tree(a[3]);tr$tip.label<-norm_pop(tr$tip.label)
 study<-c(longterm,'PACH','FRED','SC','LB','CH','RS','SAY')
 if(!all(study%in%tr$tip.label))stop('Nuclear tree must contain the 27 study populations; alternatively supply a PC table')
 tr<-ape::keep.tip(tr,study)
 pc<-cmdscale(ape::cophenetic.phylo(tr),k=2,eig=TRUE)
 pcs<-data.table(pop=rownames(pc$points),treePC1=pc$points[,1],treePC2=pc$points[,2])
 save_tsv(data.table(axis=seq_along(pc$eig),eigenvalue=pc$eig,
   fraction_positive=ifelse(pc$eig>0,pc$eig/sum(pc$eig[pc$eig>0]),NA_real_)), 'nuclear_PCoA_eigenvalues.tsv')
}else{
 pcs<-fread(a[3]);pcs<-rename_first(pcs,'pop',c('Population','population','Pop'))
 if(!all(c('pop','treePC1','treePC2')%in%names(pcs)))stop('PC table needs pop, treePC1, treePC2')
 pcs[,pop:=norm_pop(pop)]
}
if(anyDuplicated(pcs$pop))stop('Duplicate populations in PC input')
old_pc<-intersect(c('treePC1','treePC2'),names(meta));if(length(old_pc))meta[,(old_pc):=NULL]
meta<-merge(meta,pcs[,.(pop,treePC1,treePC2)],by='pop',all.x=TRUE,sort=FALSE)
if(any(!is.finite(meta$treePC1))||any(!is.finite(meta$treePC2)))stop('Missing nuclear PCs')
meta[,latitude_z:=as.numeric(scale(latitude))];meta[,longitude_z:=as.numeric(scale(longitude))]
save_tsv(meta,'populations_and_covariates.tsv')

genes0<-fread(a[4],header=FALSE)
genes<-unique(tolower(trimws(as.character(genes0[[1]]))))
genes<-setdiff(genes,c('gene','symbol','gene_symbol'))
if(length(genes)!=72L)stop('GENES input must list exactly 72 gene symbols in its first column')
d<-fread(a[1]);d<-rename_first(d,'pop',c('Pop','population','Population'))
d<-rename_first(d,'gene',c('Gene','symbol'));d<-rename_first(d,'chr',c('CHR','chromosome'))
d<-rename_first(d,'pos',c('POS','position'));d<-rename_first(d,'af',c('AF','focal_AF'))
if(!all(c('pop','gene','chr','pos','depth','af')%in%names(d)))stop('LD-pruned input needs pop, gene, chr, pos, depth, af; also deltaAF or marine_af for panel C')
d[,pop:=norm_pop(pop)];d[,gene:=tolower(as.character(gene))]
d<-d[pop%in%meta$pop & gene%in%genes & is.finite(depth) & depth>=20 & is.finite(af) & af>=0 & af<=1]
d[,snp_id:=paste(gene,chr,pos,sep='__')]
if(anyDuplicated(d[,.(snp_id,pop)]))stop('Duplicate SNP/population records; resolve these before fitting')
if('focal_count'%in%names(d)){
 d[,focal_count:=as.numeric(focal_count)]
}else d[,focal_count:=round(af*depth)]
if(any(!is.finite(d$focal_count))||any(d$focal_count<0|d$focal_count>d$depth))stop('Invalid focal counts')
d[,other_count:=depth-focal_count]
if(!'deltaAF'%in%names(d)){
 if('delta_af'%in%names(d))setnames(d,'delta_af','deltaAF')
 else if('marine_af'%in%names(d))d[,deltaAF:=af-marine_af]
 else stop('Need deltaAF, delta_af, or marine_af for signed profiles')
}
# Use metadata as the single authoritative source of model covariates.
drop<-intersect(setdiff(names(meta),'pop'),names(d));if(length(drop))d[,(drop):=NULL]
d<-merge(d,meta,by='pop',all.x=TRUE,sort=FALSE)
safe_test<-function(reduced,full){
 if(full$rank<=reduced$rank||df.residual(full)<1)return(NA_real_)
 z<-tryCatch(anova(reduced,full,test='F'),error=function(e)NULL)
 if(is.null(z))return(NA_real_)
 as.numeric(z[2,'Pr(>F)'])
}
fit_snp<-function(z,kind){
 z<-droplevels(as.data.frame(z));n<-nrow(z)
 ans<-list(n_pops=n,p=NA_real_,effect=NA_real_,status='insufficient_data')
 if(n<10L||length(unique(z$mt_lineage))<2L)return(ans)
 result<-tryCatch({
  if(kind=='geo'){
   r<-glm(cbind(focal_count,other_count)~1,family=quasibinomial(),data=z)
   f<-glm(cbind(focal_count,other_count)~latitude_z+longitude_z,family=quasibinomial(),data=z)
   pred<-predict(f,type='response');effect<-diff(range(pred))
  }else{
   r<-glm(cbind(focal_count,other_count)~treePC1+treePC2,family=quasibinomial(),data=z)
   f<-glm(cbind(focal_count,other_count)~treePC1+treePC2+mt_lineage,family=quasibinomial(),data=z)
   if(anyNA(coef(f)))stop('rank deficient')
   # Adjusted cluster frequencies: predict every cluster over the SAME observed
   # PC distribution, with pooled read-depth weights, then compare averages.
   marginal<-vapply(levels(z$mt_lineage),function(cl){
    nd<-z;nd$mt_lineage<-factor(cl,levels=levels(z$mt_lineage))
    weighted.mean(predict(f,newdata=nd,type='response'),w=z$depth)
   },numeric(1))
   effect<-diff(range(marginal))
  }
  p<-safe_test(r,f)
  list(n_pops=n,p=p,effect=effect,status=if(is.finite(p)&&f$converged)'ok'else'nonestimable')
 },error=function(e)list(n_pops=n,p=NA_real_,effect=NA_real_,status=conditionMessage(e)))
 result
}
cat('Fitting SNP models on',nrow(meta),'populations and',uniqueN(d$snp_id),'SNPs...\n')
geo<-d[,fit_snp(.SD,'geo'),by=.(gene,chr,pos,snp_id)]
mt<-d[,fit_snp(.SD,'mt'),by=.(gene,chr,pos,snp_id)]
geo[,q:=p.adjust(p,'BH')];mt[,q:=p.adjust(p,'BH')]
for(z in list(geo,mt))z[,candidate:=is.finite(q)&q<.05&is.finite(effect)&effect>=.20]
save_tsv(geo,'geography_SNP_models.tsv');save_tsv(mt,'mitochondrial_cluster_treePC_SNP_models.tsv')
gene_summary<-function(x){
 ans<-x[,.(n_tested=sum(is.finite(p)),n_FDR=sum(q<.05,na.rm=TRUE),n_candidate=sum(candidate),
   min_q=if(any(is.finite(q)))min(q,na.rm=TRUE)else NA_real_),by=gene]
 merge(data.table(gene=genes),ans,by='gene',all.x=TRUE)
}
save_tsv(gene_summary(geo),'geography_gene_summary.tsv')
save_tsv(gene_summary(mt),'mitochondrial_cluster_treePC_gene_summary.tsv')

# A single collection of region-stratified label permutations is reused across
# genes. This preserves the cluster composition of each region.
setorder(meta,pop);set.seed(10615)
labels<-as.character(meta$mt_lineage)
perm_labels<-replicate(B,{
 z<-labels
 for(reg in levels(meta$region)){
  ix<-which(meta$region==reg);z[ix]<-sample(labels[ix],length(ix),replace=FALSE)
 };z
})
permutation_distance<-function(mat){
 if(ncol(mat)<3L)return(list(diff=NA_real_,p=NA_real_,status='fewer_than_3_variable_SNPs'))
 if(any(apply(mat,1,sd)==0))return(list(diff=NA_real_,p=NA_real_,status='constant_population_profile'))
 cm<-cor(t(mat),method='pearson');dist<-1-cm
 ij<-which(upper.tri(dist),arr.ind=TRUE);v<-dist[ij]
 ii<-match(rownames(mat),meta$pop);cl<-labels[ii]
 same<-cl[ij[,1]]==cl[ij[,2]]
 if(!any(same)||all(same))return(list(diff=NA_real_,p=NA_real_,status='missing_pair_class'))
 observed<-mean(v[same])-mean(v[!same])
 null<-vapply(seq_len(B),function(b){
  lab<-perm_labels[ii,b];s<-lab[ij[,1]]==lab[ij[,2]]
  if(!any(s)||all(s))return(NA_real_)
  mean(v[s])-mean(v[!s])
 },numeric(1))
 valid<-is.finite(null)
 list(diff=observed,p=(1+sum(null[valid]<=observed+1e-12))/(1+sum(valid)),
  status='ok',mean_within=mean(v[same]),mean_between=mean(v[!same]),
  permutations_used=sum(valid),pairs=data.table(pop1=rownames(mat)[ij[,1]],pop2=rownames(mat)[ij[,2]],
   distance=v,same_cluster=same))
}
make_matrix<-function(z){
 wide<-dcast(z,pop~snp_id,value.var='deltaAF')
 if(nrow(wide)!=19L)return(NULL)
 mat<-as.matrix(wide[,-1]);rownames(mat)<-wide$pop
 keep<-colSums(is.finite(mat))==19L
 mat<-mat[,keep,drop=FALSE]
 if(ncol(mat))mat<-mat[,apply(mat,2,sd)>0,drop=FALSE]
 mat
}
profile_rows<-list()
for(g in genes){
 mat<-make_matrix(d[gene==g & is.finite(deltaAF)])
 if(is.null(mat)){profile_rows[[g]]<-data.table(gene=g,n_snps=0,diff_mean=NA_real_,p_perm_one_sided=NA_real_,status='incomplete_population_set');next}
 result<-permutation_distance(mat)
 profile_rows[[g]]<-data.table(gene=g,n_snps=ncol(mat),diff_mean=result$diff,
  p_perm_one_sided=result$p,status=result$status)
}
profiles<-rbindlist(profile_rows);profiles[,q_perm_one_sided:=p.adjust(p_perm_one_sided,'BH')]
save_tsv(profiles,'Figure6C_within_gene_profiles.tsv')

# Additional requested analysis: one shared SNP per gene, selected without
# cluster labels. This is a POPULATION-level profile across different genes.
eligible<-d[,.(n_pops=.N,mean_abs_deltaAF=mean(abs(deltaAF)),mean_depth=mean(depth)),by=.(gene,snp_id)]
eligible<-eligible[n_pops==19L & is.finite(mean_abs_deltaAF)]
setorder(eligible,gene,-mean_abs_deltaAF,-mean_depth,snp_id)
selected<-eligible[,.SD[1],by=gene]
save_tsv(selected,'cross_gene_one_SNP_per_gene_selection.tsv')
cross_data<-d[selected,on=.(gene,snp_id),nomatch=0]
cross_mat<-make_matrix(cross_data)
if(!is.null(cross_mat)){
 cr<-permutation_distance(cross_mat)
 save_tsv(data.table(n_genes=ncol(cross_mat),diff_mean=cr$diff,p_perm_one_sided=cr$p,status=cr$status),
  'cross_gene_one_SNP_profile_test.tsv')
 if(!is.null(cr$pairs))save_tsv(cr$pairs,'cross_gene_one_SNP_population_pairs.tsv')
}

theme_fig<-theme_classic(base_size=22)+theme(axis.text=element_text(colour='black',size=18),
 axis.title=element_text(size=22),plot.title=element_text(size=26,face='bold'),
 legend.position='bottom',legend.title=element_blank(),legend.text=element_text(size=16),
 plot.margin=margin(12,16,12,12))
cols<-c('Neither'='grey70','FDR only'='#4C78A8','Effect only'='#E45756','Candidate'='#C77CFF')
volcano<-function(z,title,xlab){
 z<-copy(z);z[,group:=fifelse(candidate,'Candidate',fifelse(is.finite(q)&q<.05,'FDR only',
   fifelse(is.finite(effect)&effect>=.2,'Effect only','Neither')))]
 z<-z[is.finite(q)&is.finite(effect)];z[,y:=-log10(pmax(q,.Machine$double.xmin))]
 lab<-head(z[candidate==TRUE][order(q,-effect),.SD[1],by=gene],12)
 ggplot(z,aes(effect,y))+geom_point(aes(colour=group),size=2,alpha=.75)+
 geom_hline(yintercept=-log10(.05),linetype=2,linewidth=.4)+geom_vline(xintercept=.2,linetype=2,linewidth=.4)+
 geom_text_repel(data=lab,aes(label=gene),fontface='italic',colour='black',size=5,
  max.overlaps=Inf,seed=10615,box.padding=.5,max.time=5)+
 scale_colour_manual(values=cols,drop=FALSE)+labs(title=title,x=xlab,y=expression(-log[10]('FDR q-value')))+theme_fig
}
pA<-volcano(geo,'A  Geography','Fitted allele-frequency range\nacross latitude/longitude')
pB<-volcano(mt,'B  Mitochondrial cluster','Adjusted fitted allele-frequency difference\namong mitochondrial clusters')
pp<-profiles[is.finite(diff_mean)&is.finite(p_perm_one_sided)]
pp[,group:=fifelse(diff_mean<0 & p_perm_one_sided<.05,'Negative difference, P < 0.05',
 fifelse(diff_mean<0,'Negative difference','Nonnegative difference'))]
pp[,y:=-log10(p_perm_one_sided)]
labC<-head(pp[diff_mean<0][order(p_perm_one_sided,diff_mean)],12)
pC<-ggplot(pp,aes(diff_mean,y))+geom_point(aes(colour=group),size=3)+
 geom_vline(xintercept=0,linetype=2,linewidth=.4)+geom_hline(yintercept=-log10(.05),linetype=2,linewidth=.4)+
 geom_text_repel(data=labC,aes(label=gene),fontface='italic',colour='black',size=5,
  max.overlaps=Inf,seed=10615,box.padding=.5,max.time=5)+
 scale_colour_manual(values=c('Negative difference, P < 0.05'='#E45756','Negative difference'='#4C78A8','Nonnegative difference'='grey70'))+
 labs(title='C  Within-gene signed allele-frequency profiles',x='Mean within-cluster distance - between-cluster distance',
  y=expression(-log[10]('permutation P-value')))+theme_fig
save_plot<-function(p,stem,w,h){
 ggsave(file.path(out,paste0(stem,'.pdf')),p,width=w,height=h,device=cairo_pdf,bg='white')
 ggsave(file.path(out,paste0(stem,'.png')),p,width=w,height=h,dpi=220,bg='white')
}
save_plot(pA,'Figure6A_treePC_workflow',9,9);save_plot(pB,'Figure6B_treePC_adjusted',9,9)
save_plot(pC,'Figure6C_within_gene_profiles',13,9)
draw<-function(){
 grid::grid.newpage();grid::pushViewport(grid::viewport(layout=grid::grid.layout(2,2)))
 print(pA,vp=grid::viewport(layout.pos.row=1,layout.pos.col=1),newpage=FALSE)
 print(pB,vp=grid::viewport(layout.pos.row=1,layout.pos.col=2),newpage=FALSE)
 print(pC,vp=grid::viewport(layout.pos.row=2,layout.pos.col=1:2),newpage=FALSE);grid::popViewport()
}
cairo_pdf(file.path(out,'Figure6_ABC_treePC.pdf'),width=18,height=18);draw();dev.off()
png(file.path(out,'Figure6_ABC_treePC.png'),width=18,height=18,units='in',res=180,type='cairo');draw();dev.off()
counts<-rbindlist(list(geo[,.(analysis='geography',n_tested=sum(is.finite(p)),n_FDR=sum(q<.05,na.rm=TRUE),n_candidates=sum(candidate))],
 mt[,.(analysis='mitochondrial_cluster_treePC',n_tested=sum(is.finite(p)),n_FDR=sum(q<.05,na.rm=TRUE),n_candidates=sum(candidate))]))
save_tsv(counts,'model_counts.tsv');print(counts)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
writeLines(c('All panels use the same 19 long-term freshwater populations excluding AMO.',
 'A: latitude/longitude quasibinomial model versus intercept. B: treePC1+treePC2+cluster versus treePC1+treePC2.',
 'BH correction is separate for A and B across estimable SNP tests. Candidate threshold q<0.05 and effect>=0.20.',
 'B effect: cluster-standardized fitted frequencies averaged over the same observed PC distribution with read-depth weights.',
 'C: one point per gene, multiple LD-pruned SNPs per gene; complete data across 19 populations; constant SNP columns removed.',
 'C uses signed deltaAF profiles, not PC residuals. Cluster labels are permuted within region; nominal and BH values are exported.',
 'Cross-gene supplementary test selects one shared SNP per gene using mean absolute deltaAF, WITHOUT cluster-label selection.',
 'Cross-gene pairwise distances are dependent; inference uses whole-population label permutations, not tests treating pairs as independent.',
 'The input must already be LD-pruned; this script does not infer read-backed LD.',
 paste('Permutation replicates:',B), 'Existing manuscript counts cannot be reused after changing model or filtering.'),file.path(out,'analysis_notes.txt'))
cat('[OK] Figure6 A/B/C and one-SNP cross-gene profiles saved to ',out,'\n',sep='')

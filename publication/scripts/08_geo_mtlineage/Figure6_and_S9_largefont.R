#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(ggrepel)})
a <- commandArgs(trailingOnly=TRUE)
if (length(a) != 3) stop('Usage: Rscript Figure6_and_S9_largefont.R GLM_RESULTS_DIR PARALLELISM_RESULTS_TSV OUTPUT_DIR')
glm_dir <- a[1]; parallel_file <- a[2]; outdir <- a[3]
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
geo <- fread(file.path(glm_dir,'OXPHOS72_geography_LatLong_quasibinomial_GLM_results.tsv'))
mt <- fread(file.path(glm_dir,'OXPHOS72_mtlineage_allPops_quasibinomial_GLM_results.tsv'))
par <- fread(parallel_file)
cols <- c('Not candidate'='grey70','FDR only'='#4C78A8','Large effect only'='#E45756','Candidate'='#C77CFF')
font_theme <- theme_classic(base_size=26) + theme(
 axis.title=element_text(size=28,colour='black'), axis.text=element_text(size=23,colour='black'),
 axis.line=element_line(colour='black',linewidth=0.7),
 plot.title=element_text(size=30,face='bold',hjust=0),
 legend.title=element_blank(),legend.text=element_text(size=23),
 legend.position='bottom',legend.key.size=grid::unit(0.65,'cm'),
 strip.text=element_text(size=28,face='bold'),
 plot.margin=margin(14,18,14,14))
make_volcano <- function(d,which) {
 is_geo <- which=='geo'
 eff <- if(is_geo) 'geo_fitted_AF_range' else 'max_lineage_AF_diff'
 q <- if(is_geo) 'q_geo' else 'q_mt'
 neglog <- if(is_geo) 'neglog10_q_geo' else 'neglog10_q_mt'
 sig <- if(is_geo) 'significant_geo' else 'significant_mt'
 large <- if(is_geo) 'large_effect_geo' else 'large_effect_mt'
 cand <- if(is_geo) 'candidate_geo' else 'candidate_mt'
 d <- copy(d)
 d[, group := fifelse(get(cand),'Candidate',fifelse(get(sig),'FDR only',fifelse(get(large),'Large effect only','Not candidate')))]
 d[, group:=factor(group,levels=names(cols))]
 d[, y_plot:=pmin(get(neglog),35)]
 lab <- d[get(cand)==TRUE]
 lab <- lab[order(lab[[q]],-lab[[eff]])][seq_len(min(nrow(lab),15))]
 title <- if(is_geo) 'A  Geography' else 'B  Mitochondrial lineage'
 xlabel <- if(is_geo) 'Fitted allele-frequency range\nacross latitude/longitude' else 'Maximum allele-frequency difference\namong mitochondrial lineages'
 ggplot(d,aes(x=.data[[eff]],y=y_plot))+
  geom_point(aes(colour=group),size=2.3,alpha=0.75)+
  geom_hline(yintercept=-log10(0.05),linetype='dashed',linewidth=0.5)+
  geom_vline(xintercept=0.2,linetype='dashed',linewidth=0.5)+
  geom_text_repel(data=lab,aes(label=gene),size=8.0,fontface='italic',colour='black',
    nudge_y=if(is_geo) 0 else 3,
    box.padding=0.8,point.padding=0.45,segment.size=0.35,max.overlaps=Inf,
    force=3,max.time=10,max.iter=50000,seed=106)+
  scale_colour_manual(values=cols,drop=FALSE)+
  scale_y_continuous(expand=expansion(mult=c(0.02,0.15)))+
  labs(title=title,x=xlabel,y=expression(-log[10]('FDR q-value')))+font_theme+
  guides(colour=guide_legend(nrow=2,byrow=TRUE,override.aes=list(size=4,alpha=1)))
}
pA <- make_volcano(geo,'geo'); pB <- make_volcano(mt,'mt')
par[, group:=fifelse(diff_mean<0 & p_perm_one_sided<0.05,'Parallel p < 0.05',fifelse(diff_mean<0,'Parallel','No parallel'))]
par[, group:=factor(group,levels=c('No parallel','Parallel','Parallel p < 0.05'))]
par[, y_plot:=-log10(pmax(p_perm_one_sided,.Machine$double.xmin))]
pcols <- c('No parallel'='grey70','Parallel'='#4C78A8','Parallel p < 0.05'='#E45756')
make_parallel <- function(d,followup=FALSE) {
 lab <- d[diff_mean<0][order(p_perm_one_sided,diff_mean), head(.SD,if(followup)10L else 15L), by=scope]
 p <- ggplot(d,aes(diff_mean,y_plot))+
  geom_point(aes(colour=group),size=3,alpha=0.85)+
  geom_vline(xintercept=0,linetype='dashed',colour='grey45',linewidth=0.5)+
  geom_hline(yintercept=-log10(0.05),linetype='dashed',colour='grey45',linewidth=0.5)+
  geom_text_repel(data=lab,aes(label=gene),size=8.5,fontface='italic',colour='black',
    box.padding=0.8,point.padding=0.5,segment.size=0.35,max.overlaps=Inf,
    force=3,max.time=10,max.iter=50000,seed=106)+
  scale_colour_manual(values=pcols,drop=FALSE)+
  scale_y_continuous(expand=expansion(mult=c(0.03,0.17)))+
  labs(title=if(followup)'Region-specific follow-up' else 'C  All freshwater populations',
       x='Mean within-lineage distance - between-lineage distance',
       y=expression(-log[10]('permutation p-value')))+font_theme+
  guides(colour=guide_legend(nrow=1,override.aes=list(size=4,alpha=1)))
 if(followup) p<-p+facet_wrap(~scope,nrow=1,scales='free_x')
 p
}
pC <- make_parallel(par[scope=='all'])
pS <- make_parallel(par[scope %in% c('AK','BC')],TRUE)
save_panel <- function(p,name,w,h) {
 ggsave(file.path(outdir,paste0(name,'.pdf')),p,width=w,height=h,device=grDevices::cairo_pdf,bg='white')
 ggsave(file.path(outdir,paste0(name,'.png')),p,width=w,height=h,dpi=300,bg='white')
}
save_panel(pA,'Figure6A_largefont',9,10)
save_panel(pB,'Figure6B_largefont',9,10)
save_panel(pC,'Figure6C_largefont',13,10)
save_panel(pS,'FigureS9_largefont',15,11)
# Compose with base grid: no extra layout package required.
draw_combined <- function() {
 grid::grid.newpage()
 grid::pushViewport(grid::viewport(layout=grid::grid.layout(2,2,heights=grid::unit(c(1,1),'null'))))
 print(pA,vp=grid::viewport(layout.pos.row=1,layout.pos.col=1),newpage=FALSE)
 print(pB,vp=grid::viewport(layout.pos.row=1,layout.pos.col=2),newpage=FALSE)
 print(pC,vp=grid::viewport(layout.pos.row=2,layout.pos.col=1:2),newpage=FALSE)
 grid::popViewport()
}
grDevices::cairo_pdf(file.path(outdir,'Figure6_largefont.pdf'),width=18,height=21)
draw_combined(); grDevices::dev.off()
grDevices::png(file.path(outdir,'Figure6_largefont.png'),width=18,height=21,units='in',res=300,type='cairo')
draw_combined(); grDevices::dev.off()
cat('[OK] Main Figure 6, panels A/B/C, and regional Figure S9 saved to',outdir,'\n')

#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(ggsignif);library(patchwork)})

root <- "/path/to/workspace/Kuster2026_genomewide_reanalysis/genomewide_pinpis"
dat <- fread(file.path(root,"results/genomewide_nuclear_mt13_pinpis_combined.tsv"),na.strings=c("NA","NaN",""))
tests <- fread(file.path(root,"results/pinpis_group_tests.tsv"),na.strings=c("NA","NaN",""))
outdir <- file.path(root,"figures/Figure2_pinpis_large_fonts");dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

BORDER_LWD<-0.7; AXIS_TITLE_SIZE<-30; AXIS_TEXT_SIZE<-21; N_SIZE<-6.5; STAR_SIZE<-8
ROW_TAG_SIZE<-32; POINT_SIZE<-1.9; POINT_ALPHA<-0.62; set.seed(20260810)
pal_overall<-c(mtOXPHOS="#7FA4B2",nuOXPHOS="#BC8984",`assembly factor`="#9EAD8D")
pal_rp<-c(`cyto-RP`="#7795A8",`Nmt-RP`="#C28C73")
pal_ars<-c(`cyto-ARS`="#7589A9",`Nmt-ARS`="#B77D9A")
pal_core<-c(mtOXPHOS="#7FA4B2",`nu core`="#BC8984",`nu noncore`="#B6A174")
pal_kuster<-c(`direct_n-mt`="#BC8984",`indirect_n-mt`="#948C9D",`non-n-mt`="#AAA8A1")

metrics<-c(piN="piN_mean",piS="piS_mean",piN_piS="piN_piS_ratio_of_means")
metric_lab<-list(piN=expression(pi[N]),piS=expression(pi[S]),piN_piS=expression(pi[N]/pi[S]))
background<-c(piN=median(dat[gene_source=="nuclear"]$piN_mean,na.rm=TRUE),
              piS=median(dat[gene_source=="nuclear"]$piS_mean,na.rm=TRUE),
              piN_piS=median(dat[gene_source=="nuclear"]$piN_piS_ratio_of_means,na.rm=TRUE))

base_theme<-theme_bw(base_size=AXIS_TEXT_SIZE)+theme(
 text=element_text(family="sans",face="bold",colour="black"),
 axis.title.y=element_text(size=AXIS_TITLE_SIZE,margin=margin(r=10)),
 axis.text.x=element_text(size=AXIS_TEXT_SIZE,colour="black",angle=42,hjust=1,vjust=1),
 axis.text.y=element_text(size=AXIS_TEXT_SIZE,colour="black"),
 axis.ticks=element_line(colour="black",linewidth=BORDER_LWD),
 panel.border=element_rect(colour="black",fill=NA,linewidth=BORDER_LWD),
 panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),aspect.ratio=1,
 plot.margin=margin(12,12,15,12))

panel_data<-list(
 A=list(data=dat[plot_role%in%c("mtOXPHOS","nuOXPHOS","assembly factor")],var="plot_role",
        levels=c("mtOXPHOS","nuOXPHOS","assembly factor"),labels=c(mtOXPHOS="mtOXPHOS",nuOXPHOS="nuOXPHOS",`assembly factor`="nu Assembly Factor"),pal=pal_overall,panel="A_OXPHOS_roles"),
 C=list(data=dat[plot_role%in%c("cyto-RP","Nmt-RP")],var="plot_role",levels=c("cyto-RP","Nmt-RP"),labels=c(`cyto-RP`="cyto-RP",`Nmt-RP`="Nmt-RP"),pal=pal_rp,panel="C_RP"),
 D=list(data=dat[plot_role%in%c("cyto-ARS","Nmt-ARS")],var="plot_role",levels=c("cyto-ARS","Nmt-ARS"),labels=c(`cyto-ARS`="cyto-ARS",`Nmt-ARS`="Nmt-ARS"),pal=pal_ars,panel="D_ARS"),
 E=list(data=dat[!is.na(core_plot_group)],var="core_plot_group",levels=c("mtOXPHOS","nu core","nu noncore"),labels=c(mtOXPHOS="mtOXPHOS",`nu core`="nu-core",`nu noncore`="nu-noncore"),pal=pal_core,panel="E_core_status"),
 F=list(data=dat[gene_source=="nuclear" & Kuster_class%in%c("direct_n-mt","indirect_n-mt","non-n-mt")],var="Kuster_class",levels=c("direct_n-mt","indirect_n-mt","non-n-mt"),labels=c(`direct_n-mt`="direct n-mt",`indirect_n-mt`="indirect n-mt",`non-n-mt`="non-n-mt"),pal=pal_kuster,panel="F_Kuster"))

make_simple<-function(spec,metric,point_cap=Inf){
 value<-metrics[[metric]]; d<-copy(spec$data)[is.finite(get(value))];d[,group:=factor(get(spec$var),levels=spec$levels)];d<-d[!is.na(group)]
 draw<-copy(d);if(is.finite(point_cap))draw<-draw[,if(.N>point_cap).SD[sample(.N,point_cap)] else .SD,by=group]
 metric_name<-metric
 st<-tests[panel==spec$panel & metric==metric_name & p_bh<.05]
 yr<-range(d[[value]],na.rm=TRUE);span<-diff(yr);if(!is.finite(span)||span==0)span<-1
 # Only the very large genome-wide panel is visually capped; statistics always use all values.
 cap<-if(spec$panel=="F_Kuster")as.numeric(quantile(d[[value]],.99,na.rm=TRUE)) else yr[2]
 cap<-max(cap,background[[metric]],na.rm=TRUE)
 visible_min<-min(0,min(d[[value]],na.rm=TRUE));visible_span<-cap-visible_min
 if(!is.finite(visible_span)||visible_span<=0)visible_span<-1
 base_y<-cap+.06*visible_span;span<-visible_span
 g<-ggplot(d,aes(group,.data[[value]],fill=group))+
  geom_hline(yintercept=background[[metric]],linetype="dotted",linewidth=1.05,colour="black")+
  geom_boxplot(width=.53,colour="black",linewidth=.85,outlier.shape=NA)+
  geom_point(data=draw,position=position_jitter(width=.11,seed=20260810),shape=21,size=POINT_SIZE,stroke=.3,colour="black",alpha=POINT_ALPHA)+
  scale_fill_manual(values=spec$pal,guide="none")+scale_x_discrete(labels=spec$labels,drop=FALSE)+
  labs(x=NULL,y=metric_lab[[metric]])+base_theme
 if(nrow(st)){
  yp<-base_y+seq_len(nrow(st))*.105*span
  g<-g+geom_signif(comparisons=Map(c,st$group1,st$group2),annotations=st$significance,y_position=yp,
      tip_length=.018,textsize=STAR_SIZE,size=.75,vjust=.15)
  top<-max(yp)+.1*span
 } else top<-base_y+.1*span
 g+coord_cartesian(ylim=c(min(0,yr[1]),top),clip="off")
}

make_complex<-function(metric){
 value<-metrics[[metric]];d<-copy(dat[plot_role%in%c("mtOXPHOS","nuOXPHOS","assembly factor") & !is.na(plot_complex) & is.finite(get(value))])
 d[,complex:=factor(as.character(plot_complex),levels=c("CI","CII","CIII","CIV","CV"))]
 d[,group:=factor(as.character(plot_role),levels=c("mtOXPHOS","nuOXPHOS","assembly factor"))]
 offsets<-c(mtOXPHOS=-.25,nuOXPHOS=0,`assembly factor`=.25); xbase<-setNames(seq_along(levels(d$complex)),levels(d$complex))
 ymax<-max(d[[value]],na.rm=TRUE);span<-diff(range(d[[value]],na.rm=TRUE));if(span==0)span<-1
 metric_name<-metric
 st<-tests[grepl("^B_",panel)&metric==metric_name&p_bh<.05]; seg<-list()
 if(nrow(st))for(i in seq_len(nrow(st))){cx<-sub("^B_","",st$panel[i]);lev<-sum(st$panel[1:i]==st$panel[i]);y<-ymax+(.13+.11*lev)*span
  seg[[i]]<-data.table(x1=xbase[cx]+offsets[st$group1[i]],x2=xbase[cx]+offsets[st$group2[i]],y=y,label=st$significance[i])}
 ann<-if(length(seg))rbindlist(seg) else data.table()
 g<-ggplot(d,aes(complex,.data[[value]],fill=group))+
  geom_hline(yintercept=background[[metric]],linetype="dotted",linewidth=1.05,colour="black")+
  geom_boxplot(position=position_dodge(.76),width=.62,colour="black",linewidth=.85,outlier.shape=NA)+
  geom_point(position=position_jitterdodge(jitter.width=.1,dodge.width=.76,seed=20260810),shape=21,size=1.7,stroke=.28,colour="black",alpha=POINT_ALPHA)+
  scale_fill_manual(values=pal_overall,labels=c(mtOXPHOS="mt",nuOXPHOS="nu",`assembly factor`="ass"))+
  labs(x=NULL,y=metric_lab[[metric]])+base_theme+theme(axis.text.x=element_text(size=AXIS_TEXT_SIZE,angle=0,hjust=.5),
   legend.position="top",legend.title=element_blank(),legend.text=element_text(size=20))
 if(nrow(ann))g<-g+geom_segment(data=ann,aes(x=x1,xend=x2,y=y,yend=y),inherit.aes=FALSE,linewidth=.75)+
  geom_segment(data=ann,aes(x=x1,xend=x1,y=y,yend=y-.025*span),inherit.aes=FALSE,linewidth=.75)+
  geom_segment(data=ann,aes(x=x2,xend=x2,y=y,yend=y-.025*span),inherit.aes=FALSE,linewidth=.75)+
  geom_text(data=ann,aes(x=(x1+x2)/2,y=y+.015*span,label=label),inherit.aes=FALSE,size=STAR_SIZE,fontface="bold")
 top<-if(nrow(ann))max(ann$y)+.1*span else ymax+.15*span;g+coord_cartesian(ylim=c(min(0,min(d[[value]],na.rm=TRUE)),top),clip="off")
}

make_row<-function(tag){spec<-panel_data[[tag]];cap<-if(tag=="F")1500 else Inf
 make_simple(spec,"piN",cap)|make_simple(spec,"piS",cap)|make_simple(spec,"piN_piS",cap)}
rows<-list(A=make_row("A"),B=make_complex("piN")|make_complex("piS")|make_complex("piN_piS"),
           C=make_row("C"),D=make_row("D"),E=make_row("E"),F=make_row("F"))
tag_row<-function(p,tag)p+plot_annotation(tag_levels=list(c(tag,"","")))&theme(plot.tag=element_text(size=ROW_TAG_SIZE,face="bold"))
rows<-Map(tag_row,rows,names(rows))
stubs<-c(A="Figure2A_pi_mt_nu_assembly",B="Figure2B_pi_complex",C="Figure2C_pi_RP",D="Figure2D_pi_ARS",E="Figure2E_pi_core_noncore",F="Figure2F_pi_Kuster")
for(tag in names(rows)){
 ggsave(file.path(outdir,paste0(stubs[tag],".png")),rows[[tag]],width=21,height=8.8,dpi=300,bg="white",limitsize=FALSE)
 ggsave(file.path(outdir,paste0(stubs[tag],".pdf")),rows[[tag]],width=21,height=8.8,bg="white",limitsize=FALSE)
}
combined<-wrap_plots(rows,ncol=1,heights=c(1,1.08,1,1,1,1))+plot_annotation(
 caption="Dotted lines indicate genome-wide nuclear medians. Stars show within-panel BH-adjusted Wilcoxon tests. For display only, panel F is capped at the 99th percentile; all observations were retained in statistical analyses.")
ggsave(file.path(outdir,"Figure2_pinpis_combined_A_to_F.png"),combined,width=21,height=53,dpi=300,bg="white",limitsize=FALSE)
ggsave(file.path(outdir,"Figure2_pinpis_combined_A_to_F.pdf"),combined,width=21,height=53,bg="white",limitsize=FALSE)
cat("[OK]",outdir,"\n");print(background)

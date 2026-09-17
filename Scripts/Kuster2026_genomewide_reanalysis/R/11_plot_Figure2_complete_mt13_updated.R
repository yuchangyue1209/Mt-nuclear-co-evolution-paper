#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggsignif)
  library(patchwork)
})

## ---------- paths ----------
root <- "/mnt/spareHD_2/genomewide_codeml_kuster"
nu_file <- file.path(root, "07_codeml_genomewide/codeml_master_analysis.tsv")
mt_file <- "/mnt/spareHD_2/genomewide_codeml_kuster/10_erc_mt13_codeml/mt13_M0_updated.tsv"
outdir <- file.path(root, "09_figures/Figure2_complete_updated_mt13_rerun")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

## ---------- visual constants ----------
BORDER_LWD <- 0.6
AXIS_TITLE_SIZE <- 30
AXIS_TEXT_SIZE <- 21
STAR_SIZE <- 8
ROW_TAG_SIZE <- 32
POINT_SIZE <- 1.8
POINT_ALPHA <- 0.60
SIGNIF_ONLY <- TRUE
ALPHA <- 0.05
set.seed(20260810)

pal_overall <- c(mt="#7FA4B2", nu="#BC8984", ass="#9EAD8D")
lab_overall <- c(mt="mtOXPHOS", nu="nuOXPHOS", ass="nu Assembly Factor")
pal_rp <- c(`cyto-ribo`="#7795A8", `Nmt-ribo`="#C28C73")
lab_rp <- c(`cyto-ribo`="cyto-RP", `Nmt-ribo`="Nmt-RP")
pal_ars <- c(`cyto-ARS`="#7589A9", `Nmt-ARS`="#B77D9A")
pal_core <- c(mt="#7FA4B2", nu_core="#BC8984", nu_noncore="#B6A174")
lab_core <- c(mt="mtOXPHOS", nu_core="nu-core", nu_noncore="nu-noncore")
pal_kuster <- c(`direct_n-mt`="#BC8984", `indirect_n-mt`="#948C9D", `non-n-mt`="#AAA8A1")
lab_kuster <- c(`direct_n-mt`="direct n-mt", `indirect_n-mt`="indirect n-mt", `non-n-mt`="non-n-mt")

## ---------- data ----------
nu <- read.delim(nu_file, check.names=FALSE)
nu$dN <- as.numeric(nu$tree_length_dN)
nu$dS <- as.numeric(nu$tree_length_dS)
nu$omega_plot <- as.numeric(nu$omega_ES1)
stopifnot(nrow(nu) == 17965)

mt0 <- read.delim(mt_file, check.names=FALSE)
mt0$role <- tolower(trimws(mt0$role))
mt <- mt0[mt0$model == "M0" & mt0$role == "mt", ]
mt$dN <- as.numeric(mt$dN); mt$dS <- as.numeric(mt$dS)
mt$omega_plot <- ifelse(as.numeric(mt$omega) < 999 & mt$dS >= 0.001, as.numeric(mt$omega), NA_real_)

normalize_complex <- function(z) {
  z <- toupper(gsub("[[:space:]]+", "", trimws(as.character(z))))
  out <- rep(NA_character_, length(z))
  out[z %in% c("CI","I")] <- "I"; out[z %in% c("CII","II")] <- "II"
  out[z %in% c("CIII","III")] <- "III"; out[z %in% c("CIV","IV")] <- "IV"
  out[z %in% c("CV","V")] <- "V"; out
}
nu$Complex <- normalize_complex(nu$own_complex)
mt$Complex <- normalize_complex(mt$complex)

background <- c(
  dN=median(nu$dN, na.rm=TRUE),
  dS=median(nu$dS, na.rm=TRUE),
  omega_plot=median(nu$omega_plot, na.rm=TRUE)
)

mk <- function(group, Complex=NA_character_, dN, dS, omega_plot) {
  data.frame(group, Complex, dN, dS, omega_plot, stringsAsFactors=FALSE)
}
mt_dat <- mk("mt", mt$Complex, mt$dN, mt$dS, mt$omega_plot)
nuclear_module <- ifelse(nu$own_role=="subunit", "nu", ifelse(nu$own_role=="assembly_factor", "ass", NA))
nu_module <- mk(nuclear_module, nu$Complex, nu$dN, nu$dS, nu$omega_plot)
A <- rbind(mt_dat, nu_module[!is.na(nu_module$group),]); B <- A[!is.na(A$Complex),]
B$Complex <- factor(B$Complex, levels=c("I","II","III","IV","V"))
C <- mk(nu$own_role, NA, nu$dN, nu$dS, nu$omega_plot); C <- C[C$group %in% c("cyto-ribo","Nmt-ribo"),]
D <- mk(nu$own_role, NA, nu$dN, nu$dS, nu$omega_plot); D <- D[D$group %in% c("cyto-ARS","Nmt-ARS"),]
core_group <- ifelse(nu$own_role=="subunit" & nu$core_status=="nu_core", "nu_core", ifelse(nu$own_role=="subunit" & nu$core_status=="nu_noncore", "nu_noncore", NA))
E <- rbind(mt_dat, mk(core_group, nu$Complex, nu$dN, nu$dS, nu$omega_plot)[!is.na(core_group),])
F <- mk(nu$Kuster_class, NA, nu$dN, nu$dS, nu$omega_plot)

## ---------- statistics and plotting ----------
stars <- function(p) ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*","n.s.")))
metric_label <- list(dN=expression(d[N]), dS=expression(d[S]), omega_plot=expression(d[N]/d[S]))

base_theme <- theme_bw(base_size=AXIS_TEXT_SIZE) + theme(
  text=element_text(family="sans", face="bold", colour="black"),
  axis.title.y=element_text(size=AXIS_TITLE_SIZE, margin=margin(r=10)),
  axis.text.x=element_text(size=AXIS_TEXT_SIZE, colour="black", angle=45, hjust=1),
  axis.text.y=element_text(size=AXIS_TEXT_SIZE, colour="black"),
  axis.ticks=element_line(colour="black", linewidth=BORDER_LWD),
  panel.border=element_rect(colour="black", fill=NA, linewidth=BORDER_LWD),
  panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(), aspect.ratio=1
)

simple_panel <- function(dat, metric, levels, palette, labels, row_name, cap=Inf) {
  d <- dat[is.finite(dat[[metric]]),]; d$group <- factor(d$group, levels=levels)
  pairs <- combn(levels, 2, simplify=FALSE)
  st <- do.call(rbind, lapply(pairs, function(pr) {
    a <- d[[metric]][d$group==pr[1]]; b <- d[[metric]][d$group==pr[2]]
    wt <- suppressWarnings(wilcox.test(a,b,exact=FALSE))
    data.frame(row=row_name, metric=metric, Complex=NA, group1=pr[1], group2=pr[2], n1=length(a), n2=length(b), median1=median(a), median2=median(b), mean1=mean(a), mean2=mean(b), p_raw=wt$p.value)
  }))
  st$p_adj <- p.adjust(st$p_raw,"BH"); st$stars <- stars(st$p_adj)
  draw <- if(SIGNIF_ONLY) st[st$p_adj<ALPHA,] else st
  points <- d
  if(is.finite(cap)) points <- do.call(rbind,lapply(split(d,d$group),function(z) if(nrow(z)>cap) z[sample(nrow(z),cap),] else z))
  yr <- range(d[[metric]],na.rm=TRUE); span <- diff(yr); if(span==0) span <- 1
  g <- ggplot(d,aes(group,.data[[metric]],fill=group)) +
    geom_hline(yintercept=background[[metric]],linetype="dotted",linewidth=1,colour="black") +
    geom_boxplot(width=.52,colour="black",linewidth=.8,outlier.shape=NA) +
    geom_point(data=points,position=position_jitter(width=.11,seed=20260810),shape=21,size=POINT_SIZE,stroke=.3,colour="black",alpha=POINT_ALPHA) +
    scale_fill_manual(values=palette,guide="none") + scale_x_discrete(labels=labels,drop=FALSE) +
    labs(x=NULL,y=metric_label[[metric]]) + base_theme
  if(nrow(draw)>0) {
    ypos <- yr[2]+seq_len(nrow(draw))*.13*span
    g <- g + geom_signif(comparisons=Map(c,draw$group1,draw$group2),annotations=draw$stars,y_position=ypos,tip_length=.02,textsize=STAR_SIZE,size=.7) +
      coord_cartesian(ylim=c(min(0,yr[1]),max(ypos)+.15*span),clip="off")
  }
  list(plot=g,stats=st)
}

complex_panel <- function(metric) {
  d <- B[is.finite(B[[metric]]),]; d$group <- factor(d$group,levels=c("mt","nu","ass"))
  st <- list(); k <- 1
  for(cx in levels(d$Complex)) {
    z <- d[d$Complex==cx,]; present <- c("mt","nu","ass")[c("mt","nu","ass") %in% unique(z$group)]
    if(length(present)<2) next
    tmp <- do.call(rbind,lapply(combn(present,2,simplify=FALSE),function(pr){
      a<-z[[metric]][z$group==pr[1]]; b<-z[[metric]][z$group==pr[2]]; wt<-suppressWarnings(wilcox.test(a,b,exact=FALSE))
      data.frame(row="B_complex",metric=metric,Complex=cx,group1=pr[1],group2=pr[2],n1=length(a),n2=length(b),median1=median(a),median2=median(b),mean1=mean(a),mean2=mean(b),p_raw=wt$p.value)
    })); tmp$p_adj<-p.adjust(tmp$p_raw,"BH"); tmp$stars<-stars(tmp$p_adj); st[[k]]<-tmp;k<-k+1
  }
  st <- do.call(rbind,st)
  g <- ggplot(d,aes(Complex,.data[[metric]],fill=group)) +
    geom_hline(yintercept=background[[metric]],linetype="dotted",linewidth=1,colour="black") +
    geom_boxplot(position=position_dodge(.76),width=.62,colour="black",linewidth=.8,outlier.shape=NA) +
    geom_point(position=position_jitterdodge(jitter.width=.1,dodge.width=.76,seed=20260810),shape=21,size=1.5,stroke=.25,colour="black",alpha=POINT_ALPHA) +
    scale_fill_manual(values=pal_overall,labels=c(mt="mt",nu="nu",ass="ass")) + labs(x=NULL,y=metric_label[[metric]]) + base_theme +
    theme(axis.text.x=element_text(size=AXIS_TEXT_SIZE,angle=0,hjust=.5),legend.position="top",legend.title=element_blank(),legend.text=element_text(size=20))
  list(plot=g,stats=st)
}

make_row <- function(dat,levels,pal,labs,name,cap=Inf) {
  r <- lapply(c("dN","dS","omega_plot"),function(m) simple_panel(dat,m,levels,pal,labs,name,cap)); names(r)<-c("dN","dS","omega")
  list(plot=r$dN$plot|r$dS$plot|r$omega$plot,stats=do.call(rbind,lapply(r,"[[","stats")))
}

rowA<-make_row(A,c("mt","nu","ass"),pal_overall,lab_overall,"A_overall")
rb<-lapply(c("dN","dS","omega_plot"),complex_panel); rowB<-list(plot=rb[[1]]$plot|rb[[2]]$plot|rb[[3]]$plot,stats=do.call(rbind,lapply(rb,"[[","stats")))
rowC<-make_row(C,c("cyto-ribo","Nmt-ribo"),pal_rp,lab_rp,"C_RP")
rowD<-make_row(D,c("cyto-ARS","Nmt-ARS"),pal_ars,c(`cyto-ARS`="cyto-ARS",`Nmt-ARS`="Nmt-ARS"),"D_ARS")
rowE<-make_row(E,c("mt","nu_core","nu_noncore"),pal_core,lab_core,"E_core")
rowF<-make_row(F,c("direct_n-mt","indirect_n-mt","non-n-mt"),pal_kuster,lab_kuster,"F_Kuster",1500)

tag_row <- function(p,tag) p + plot_annotation(tag_levels=list(c(tag,"",""))) & theme(plot.tag=element_text(size=ROW_TAG_SIZE,face="bold"))
plots <- list(A=tag_row(rowA$plot,"A"),B=tag_row(rowB$plot,"B"),C=tag_row(rowC$plot,"C"),D=tag_row(rowD$plot,"D"),E=tag_row(rowE$plot,"E"),F=tag_row(rowF$plot,"F"))

for(tag in names(plots)) {
  stub <- c(A="Figure2A_mt_nu_assembly",B="Figure2B_complex",C="Figure2C_RP",D="Figure2D_ARS",E="Figure2E_core_noncore",F="Figure2F_Kuster")[[tag]]
  ggsave(file.path(outdir,paste0(stub,".png")),plots[[tag]],width=21,height=8.5,dpi=300,bg="white")
  ggsave(file.path(outdir,paste0(stub,".pdf")),plots[[tag]],width=21,height=8.5,device=cairo_pdf,bg="white")
}

combined <- wrap_plots(plots,ncol=1,heights=c(1,1.08,1,1,1,1)) + plot_annotation(caption="Dotted lines show genome-wide medians. Invariant genes contribute dN=dS=0; omega includes identifiable ES>=1 estimates only.")
ggsave(file.path(outdir,"Figure2_combined_A_to_F.png"),combined,width=21,height=52,dpi=300,bg="white",limitsize=FALSE)
ggsave(file.path(outdir,"Figure2_combined_A_to_F.pdf"),combined,width=21,height=52,device=cairo_pdf,bg="white",limitsize=FALSE)

all_stats <- rbind(rowA$stats,rowB$stats,rowC$stats,rowD$stats,rowE$stats,rowF$stats)
write.table(all_stats,file.path(outdir,"Figure2_complete_statistics.tsv"),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
cat("[OK] Figure 2 ->",outdir,"\n"); print(background)


#!/usr/bin/env Rscript
# Rebuild manuscript panels from the archived statistical/plotting workflows.
# Arguments: nuclear codeml table, mitochondrial codeml table, pi workflow root,
# output directory. pi workflow root contains results/*combined.tsv and *tests.tsv.
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
a <- commandArgs(trailingOnly=TRUE)
if (length(a)!=4L) stop('Usage: Rscript 12_export_manuscript_Figure2_S2_S5.R NUCLEAR_CODEML_TSV MT_CODEML_TSV PINPIS_ROOT OUTPUT_DIR')
if (!all(file.exists(a[1:3]))) stop('Missing input file or pi workflow root')
out <- a[4]; dir.create(out,recursive=TRUE,showWarnings=FALSE)
script_arg <- grep('^--file=',commandArgs(),value=TRUE)
here <- dirname(normalizePath(sub('^--file=','',script_arg[1])))
dn_script <- file.path(here,'..','server_final','11_plot_figure2.R')
pi_script <- file.path(here,'pinpis','07_plot_pinpis_Figure2_large_fonts.R')
load_rows <- function(path,kind) {
 e <- new.env(parent=globalenv())
 for (expr in parse(path)) {
  assignment <- is.call(expr) && as.character(expr[[1]])[1] %in% c('<-','=')
  lhs <- if(assignment && is.symbol(expr[[2]])) as.character(expr[[2]]) else ''
  # Stop before legacy tags and exports; source calculations/row builders stay intact.
  if(kind=='dn' && lhs=='plot_A') break
  if(kind=='pi' && lhs=='rows' && is.call(expr[[3]]) && identical(expr[[3]][[1]],as.name('Map'))) break
  overrides <- if(kind=='dn') list(nuclear_file=a[1],old_mt_file=a[2],outdir=out,BORDER_LWD=.8) else list(root=a[3],outdir=out)
  if(lhs %in% names(overrides)) assign(lhs,overrides[[lhs]],envir=e) else eval(expr,envir=e)
 }
 if(kind=='dn') {
  setNames(lapply(LETTERS[1:6],function(tag) get(paste0('row_',tag),e)$plot),LETTERS[1:6])
 } else e$rows
}
dn <- load_rows(dn_script,'dn'); pi <- load_rows(pi_script,'pi')
frame <- theme(panel.border=element_rect(colour='black',fill=NA,linewidth=.8),axis.line=element_blank())
tagged <- function(p,tag) {
 p <- (p & frame)+plot_annotation(tag_levels=list(c(tag,'','')))
 p <- p & theme(plot.tag=element_text(size=32,face='bold',colour='black'))
 wrap_elements(full=p)
}
export <- function(p,name,height) {
 ggsave(file.path(out,paste0(name,'.pdf')),p,width=21,height=height,device=cairo_pdf,bg='white',limitsize=FALSE)
 ggsave(file.path(out,paste0(name,'.png')),p,width=21,height=height,dpi=300,bg='white',limitsize=FALSE)
}
export(wrap_plots(Map(tagged,dn[c('A','C','D')],c('A','B','C')),ncol=1),'Figure2_ACD_largefont',26.4)
export(wrap_plots(Map(tagged,pi[c('A','C','D')],c('A','B','C')),ncol=1),'FigureS2_ACD_largefont',26.4)
for (k in seq_along(c('B','E','F'))) {
 old <- c('B','E','F')[k]
 p <- wrap_plots(list(tagged(dn[[old]],'A'),tagged(pi[[old]],'B')),ncol=1)
 export(p,c('FigureS3_complex','FigureS4_core_noncore','FigureS5_direct_indirect_non')[k],17.6)
}
fwrite_map <- data.frame(figure=c('Figure2','FigureS2','FigureS3','FigureS4','FigureS5'),source_rows=c('dN A,C,D','pi A,C,D','dN B / pi B','dN E / pi E','dN F / pi F'))
write.table(fwrite_map,file.path(out,'manuscript_panel_map.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
cat('[OK] Manuscript figures saved to ',out,'\n',sep='')

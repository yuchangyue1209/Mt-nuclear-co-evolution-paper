平行
> library(ggplot2)

DELTA_FILE <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_SNPlevel_vs_mtPC_withAMO/deltaAF_long.withAMO.tsv.gz"
OUTDIR <- "/mnt/spareHD_2/nu_287/q2_parallelism/figs_deltaAF_pop_extreme"
dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

DEL <- read.table(gzfile(DELTA_FILE), header=TRUE, sep="\t")

POPSTAT <- aggregate(deltaAF ~ region + pop,
                     data=DEL,
                     FUN=function(x){
                       c(mean=mean(x, na.rm=TRUE),
                         mean_abs=mean(abs(x), na.rm=TRUE),
                         sd=sd(x, na.rm=TRUE),
                         n=length(x))
                     })

# 展开 list 列
POPSTAT <- do.call(data.frame, POPSTAT)
colnames(POPSTAT) <- c("region","pop",
                       "mean_deltaAF",
                       "mean_abs_deltaAF",
                       "sd_deltaAF",
                       "n_snps")

POPSTAT <- POPSTAT[order(-POPSTAT$mean_abs_deltaAF), ]

write.table(POPSTAT,
            file=file.path(OUTDIR,"POP_deltaAF_summary.tsv"),
            sep="\t",
            row.names=FALSE,
            quote=FALSE)

head(POPSTAT)
   region  pop  mean_deltaAF mean_abs_deltaAF sd_deltaAF n_snps
7      BC  JOE -0.0002997228      0.002756867 0.03623437  42634
4      BC ECHO -0.0002287575      0.002599823 0.03662072  42592
18     AK   TL  0.0004880247      0.002376227 0.02936164  42556
21     AK   WT  0.0004924862      0.002365655 0.02939016  42649
16     BC  SWA -0.0001374537      0.002360428 0.03267996  42549
14     AK   SL  0.0006115779      0.002333856 0.02884917  42668
> library(ggplot2)

DELTA_FILE <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_SNPlevel_vs_mtPC_withAMO/deltaAF_long.withAMO.tsv.gz"
OUTDIR <- "/mnt/spareHD_2/nu_287/q2_parallelism/figs_deltaAF_pop_extreme"
dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

DEL <- read.table(gzfile(DELTA_FILE), header=TRUE, sep="\t")

POPSTAT <- aggregate(deltaAF ~ region + pop,
                     data=DEL,
                     FUN=function(x){
                       c(mean=mean(x, na.rm=TRUE),
                         mean_abs=mean(abs(x), na.rm=TRUE),
                         sd=sd(x, na.rm=TRUE),
                         n=length(x))
                     })

# 展开 list 列
POPSTAT <- do.call(data.frame, POPSTAT)
colnames(POPSTAT) <- c("region","pop",
                       "mean_deltaAF",
                       "mean_abs_deltaAF",
                       "sd_deltaAF",
                       "n_snps")

POPSTAT <- POPSTAT[order(-POPSTAT$mean_abs_deltaAF), ]

write.table(POPSTAT,
            file=file.path(OUTDIR,"POP_deltaAF_summary.tsv"),
            sep="\t",
            row.names=FALSE,
            quote=FALSE)

head(POPSTAT)
   region  pop  mean_deltaAF mean_abs_deltaAF sd_deltaAF n_snps
7      BC  JOE -0.0002997228      0.002756867 0.03623437  42634
4      BC ECHO -0.0002287575      0.002599823 0.03662072  42592
18     AK   TL  0.0004880247      0.002376227 0.02936164  42556
21     AK   WT  0.0004924862      0.002365655 0.02939016  42649
16     BC  SWA -0.0001374537      0.002360428 0.03267996  42549
14     AK   SL  0.0006115779      0.002333856 0.02884917  42668
> aggregate(mean_abs_deltaAF ~ region, data=POPSTAT, mean)
  region mean_abs_deltaAF
1     AK      0.002213191
2     BC      0.002168591
> POPSTAT$z_strength <- scale(POPSTAT$mean_abs_deltaAF)
POPSTAT$z_sd <- scale(POPSTAT$sd_deltaAF)
> 

Save workspace image? [y/n/c]: n

(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism$ zcat /mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_SNPlevel_vs_mtPC_noAMO/deltaAF_long.noAMO.tsv.gz \
| head -n 1 \
| tr '\t' '\n' \
| nl -ba
     1	region
     2	snp
     3	chr
     4	pos
     5	gene
     6	pop
     7	deltaAF
     8	af
     9	marine_af
    10	treePC1
    11	treePC2
    12	mitoPC1
    13	mitoPC2
    14	mitoPC3
    15	mitoPC4
    16	mitoPC5
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism$ cat > /mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_from_af72.R <<'RS'
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

# ========= Inputs =========
AF_FILE <- Sys.getenv("AF_FILE", unset="/mnt/spareHD_2/nu_287/q2_parallelism/af_long_final_72genes_subunit_with_si.tsv.gz")
OUTDIR  <- Sys.getenv("OUTDIR",  unset="/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes")
dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

# ========= Params (per your Methods) =========
MIN_DEPTH <- as.integer(Sys.getenv("MIN_DEPTH", unset="10"))
EPS_ZERO  <- as.numeric(Sys.getenv("EPS_ZERO",  unset="0.01"))

MIN_AVAIL_AK <- as.integer(Sys.getenv("MIN_AVAIL_AK", unset="5"))  # 5/9
MIN_AVAIL_BC <- as.integer(Sys.getenv("MIN_AVAIL_BC", unset="8"))  # 8/13

MIN_NONZERO_AK <- as.integer(Sys.getenv("MIN_NONZERO_AK", unset="5"))
MIN_NONZERO_BC <- as.integer(Sys.getenv("MIN_NONZERO_BC", unset="8"))

GENE_MIN_SNPS <- as.integer(Sys.getenv("GENE_MIN_SNPS", unset="10"))
SNP_CONC_THR  <- as.numeric(Sys.getenv("SNP_CONC_THR",  unset="0.70"))
GENE_PROP_THR <- as.numeric(Sys.getenv("GENE_PROP_THR", unset="0.30"))

# ========= Pop groups =========
AK_fresh  <- c("FG","LG","SR","SL","TL","WB","WT","WK","LB")
BC_fresh  <- c("SWA","THE","JOE","BEA","MUC","PYE","ROS","AMO","BOOT","ECHO","LAW","GOS","ROB")
AK_marine <- "RS"
BC_marine <- "SAY"
keep_pops <- unique(c(AK_fresh, AK_marine, BC_fresh, BC_marine))

normalize_pop <- function(x){
  x <- toupper(x)
  gsub("^(\\d+_)?([A-Z]+)(?:_S\\d+)?$", "\\2", x, perl=TRUE)
}

RSt("[OK] wrote outputs to: ", OUTDIR, "\n", sep="")ional_parallel_genes_AK_intersect_BC.tsv"))h concordance>=", SNP_CONC_THR),]
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism$ AF_FILE="/mnt/spareHD_2/nu_287/q2_parallelism/af_long_final_72genes_subunit_with_si.tsv.gz" \
OUTDIR="/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes" \
Rscript /mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_from_af72.R
[OK] wrote outputs to: /mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism$ cut -f1-6 /mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes/perGene_withinRegion_parallel.tsv | head
wc -l /mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes/crossRegional_parallel_genes_AK_intersect_BC.tsv
region	gene	n_snps	n_snps_conc70	prop_conc70	gene_maj_sign
AK	ATP5F1B	26	15	0.576923076923077	+
BC	ATP5F1B	38	26	0.68421052631579	+
BC	COX4I1	15	14	0.933333333333333	+
AK	NDUFS3	20	11	0.55	+
BC	NDUFS3	19	13	0.68421052631579	+
BC	UQCRFS1	16	16	1	+
AK	SDHC	10	8	0.8	+
BC	SDHC	11	10	0.909090909090909	+
BC	ATP5F1C	15	14	0.933333333333333	+
14 /mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes/crossRegional_parallel_genes_AK_intersect_BC.tsv
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism$ zcat perSNP_withinRegion_concordance.tsv.gz \
| awk -F'\t' 'NR>1{print $10}' \
| tr ',' '\n' \
| sort | uniq -c | sort -nr
gzip: perSNP_withinRegion_concordance.tsv.gz: No such file or directory
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism$ zcat /mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes/perSNP_withinRegion_concordance.tsv.gz \
| awk -F'\t' 'NR>1{print $10}' \
| tr ',' '\n' \
| sort | uniq -c | sort -nr
    806 BOOT
    803 ECHO
    802 MUC
    800 GOS
    796 SWA
    794 PYE
    792 LAW
    792 BEA
    792 AMO
    791 THE
    788 JOE
    781 ROB
    394 WK
    392 TL
    391 SL
    390 FG
    387 SR
    382 WT
    382 LB
    381 WB
    381 LG
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism$  zcat /mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes/perSNP_withinRegion_concordance.tsv.gz \
| awk -F'\t' 'NR>1{print $10}' \
| tr ',' '\n' \
| sort | uniq -c | sort -nr
    806 BOOT
    803 ECHO
    802 MUC
    800 GOS
    796 SWA
    794 PYE
    792 LAW
    792 BEA
    792 AMO
    791 THE
    788 JOE
    781 ROB
    394 WK
    392 TL
    391 SL
    390 FG
    387 SR
    382 WT
    382 LB
    381 WB
    381 LG
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism$ 
    806 BOOT
    803 ECHO
    802 MUC
    800 GOS
    796 SWA
    794 PYE
    792 LAW
    792 BEA
    792 AMO
    791 THE
    788 JOE
    781 ROB
    394 WK
    392 TL
    391 SL
    390 FG
    387 SR
    382 WT
    382 LB
    381 WB
    381 LG
806: command not found
803: command not found
802: command not found
800: command not found
796: command not found
794: command not found
792: command not found
792: command not found
792: command not found
791: command not found
Command '788' not found, did you mean:
  command 'z88' from deb z88 (13.0.0+dfsg2-6)
Try: apt install <deb name>
781: command not found
394: command not found
392: command not found
391: command not found
390: command not found
387: command not found
382: command not found
382: command not found
381: command not found
381: command not found
-bash: syntax error near unexpected token `cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism$'
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism$ cd SNPF="/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes/perSNP_withinRegion_concordance.tsv.gz"
GENEF="/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes/perGene_withinRegion_parallel.tsv"
CROSS="/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes/crossRegional_parallel_genes_AK_intersect_BC.tsv"
-bash: cd: SNPF=/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes/perSNP_withinRegion_concordance.tsv.gz: No such file or directory
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism$ cd /mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes/
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ SNPF="/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes/perSNP_withinRegion_concordance.tsv.gz"
GENEF="/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes/perGene_withinRegion_parallel.tsv"
CROSS="/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes/crossRegional_parallel_genes_AK_intersect_BC.tsv"
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ zcat "$SNPF" \
| awk -F'\t' 'NR>1 && $1=="AK"{print $10}' \
| tr ',' '\n' \
| sort | uniq -c | sort -nr
    394 WK
    392 TL
    391 SL
    390 FG
    387 SR
    382 WT
    382 LB
    381 WB
    381 LG
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ zcat "$SNPF" \
| awk -F'\t' 'NR>1 && $1=="BC"{print $10}' \
| tr ',' '\n' \
| sort | uniq -c | sort -nr
    806 BOOT
    803 ECHO
    802 MUC
    800 GOS
    796 SWA
    794 PYE
    792 LAW
    792 BEA
    792 AMO
    791 THE
    788 JOE
    781 ROB
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ awk -F'\t' 'NR>1 && $1=="AK" && $7=="TRUE"{print $2}' "$GENEF" | sort -u > /mnt/spareHD_2/nu_287/q2_parallelism/AK_parallel_genes.list
wc -l /mnt/spareHD_2/nu_287/q2_parallelism/AK_parallel_genes.list
0 /mnt/spareHD_2/nu_287/q2_parallelism/AK_parallel_genes.list
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ awk -F'\t' 'NR>1 && $1=="BC" && $7=="TRUE"{print $2}' "$GENEF" | sort -u > /mnt/spareHD_2/nu_287/q2_parallelism/BC_parallel_genes.list
wc -l /mnt/spareHD_2/nu_287/q2_parallelism/BC_parallel_genes.list
0 /mnt/spareHD_2/nu_287/q2_parallelism/BC_parallel_genes.list
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$  awk -F'\t' 'NR>1 && $1=="AK" && $7=="TRUE"{print $2}' "$GENEF" | sort -u > /mnt/spareHD_2/nu_287/q2_parallelism/AK_parallel_genes.list
wc -l /mnt/spareHD_2/nu_287/q2_parallelism/AK_parallel_genes.list
0 /mnt/spareHD_2/nu_287/q2_parallelism/AK_parallel_genes.list
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ awk -F'\t' 'NR>1 && $1=="BC" && $7=="TRUE"{print $2}' "$GENEF" | sort -u > /mnt/spareHD_2/nu_287/q2_parallelism/BC_parallel_genes.list
wc -l /mnt/spareHD_2/nu_287/q2_parallelism/BC_parallel_genes.list
0 /mnt/spareHD_2/nu_287/q2_parallelism/BC_parallel_genes.list
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ 
0 /mnt/spareHD_2/nu_287/q2_parallelism/AK_parallel_genes.list
0: command not found
-bash: syntax error near unexpected token `cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$'
0 /mnt/spareHD_2/nu_287/q2_parallelism/BC_parallel_genes.list
0: command not found
-bash: syntax error near unexpected token `cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$'
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ GENEF="/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes/perGene_withinRegion_parallel.tsv"
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ awk -F'\t' 'NR>1 && $1=="AK" && $3>=10 && $5>=0.30 {print $2}' "$GENEF" \
| sort -u > /mnt/spareHD_2/nu_287/q2_parallelism/AK_parallel_genes.list

wc -l /mnt/spareHD_2/nu_287/q2_parallelism/AK_parallel_genes.list
head /mnt/spareHD_2/nu_287/q2_parallelism/AK_parallel_genes.list
13 /mnt/spareHD_2/nu_287/q2_parallelism/AK_parallel_genes.list
ATP5F1A
ATP5F1B
CYC1
HCCSB
NDUFA10
NDUFS1
NDUFS2
NDUFS3
NDUFS8
SDHA
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ awk -F'\t' 'NR>1 && $1=="BC" && $3>=10 && $5>=0.30 {print $2}' "$GENEF" \
| sort -u > /mnt/spareHD_2/nu_287/q2_parallelism/BC_parallel_genes.list

wc -l /mnt/spareHD_2/nu_287/q2_parallelism/BC_parallel_genes.list
head /mnt/spareHD_2/nu_287/q2_parallelism/BC_parallel_genes.list
26 /mnt/spareHD_2/nu_287/q2_parallelism/BC_parallel_genes.list
ATP5F1A
ATP5F1B
ATP5F1C
COX4I1
CYC1
DMAC2L
HCCSB
NDUFA10
NDUFA4
NDUFAB1B
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ comm -12 /mnt/spareHD_2/nu_287/q2_parallelism/AK_parallel_genes.list \
         /mnt/spareHD_2/nu_287/q2_parallelism/BC_parallel_genes.list \
> /mnt/spareHD_2/nu_287/q2_parallelism/AK_BC_intersect_genes.list

wc -l /mnt/spareHD_2/nu_287/q2_parallelism/AK_BC_intersect_genes.list
head /mnt/spareHD_2/nu_287/q2_parallelism/AK_BC_intersect_genes.list
13 /mnt/spareHD_2/nu_287/q2_parallelism/AK_BC_intersect_genes.list
ATP5F1A
ATP5F1B
CYC1
HCCSB
NDUFA10
NDUFS1
NDUFS2
NDUFS3
NDUFS8
SDHA
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ SNPF="/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes/perSNP_withinRegion_concordance.tsv.gz"
LIST="/mnt/spareHD_2/nu_287/q2_parallelism/AK_BC_intersect14_genes.list"
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ zcat "$SNPF" \
| awk -F'\t' 'NR==1{next} $1=="AK"{print $4"\t"$10}' \
| awk -F'\t' 'BEGIN{while((getline<L)>0) g[$1]=1} (g[$1]==1){print $2}' L="$LIST" \
| tr ',' '\n' | sort | uniq -c | sort -nr
awk: cmd. line:1: fatal: expression for `<' redirection has null string value
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ zcat "$SNPF" \
| awk -F'\t' 'NR==1{next} $1=="AK"{print $4"\t"$10}' \
| awk -F'\t' 'BEGIN{while((getline<L)>0) g[$1]=1} (g[$1]==1){print $2}' L="$LIST" \
| tr ',' '\n' | sort | uniq -c | sort -nr
awk: cmd. line:1: fatal: expression for `<' redirection has null string value
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ zcat /mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes/perSNP_withinRegion_concordance.tsv.gz \
| awk -F'\t' 'NR==1{next} $1=="AK"{print $4"\t"$10}' \
| awk -F'\t' 'BEGIN{
    while((getline < "/mnt/spareHD_2/nu_287/q2_parallelism/AK_BC_intersect14_genes.list")>0)
        g[$1]=1
} (g[$1]==1){print $2}' \
| tr ',' '\n' | sort | uniq -c | sort -nr
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ zcat /mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes/perSNP_withinRegion_concordance.tsv.gz \
| awk -F'\t' 'NR==1{next} $1=="BC"{print $4"\t"$10}' \
| awk -F'\t' 'BEGIN{
    while((getline < "/mnt/spareHD_2/nu_287/q2_parallelism/AK_BC_intersect14_genes.list")>0)
        g[$1]=1
} (g[$1]==1){print $2}' \
| tr ',' '\n' | sort | uniq -c | sort -nr
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ zcat /mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes/perSNP_withinRegion_concordance.tsv.gz \
| head -n 1 \
| tr '\t' '\n' \
| nl -ba
     1	region
     2	chr
     3	pos
     4	gene
     5	snp
     6	n_nonzero
     7	maj_sign
     8	concordance
     9	mean_abs_delta
    10	pops_nonzero
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ tr -d '\r' < /mnt/spareHD_2/nu_287/q2_parallelism/AK_BC_intersect14_genes.list \
| awk '{gsub(/^[ \t]+|[ \t]+$/,""); print toupper($0)}' \
| sort -u \
> /mnt/spareHD_2/nu_287/q2_parallelism/AK_BC_intersect14_genes.clean.list

wc -l /mnt/spareHD_2/nu_287/q2_parallelism/AK_BC_intersect14_genes.clean.list
head /mnt/spareHD_2/nu_287/q2_parallelism/AK_BC_intersect14_genes.clean.list
-bash: /mnt/spareHD_2/nu_287/q2_parallelism/AK_BC_intersect14_genes.list: No such file or directory
0 /mnt/spareHD_2/nu_287/q2_parallelism/AK_BC_intersect14_genes.clean.list
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ tr -d '\r' < /mnt/spareHD_2/nu_287/q2_parallelism/AK_BC_intersect14_genes.list \
| awk '{gsub(/^[ \t]+|[ \t]+$/,""); print toupper($0)}' \
| sort -u \
> /mnt/spareHD_2/nu_287/q2_parallelism/AK_BC_intersect14_genes.clean.list

wc -l /mnt/spareHD_2/nu_287/q2_parallelism/AK_BC_intersect14_genes.clean.list
head /mnt/spareHD_2/nu_287/q2_parallelism/AK_BC_intersect14_genes.clean.list
-bash: /mnt/spareHD_2/nu_287/q2_parallelism/AK_BC_intersect14_genes.list: No such file or directory
0 /mnt/spareHD_2/nu_287/q2_parallelism/AK_BC_intersect14_genes.clean.list
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ CROSS="/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes/crossRegional_parallel_genes_AK_intersect_BC.tsv"

awk -F'\t' 'NR>1{print toupper($1)}' "$CROSS" \
| tr -d '\r' \
| awk '{gsub(/^[ \t]+|[ \t]+$/,""); print}' \
| sort -u \
> /mnt/spareHD_2/nu_287/q2_parallelism/AK_BC_intersect14_genes.clean.list

wc -l /mnt/spareHD_2/nu_287/q2_parallelism/AK_BC_intersect14_genes.clean.list
head /mnt/spareHD_2/nu_287/q2_parallelism/AK_BC_intersect14_genes.clean.list
13 /mnt/spareHD_2/nu_287/q2_parallelism/AK_BC_intersect14_genes.clean.list
ATP5F1A
ATP5F1B
CYC1
HCCSB
NDUFA10
NDUFS1
NDUFS2
NDUFS3
NDUFS8
SDHA
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ SNPF="/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes/perSNP_withinRegion_concordance.tsv.gz"
LIST="/mnt/spareHD_2/nu_287/q2_parallelism/AK_BC_intersect14_genes.clean.list"

zcat "$SNPF" \
| awk -F'\t' 'NR>1 && $1=="AK"{g=toupper($4); sub(/\r$/,"",g); print g"\t"$10}' \
| awk -F'\t' 'BEGIN{while((getline < L)>0) keep[$1]=1} (keep[$1]==1){print $2}' L="$LIST" \
| tr ',' '\n' | sort | uniq -c | sort -nr
awk: cmd. line:1: fatal: expression for `<' redirection has null string value
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ zcat "$SNPF" \
| awk -F'\t' 'NR>1 && $1=="BC"{g=toupper($4); sub(/\r$/,"",g); print g"\t"$10}' \
| awk -F'\t' 'BEGIN{while((getline < L)>0) keep[$1]=1} (keep[$1]==1){print $2}' L="$LIST" \
| tr ',' '\n' | sort | uniq -c | sort -nr
awk: cmd. line:1: fatal: expression for `<' redirection has null string value
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ ls -l /mnt/spareHD_2/nu_287/q2_parallelism/AK_BC_intersect14_genes.clean.list
-rw-rw-r-- 1 cyu cyu 88 Feb 24 21:30 /mnt/spareHD_2/nu_287/q2_parallelism/AK_BC_intersect14_genes.clean.list
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ zcat /mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes/perSNP_withinRegion_concordance.tsv.gz \
| awk -F'\t' 'NR>1 && $1=="AK"{g=toupper($4); sub(/\r$/,"",g); print g"\t"$10}' \
| awk -F'\t' 'BEGIN{
    while((getline < "/mnt/spareHD_2/nu_287/q2_parallelism/AK_BC_intersect14_genes.clean.list")>0)
        keep[$1]=1
} (keep[$1]==1){print $2}' \
| tr ',' '\n' | sort | uniq -c | sort -nr
    208 WK
    208 TL
    208 SL
    205 SR
    204 FG
    201 WT
    201 WB
    201 LG
    199 LB
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ zcat /mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes/perSNP_withinRegion_concordance.tsv.gz \
| awk -F'\t' 'NR>1 && $1=="AK"{g=toupper($4); sub(/\r$/,"",g); print g"\t"$10}' \
| awk -F'\t' 'BEGIN{
    while((getline < "/mnt/spareHD_2/nu_287/q2_parallelism/AK_BC_intersect14_genes.clean.list")>0)
        keep[$1]=1
} (keep[$1]==1){print $2}' \
| tr ',' '\n' | sort | uniq -c | sort -nr
    208 WK
    208 TL
    208 SL
    205 SR
    204 FG
    201 WT
    201 WB
    201 LG
    199 LB
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ zcat /mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes/perSNP_withinRegion_concordance.tsv.gz \
| awk -F'\t' 'NR>1 && $1=="BC"{g=toupper($4); sub(/\r$/,"",g); print g"\t"$10}' \
| awk -F'\t' 'BEGIN{
    while((getline < "/mnt/spareHD_2/nu_287/q2_parallelism/AK_BC_intersect14_genes.clean.list")>0)
        keep[$1]=1
} (keep[$1]==1){print $2}' \
| tr ',' '\n' | sort | uniq -c | sort -nr
    380 MUC
    380 ECHO
    380 BOOT
    379 GOS
    376 SWA
    376 LAW
    375 PYE
    373 BEA
    373 AMO
    371 THE
    370 ROB
    370 JOE
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ SNPF="/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes/perSNP_withinRegion_concordance.tsv.gz"
GENEF="/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes/perGene_withinRegion_parallel.tsv"
CROSS="/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes/crossRegional_parallel_genes_AK_intersect_BC.tsv"
OUTD="/mnt/spareHD_2/nu_287/q2_parallelism/q2_tables_gene_snp_pop"
mkdir -p "$OUTD"
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ awk -F'\t' 'NR>1 && $1=="AK" && $3>=10 && $5>=0.30 {print toupper($2)}' "$GENEF" \
| sort -u > "$OUTD/AK_parallel_genes.list"

wc -l "$OUTD/AK_parallel_genes.list"
13 /mnt/spareHD_2/nu_287/q2_parallelism/q2_tables_gene_snp_pop/AK_parallel_genes.list
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ awk -F'\t' 'NR>1 && $1=="BC" && $3>=10 && $5>=0.30 {print toupper($2)}' "$GENEF" \
| sort -u > "$OUTD/BC_parallel_genes.list"

wc -l "$OUTD/BC_parallel_genes.list"
26 /mnt/spareHD_2/nu_287/q2_parallelism/q2_tables_gene_snp_pop/BC_parallel_genes.list
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ awk -F'\t' 'NR>1{print toupper($1)}' "$CROSS" | tr -d '\r' | sort -u \
> "$OUTD/INTERSECT_genes.list"

wc -l "$OUTD/INTERSECT_genes.list"
head "$OUTD/INTERSECT_genes.list"
13 /mnt/spareHD_2/nu_287/q2_parallelism/q2_tables_gene_snp_pop/INTERSECT_genes.list
ATP5F1A
ATP5F1B
CYC1
HCCSB
NDUFA10
NDUFS1
NDUFS2
NDUFS3
NDUFS8
SDHA
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ zcat "$SNPF" \
| awk -F'\t' 'NR==1{print "region\tgene\tsnp\tchr\tpos\tmaj_sign\tconcordance\tmean_abs_delta\tpops_nonzero"; next}
             $1=="AK"{g=toupper($4); sub(/\r$/,"",g); print $1"\t"g"\t"$5"\t"$2"\t"$3"\t"$7"\t"$8"\t"$9"\t"$10}' \
| awk -F'\t' 'BEGIN{while((getline<LIST)>0) keep[$1]=1}
             NR==1 || keep[$2]{print}' LIST="$OUTD/AK_parallel_genes.list" \
> "$OUTD/AK_gene_snp_pops.tsv"

wc -l "$OUTD/AK_gene_snp_pops.tsv"
head "$OUTD/AK_gene_snp_pops.tsv"
awk: cmd. line:1: fatal: expression for `<' redirection has null string value
0 /mnt/spareHD_2/nu_287/q2_parallelism/q2_tables_gene_snp_pop/AK_gene_snp_pops.tsv
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ zcat /mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes/perSNP_withinRegion_concordance.tsv.gz \
| awk -F'\t' 'NR==1{print "region\tgene\tsnp\tchr\tpos\tmaj_sign\tconcordance\tmean_abs_delta\tpops_nonzero"; next}
             $1=="AK"{g=toupper($4); sub(/\r$/,"",g); print $1"\t"g"\t"$5"\t"$2"\t"$3"\t"$7"\t"$8"\t"$9"\t"$10}' \
| awk -F'\t' 'BEGIN{
    while((getline < "/mnt/spareHD_2/nu_287/q2_parallelism/q2_tables_gene_snp_pop/AK_parallel_genes.list")>0)
      keep[$1]=1
} NR==1 || keep[$2]{print}' \
> /mnt/spareHD_2/nu_287/q2_parallelism/q2_tables_gene_snp_pop/AK_gene_snp_pops.tsv

wc -l /mnt/spareHD_2/nu_287/q2_parallelism/q2_tables_gene_snp_pop/AK_gene_snp_pops.tsv
head /mnt/spareHD_2/nu_287/q2_parallelism/q2_tables_gene_snp_pop/AK_gene_snp_pops.tsv
219 /mnt/spareHD_2/nu_287/q2_parallelism/q2_tables_gene_snp_pop/AK_gene_snp_pops.tsv
region	gene	snp	chr	pos	maj_sign	concordance	mean_abs_delta	pops_nonzero
AK	ATP5F1B	chrI:27378199:ATP5F1B	chrI	27378199	+	0.875	0.058100375	FG,LB,LG,SR,TL,WB,WT,WK
AK	ATP5F1B	chrI:27378252:ATP5F1B	chrI	27378252	+	1	0.110687555555556	FG,LB,LG,SR,SL,TL,WB,WT,WK
AK	ATP5F1B	chrI:27378275:ATP5F1B	chrI	27378275	+	1	0.0231479999999999	FG,LB,LG,SR,SL,TL,WB,WT,WK
AK	ATP5F1B	chrI:27379018:ATP5F1B	chrI	27379018	+	0.555555555555556	0.144910333333333	FG,LB,LG,SR,SL,TL,WB,WT,WK
AK	ATP5F1B	chrI:27379326:ATP5F1B	chrI	27379326	+	0.875	0.270958875	LB,LG,SR,SL,TL,WB,WT,WK
AK	ATP5F1B	chrI:27379389:ATP5F1B	chrI	27379389	-	0.666666666666667	0.143726777777778	FG,LB,LG,SR,SL,TL,WB,WT,WK
AK	ATP5F1B	chrI:27379410:ATP5F1B	chrI	27379410	-	0.625	0.2723375	FG,LB,LG,SR,TL,WB,WT,WK
AK	ATP5F1B	chrI:27379476:ATP5F1B	chrI	27379476	+	1	0.151396444444444	FG,LB,LG,SR,SL,TL,WB,WT,WK
AK	ATP5F1B	chrI:27379690:ATP5F1B	chrI	27379690	+	1	0.012346	FG,LB,SR,SL,TL,WT,WK
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ zcat /mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes/perSNP_withinRegion_concordance.tsv.gz \
| awk -F'\t' 'NR==1{print "region\tgene\tsnp\tchr\tpos\tmaj_sign\tconcordance\tmean_abs_delta\tpops_nonzero"; next}
             $1=="BC"{g=toupper($4); sub(/\r$/,"",g); print $1"\t"g"\t"$5"\t"$2"\t"$3"\t"$7"\t"$8"\t"$9"\t"$10}' \
| awk -F'\t' 'BEGIN{
    while((getline < "/mnt/spareHD_2/nu_287/q2_parallelism/q2_tables_gene_snp_pop/BC_parallel_genes.list")>0)
      keep[$1]=1
} NR==1 || keep[$2]{print}' \
> /mnt/spareHD_2/nu_287/q2_parallelism/q2_tables_gene_snp_pop/BC_gene_snp_pops.tsv

wc -l /mnt/spareHD_2/nu_287/q2_parallelism/q2_tables_gene_snp_pop/BC_gene_snp_pops.tsv
head /mnt/spareHD_2/nu_287/q2_parallelism/q2_tables_gene_snp_pop/BC_gene_snp_pops.tsv
576 /mnt/spareHD_2/nu_287/q2_parallelism/q2_tables_gene_snp_pop/BC_gene_snp_pops.tsv
region	gene	snp	chr	pos	maj_sign	concordance	mean_abs_delta	pops_nonzero
BC	ATP5F1B	chrI:27378199:ATP5F1B	chrI	27378199	+	1	0.0553936666666667	THE,JOE,BEA,MUC,PYE,AMO,GOS,ROB,BOOT,ECHO,LAW,SWA
BC	ATP5F1B	chrI:27378204:ATP5F1B	chrI	27378204	+	1	0.013889	THE,JOE,BEA,MUC,PYE,AMO,GOS,ROB,BOOT,ECHO,LAW,SWA
BC	ATP5F1B	chrI:27378216:ATP5F1B	chrI	27378216	+	1	0.0438988333333334	THE,JOE,BEA,MUC,PYE,AMO,GOS,ROB,BOOT,ECHO,LAW,SWA
BC	ATP5F1B	chrI:27378252:ATP5F1B	chrI	27378252	+	0.666666666666667	0.0374903333333333	THE,JOE,BEA,MUC,PYE,AMO,GOS,ROB,BOOT,ECHO,LAW,SWA
BC	ATP5F1B	chrI:27379018:ATP5F1B	chrI	27379018	+	0.666666666666667	0.216926	THE,JOE,BEA,MUC,PYE,AMO,GOS,ROB,BOOT,ECHO,LAW,SWA
BC	ATP5F1B	chrI:27379326:ATP5F1B	chrI	27379326	-	0.583333333333333	0.249932	THE,JOE,BEA,MUC,PYE,AMO,GOS,ROB,BOOT,ECHO,LAW,SWA
BC	ATP5F1B	chrI:27379389:ATP5F1B	chrI	27379389	+	0.727272727272727	0.309179545454545	THE,JOE,BEA,MUC,PYE,AMO,GOS,BOOT,ECHO,LAW,SWA
BC	ATP5F1B	chrI:27379410:ATP5F1B	chrI	27379410	+	0.583333333333333	0.242088416666667	THE,JOE,BEA,MUC,PYE,AMO,GOS,ROB,BOOT,ECHO,LAW,SWA
BC	ATP5F1B	chrI:27379437:ATP5F1B	chrI	27379437	+	1	0.027397	THE,JOE,BEA,MUC,PYE,AMO,GOS,ROB,BOOT,ECHO,LAW,SWA
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ zcat /mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes/perSNP_withinRegion_concordance.tsv.gz \
| awk -F'\t' 'NR==1{print "region\tgene\tsnp\tchr\tpos\tmaj_sign\tconcordance\tmean_abs_delta\tpops_nonzero"; next}
             ($1=="AK" || $1=="BC"){g=toupper($4); sub(/\r$/,"",g); print $1"\t"g"\t"$5"\t"$2"\t"$3"\t"$7"\t"$8"\t"$9"\t"$10}' \
| awk -F'\t' 'BEGIN{
    while((getline < "/mnt/spareHD_2/nu_287/q2_parallelism/q2_tables_gene_snp_pop/INTERSECT_genes.list")>0)
      keep[$1]=1
} NR==1 || keep[$2]{print}' \
> /mnt/spareHD_2/nu_287/q2_parallelism/q2_tables_gene_snp_pop/INTERSECT_gene_snp_pops.tsv

wc -l /mnt/spareHD_2/nu_287/q2_parallelism/q2_tables_gene_snp_pop/INTERSECT_gene_snp_pops.tsv
head /mnt/spareHD_2/nu_287/q2_parallelism/q2_tables_gene_snp_pop/INTERSECT_gene_snp_pops.tsv
603 /mnt/spareHD_2/nu_287/q2_parallelism/q2_tables_gene_snp_pop/INTERSECT_gene_snp_pops.tsv
region	gene	snp	chr	pos	maj_sign	concordance	mean_abs_delta	pops_nonzero
AK	ATP5F1B	chrI:27378199:ATP5F1B	chrI	27378199	+	0.875	0.058100375	FG,LB,LG,SR,TL,WB,WT,WK
BC	ATP5F1B	chrI:27378199:ATP5F1B	chrI	27378199	+	1	0.0553936666666667	THE,JOE,BEA,MUC,PYE,AMO,GOS,ROB,BOOT,ECHO,LAW,SWA
BC	ATP5F1B	chrI:27378204:ATP5F1B	chrI	27378204	+	1	0.013889	THE,JOE,BEA,MUC,PYE,AMO,GOS,ROB,BOOT,ECHO,LAW,SWA
BC	ATP5F1B	chrI:27378216:ATP5F1B	chrI	27378216	+	1	0.0438988333333334	THE,JOE,BEA,MUC,PYE,AMO,GOS,ROB,BOOT,ECHO,LAW,SWA
AK	ATP5F1B	chrI:27378252:ATP5F1B	chrI	27378252	+	1	0.110687555555556	FG,LB,LG,SR,SL,TL,WB,WT,WK
BC	ATP5F1B	chrI:27378252:ATP5F1B	chrI	27378252	+	0.666666666666667	0.0374903333333333	THE,JOE,BEA,MUC,PYE,AMO,GOS,ROB,BOOT,ECHO,LAW,SWA
AK	ATP5F1B	chrI:27378275:ATP5F1B	chrI	27378275	+	1	0.0231479999999999	FG,LB,LG,SR,SL,TL,WB,WT,WK
AK	ATP5F1B	chrI:27379018:ATP5F1B	chrI	27379018	+	0.555555555555556	0.144910333333333	FG,LB,LG,SR,SL,TL,WB,WT,WK
BC	ATP5F1B	chrI:27379018:ATP5F1B	chrI	27379018	+	0.666666666666667	0.216926	THE,JOE,BEA,MUC,PYE,AMO,GOS,ROB,BOOT,ECHO,LAW,SWA
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ zcat /mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes/perSNP_withinRegion_concordance.tsv.gz \
| awk -F'\t' 'NR==1{print "region\tgene\tsnp\tchr\tpos\tmaj_sign\tconcordance\tmean_abs_delta\tpops_nonzero"; next}
             $1=="AK"{g=toupper($4); sub(/\r$/,"",g); print $1"\t"g"\t"$5"\t"$2"\t"$3"\t"$7"\t"$8"\t"$9"\t"$10}' \
| awk -F'\t' 'BEGIN{
    while((getline < "/mnt/spareHD_2/nu_287/q2_parallelism/q2_tables_gene_snp_pop/AK_parallel_genes.list")>0)
      keep[$1]=1
} NR==1 || keep[$2]{print}' \
> /mnt/spareHD_2/nu_287/q2_parallelism/q2_tables_gene_snp_pop/AK_gene_snp_pops.tsv

wc -l /mnt/spareHD_2/nu_287/q2_parallelism/q2_tables_gene_snp_pop/AK_gene_snp_pops.tsv
head /mnt/spareHD_2/nu_287/q2_parallelism/q2_tables_gene_snp_pop/AK_gene_snp_pops.tsv
219 /mnt/spareHD_2/nu_287/q2_parallelism/q2_tables_gene_snp_pop/AK_gene_snp_pops.tsv
region	gene	snp	chr	pos	maj_sign	concordance	mean_abs_delta	pops_nonzero
AK	ATP5F1B	chrI:27378199:ATP5F1B	chrI	27378199	+	0.875	0.058100375	FG,LB,LG,SR,TL,WB,WT,WK
AK	ATP5F1B	chrI:27378252:ATP5F1B	chrI	27378252	+	1	0.110687555555556	FG,LB,LG,SR,SL,TL,WB,WT,WK
AK	ATP5F1B	chrI:27378275:ATP5F1B	chrI	27378275	+	1	0.0231479999999999	FG,LB,LG,SR,SL,TL,WB,WT,WK
AK	ATP5F1B	chrI:27379018:ATP5F1B	chrI	27379018	+	0.555555555555556	0.144910333333333	FG,LB,LG,SR,SL,TL,WB,WT,WK
AK	ATP5F1B	chrI:27379326:ATP5F1B	chrI	27379326	+	0.875	0.270958875	LB,LG,SR,SL,TL,WB,WT,WK
AK	ATP5F1B	chrI:27379389:ATP5F1B	chrI	27379389	-	0.666666666666667	0.143726777777778	FG,LB,LG,SR,SL,TL,WB,WT,WK
AK	ATP5F1B	chrI:27379410:ATP5F1B	chrI	27379410	-	0.625	0.2723375	FG,LB,LG,SR,TL,WB,WT,WK
AK	ATP5F1B	chrI:27379476:ATP5F1B	chrI	27379476	+	1	0.151396444444444	FG,LB,LG,SR,SL,TL,WB,WT,WK
AK	ATP5F1B	chrI:27379690:ATP5F1B	chrI	27379690	+	1	0.012346	FG,LB,SR,SL,TL,WT,WK
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ zcat /mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes/perSNP_withinRegion_concordance.tsv.gz \
| awk -F'\t' 'NR==1{print "region\tgene\tsnp\tchr\tpos\tmaj_sign\tconcordance\tmean_abs_delta\tpops_nonzero"; next}
             $1=="BC"{g=toupper($4); sub(/\r$/,"",g); print $1"\t"g"\t"$5"\t"$2"\t"$3"\t"$7"\t"$8"\t"$9"\t"$10}' \
| awk -F'\t' 'BEGIN{
    while((getline < "/mnt/spareHD_2/nu_287/q2_parallelism/q2_tables_gene_snp_pop/BC_parallel_genes.list")>0)
      keep[$1]=1
} NR==1 || keep[$2]{print}' \
> /mnt/spareHD_2/nu_287/q2_parallelism/q2_tables_gene_snp_pop/BC_gene_snp_pops.tsv

wc -l /mnt/spareHD_2/nu_287/q2_parallelism/q2_tables_gene_snp_pop/BC_gene_snp_pops.tsv
head /mnt/spareHD_2/nu_287/q2_parallelism/q2_tables_gene_snp_pop/BC_gene_snp_pops.tsv
576 /mnt/spareHD_2/nu_287/q2_parallelism/q2_tables_gene_snp_pop/BC_gene_snp_pops.tsv
region	gene	snp	chr	pos	maj_sign	concordance	mean_abs_delta	pops_nonzero
BC	ATP5F1B	chrI:27378199:ATP5F1B	chrI	27378199	+	1	0.0553936666666667	THE,JOE,BEA,MUC,PYE,AMO,GOS,ROB,BOOT,ECHO,LAW,SWA
BC	ATP5F1B	chrI:27378204:ATP5F1B	chrI	27378204	+	1	0.013889	THE,JOE,BEA,MUC,PYE,AMO,GOS,ROB,BOOT,ECHO,LAW,SWA
BC	ATP5F1B	chrI:27378216:ATP5F1B	chrI	27378216	+	1	0.0438988333333334	THE,JOE,BEA,MUC,PYE,AMO,GOS,ROB,BOOT,ECHO,LAW,SWA
BC	ATP5F1B	chrI:27378252:ATP5F1B	chrI	27378252	+	0.666666666666667	0.0374903333333333	THE,JOE,BEA,MUC,PYE,AMO,GOS,ROB,BOOT,ECHO,LAW,SWA
BC	ATP5F1B	chrI:27379018:ATP5F1B	chrI	27379018	+	0.666666666666667	0.216926	THE,JOE,BEA,MUC,PYE,AMO,GOS,ROB,BOOT,ECHO,LAW,SWA
BC	ATP5F1B	chrI:27379326:ATP5F1B	chrI	27379326	-	0.583333333333333	0.249932	THE,JOE,BEA,MUC,PYE,AMO,GOS,ROB,BOOT,ECHO,LAW,SWA
BC	ATP5F1B	chrI:27379389:ATP5F1B	chrI	27379389	+	0.727272727272727	0.309179545454545	THE,JOE,BEA,MUC,PYE,AMO,GOS,BOOT,ECHO,LAW,SWA
BC	ATP5F1B	chrI:27379410:ATP5F1B	chrI	27379410	+	0.583333333333333	0.242088416666667	THE,JOE,BEA,MUC,PYE,AMO,GOS,ROB,BOOT,ECHO,LAW,SWA
BC	ATP5F1B	chrI:27379437:ATP5F1B	chrI	27379437	+	1	0.027397	THE,JOE,BEA,MUC,PYE,AMO,GOS,ROB,BOOT,ECHO,LAW,SWA
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ IN="/mnt/spareHD_2/nu_287/q2_parallelism/q2_tables_gene_snp_pop/INTERSECT_gene_snp_pops.tsv"
OUT="/mnt/spareHD_2/nu_287/q2_parallelism/q2_tables_gene_snp_pop/INTERSECT_sharedSNP_AK_BC_twoLines.tsv"

awk -F'\t' '
NR==1{hdr=$0; next}
{
  key=$2 FS $3   # gene + snp
  seen[key][$1]=1
  line[key,$1]=$0
}
END{
  print hdr
  for(k in seen){
    # only keep if both AK and BC exist
    if(seen[k]["AK"] && seen[k]["BC"]){
      print line[k,"AK"]
      print line[k,"BC"]
    }
  }
' "$IN" | sort -t$'\t' -k2,2 -k3,3 -k1,1 > "$OUT"

wc -l "$OUT"
head "$OUT"
awk: cmd. line:16:   }
awk: cmd. line:16:    ^ unexpected newline or end of string
0 /mnt/spareHD_2/nu_287/q2_parallelism/q2_tables_gene_snp_pop/INTERSECT_sharedSNP_AK_BC_twoLines.tsv
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ IN="/mnt/spareHD_2/nu_287/q2_parallelism/q2_tables_gene_snp_pop/INTERSECT_gene_snp_pops.tsv"
OUT="/mnt/spareHD_2/nu_287/q2_parallelism/q2_tables_gene_snp_pop/INTERSECT_sharedSNP_AK_BC_twoLines.tsv"

awk -F'\t' '
NR==1{hdr=$0; next}
{
  key = $2 FS $3              # gene + snp
  if($1=="AK"){ hasAK[key]=1; lineAK[key]=$0 }
  else if($1=="BC"){ hasBC[key]=1; lineBC[key]=$0 }
}
END{
  print hdr
  for(key in hasAK){
    if(hasBC[key]){
      print lineAK[key]
      print lineBC[key]
    }
  }
}
' "$IN" | sort -t$'\t' -k2,2 -k3,3 -k1,1 > "$OUT"

wc -l "$OUT"
head "$OUT"
327 /mnt/spareHD_2/nu_287/q2_parallelism/q2_tables_gene_snp_pop/INTERSECT_sharedSNP_AK_BC_twoLines.tsv
AK	ATP5F1A	chrXIV:3623321:ATP5F1A	chrXIV	3623321	+	1	0.0562958888888889	FG,LB,LG,SR,SL,TL,WB,WT,WK
BC	ATP5F1A	chrXIV:3623321:ATP5F1A	chrXIV	3623321	+	1	0.037037	THE,JOE,BEA,MUC,PYE,AMO,GOS,ROB,BOOT,ECHO,LAW,SWA
AK	ATP5F1A	chrXIV:3623363:ATP5F1A	chrXIV	3623363	+	0.875	0.2071035	FG,LB,SR,SL,TL,WB,WT,WK
BC	ATP5F1A	chrXIV:3623363:ATP5F1A	chrXIV	3623363	+	0.75	0.240047333333333	THE,JOE,BEA,MUC,PYE,AMO,GOS,ROB,BOOT,ECHO,LAW,SWA
AK	ATP5F1A	chrXIV:3623507:ATP5F1A	chrXIV	3623507	+	0.777777777777778	0.151566888888889	FG,LB,LG,SR,SL,TL,WB,WT,WK
BC	ATP5F1A	chrXIV:3623507:ATP5F1A	chrXIV	3623507	+	0.75	0.111673416666667	THE,JOE,BEA,MUC,PYE,AMO,GOS,ROB,BOOT,ECHO,LAW,SWA
AK	ATP5F1A	chrXIV:3623588:ATP5F1A	chrXIV	3623588	+	0.888888888888889	0.173290222222222	FG,LB,LG,SR,SL,TL,WB,WT,WK
BC	ATP5F1A	chrXIV:3623588:ATP5F1A	chrXIV	3623588	+	0.75	0.241112	THE,JOE,BEA,MUC,PYE,AMO,GOS,ROB,BOOT,ECHO,LAW,SWA
AK	ATP5F1A	chrXIV:3623805:ATP5F1A	chrXIV	3623805	+	0.625	0.2320805	FG,LB,LG,SR,SL,TL,WB,WK
BC	ATP5F1A	chrXIV:3623805:ATP5F1A	chrXIV	3623805	-	0.818181818181818	0.286829909090909	THE,JOE,BEA,MUC,PYE,GOS,ROB,BOOT,ECHO,LAW,SWA
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ echo "shared SNP count =" $(( ( $(wc -l < "$OUT") - 1 ) / 2 ))
shared SNP count = 163
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ awk -F'\t' 'NR>1{print $2}' /mnt/spareHD_2/nu_287/q2_parallelism/q2_tables_gene_snp_pop/INTERSECT_sharedSNP_AK_BC_twoLines.tsv \
| sort -u | wc -l

awk -F'\t' 'NR>1{print $2}' /mnt/spareHD_2/nu_287/q2_parallelism/q2_tables_gene_snp_pop/INTERSECT_sharedSNP_AK_BC_twoLines.tsv \
| sort -u
14
ATP5F1A
ATP5F1B
CYC1
HCCSB
NDUFA10
NDUFS1
NDUFS2
NDUFS3
NDUFS8
SDHA
SDHC
UQCRC1
UQCRC2B
gene
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_withinRegion_72genes$ 












#表格 
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(VennDiagram) # 如果没有请先 install.packages("VennDiagram")
})

# ==========================================
# Path (使用你系统里的中间文件)
# ==========================================
IN_LM <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_mtCluster_manual_perm0/gene_enrichment/LM_with_q_and_driver.tsv.gz"
OUTDIR <- "/mnt/spareHD_2/nu_287/q2_parallelism/comparison_AK_BC"
dir.create(OUTDIR, showWarnings = FALSE)

# 阈值设置
Q_CUT <- 0.10

# ==========================================
# 1. 读取并筛选显著位点
# ==========================================
cat("[read] Loading LM results...\n")
LM <- fread(IN_LM)

# 提取 AK 和 BC 的显著 SNP 集合 (q < Q_CUT)
ak_sig_snps <- LM[region == "AK" & q_cluster < Q_CUT, snp]
bc_sig_snps <- LM[region == "BC" & q_cluster < Q_CUT, snp]

# 提取基因集合
ak_sig_genes <- unique(LM[region == "AK" & q_cluster < Q_CUT, gene])
bc_sig_genes <- unique(LM[region == "BC" & q_cluster < Q_CUT, gene])

# ==========================================
# 2. 计算重合度与超几何检验 (Hypergeometric Test)
# ==========================================
# 总体背景（Universe）：两个 region 都测到的基因
universe_genes <- intersect(unique(LM[region=="AK", gene]), unique(LM[region=="BC", gene]))
intersect_genes <- intersect(ak_sig_genes, bc_sig_genes)

# 超几何检验：在 AK 显著的基因中，有多少也是 BC 显著的，是否显著高于随机预期？
# phyper(q, m, n, k)
p_overlap <- 1 - phyper(length(intersect_genes) - 1, 
                         length(ak_sig_genes), 
                         length(universe_genes) - length(ak_sig_genes), 
                         length(bc_sig_genes))

cat("\n=== Gene-level Comparison ===\n")
cat("AK Sig Genes: ", length(ak_sig_genes), "\n")
cat("BC Sig Genes: ", length(bc_sig_genes), "\n")
cat("Shared Genes: ", length(intersect_genes), "\n")
cat("Hypergeometric p-value: ", p_overlap, "\n")

# ==========================================
# 3. 输出重合的基因列表
# ==========================================
shared_df <- LM[gene %in% intersect_genes & q_cluster < Q_CUT, 
                .(region, snp, gene, driver_cluster, q_cluster)]
fwrite(shared_df, file.path(OUTDIR, "Shared_Sig_Genes_AK_BC.tsv"), sep="\t")

# ==========================================
# 4. 可视化 (Venn Diagram)
# ==========================================
venn.plot <- venn.diagram(
  x = list(AK = ak_sig_genes, BC = bc_sig_genes),
  filename = file.path(OUTDIR, "Venn_AK_BC_Overlap.png"),
  disable.logging = TRUE,
  main = paste0("Overlap of Cluster-associated Genes (q < ", Q_CUT, ")"),
  fill = c("skyblue", "pink"),
  alpha = 0.5
)

cat("\n[OK] Results saved in: ", OUTDIR, "\n")



# ==========================================
# 分析 40 个共享基因的 Driver Cluster 一致性
# ==========================================
cat("[analysis] Comparing Driver Clusters for shared genes...\n")

# 1. 提取共享基因在两地的 Driver
shared_genes <- intersect(LM[region=="AK" & q_cluster < 0.1, gene], 
                          LM[region=="BC" & q_cluster < 0.1, gene])

drivers_comp <- LM[gene %in% shared_genes & q_cluster < 0.1, 
                   .(driver = unique(driver_cluster)), 
                   by = .(region, gene)]

# 2. 转换为宽表对比
comp_table <- dcast(drivers_comp, gene ~ region, value.var = "driver")

# 3. 统计驱动模式对 (AK_driver vs BC_driver)
driver_summary <- comp_table[, .N, by = .(AK, BC)][order(-N)]

cat("\n=== Driver Cluster Match (AK vs BC) ===\n")
print(driver_summary)

# 4. 可视化驱动者匹配矩阵
p_driver <- ggplot(driver_summary, aes(x = AK, y = BC, fill = N)) +
  geom_tile() +
  geom_text(aes(label = N), color = "white") +
  scale_fill_gradient(low = "steelblue", high = "red") +
  labs(title = "Consistency of Mitochondrial Drivers",
       subtitle = "Number of shared genes driven by specific Cluster pairs",
       x = "Driver Cluster in Alaska",
       y = "Driver Cluster in British Columbia") +
  theme_minimal()

ggsave(file.path(final_out, "Fig_Driver_Consistency_Matrix.png"), p_driver, width = 6, height = 5)

cat("[OK] Driver analysis saved to: ", final_out, "\n")

#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(VennDiagram) # 如果没有请先 install.packages("VennDiagram")
})

# ==========================================
# Path (使用你系统里的中间文件)
# ==========================================
IN_LM <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_mtCluster_manual_perm0/gene_enrichment/LM_with_q_and_driver.tsv.gz"
OUTDIR <- "/mnt/spareHD_2/nu_287/q2_parallelism/comparison_AK_BC"
dir.create(OUTDIR, showWarnings = FALSE)

# 阈值设置
Q_CUT <- 0.10

# ==========================================
# 1. 读取并筛选显著位点
# ==========================================
cat("[read] Loading LM results...\n")
LM <- fread(IN_LM)

# 提取 AK 和 BC 的显著 SNP 集合 (q < Q_CUT)
ak_sig_snps <- LM[region == "AK" & q_cluster < Q_CUT, snp]
bc_sig_snps <- LM[region == "BC" & q_cluster < Q_CUT, snp]

# 提取基因集合
ak_sig_genes <- unique(LM[region == "AK" & q_cluster < Q_CUT, gene])
bc_sig_genes <- unique(LM[region == "BC" & q_cluster < Q_CUT, gene])

# ==========================================
# 2. 计算重合度与超几何检验 (Hypergeometric Test)
# ==========================================
# 总体背景（Universe）：两个 region 都测到的基因
universe_genes <- intersect(unique(LM[region=="AK", gene]), unique(LM[region=="BC", gene]))
intersect_genes <- intersect(ak_sig_genes, bc_sig_genes)

# 超几何检验：在 AK 显著的基因中，有多少也是 BC 显著的，是否显著高于随机预期？
# phyper(q, m, n, k)
p_overlap <- 1 - phyper(length(intersect_genes) - 1, 
                         length(ak_sig_genes), 
                         length(universe_genes) - length(ak_sig_genes), 
                         length(bc_sig_genes))

cat("\n=== Gene-level Comparison ===\n")
cat("AK Sig Genes: ", length(ak_sig_genes), "\n")
cat("BC Sig Genes: ", length(bc_sig_genes), "\n")
cat("Shared Genes: ", length(intersect_genes), "\n")
cat("Hypergeometric p-value: ", p_overlap, "\n")

# ==========================================
# 3. 输出重合的基因列表
# ==========================================
shared_df <- LM[gene %in% intersect_genes & q_cluster < Q_CUT, 
                .(region, snp, gene, driver_cluster, q_cluster)]
fwrite(shared_df, file.path(OUTDIR, "Shared_Sig_Genes_AK_BC.tsv"), sep="\t")

# ==========================================
# 4. 可视化 (Venn Diagram)
# ==========================================
venn.plot <- venn.diagram(
  x = list(AK = ak_sig_genes, BC = bc_sig_genes),
  filename = file.path(OUTDIR, "Venn_AK_BC_Overlap.png"),
  disable.logging = TRUE,
  main = paste0("Overlap of Cluster-associated Genes (q < ", Q_CUT, ")"),
  fill = c("skyblue", "pink"),
  alpha = 0.5
)

cat("\n[OK] Results saved in: ", OUTDIR, "\n")

You said
=== Gene-level Comparison ===

AK Sig Genes:  69 

BC Sig Genes:  40 

Shared Genes:  40 

Hypergeometric p-value:  0.08316566 


# ==========================================
# 分析 40 个共享基因的 Driver Cluster 一致性
# ==========================================
cat("[analysis] Comparing Driver Clusters for shared genes...\n")

# 1. 提取共享基因在两地的 Driver
shared_genes <- intersect(LM[region=="AK" & q_cluster < 0.1, gene], 
                          LM[region=="BC" & q_cluster < 0.1, gene])

drivers_comp <- LM[gene %in% shared_genes & q_cluster < 0.1, 
                   .(driver = unique(driver_cluster)), 
                   by = .(region, gene)]

# 2. 转换为宽表对比
comp_table <- dcast(drivers_comp, gene ~ region, value.var = "driver")

# 3. 统计驱动模式对 (AK_driver vs BC_driver)
driver_summary <- comp_table[, .N, by = .(AK, BC)][order(-N)]

cat("\n=== Driver Cluster Match (AK vs BC) ===\n")
print(driver_summary)

# 4. 可视化驱动者匹配矩阵
p_driver <- ggplot(driver_summary, aes(x = AK, y = BC, fill = N)) +
  geom_tile() +
  geom_text(aes(label = N), color = "white") +
  scale_fill_gradient(low = "steelblue", high = "red") +
  labs(title = "Consistency of Mitochondrial Drivers",
       subtitle = "Number of shared genes driven by specific Cluster pairs",
       x = "Driver Cluster in Alaska",
       y = "Driver Cluster in British Columbia") +
  theme_minimal()

ggsave(file.path(final_out, "Fig_Driver_Consistency_Matrix.png"), p_driver, width = 6, height = 5)

cat("[OK] Driver analysis saved to: ", final_out, "\n")
=== Driver Cluster Match (AK vs BC) ===

      AK    BC     N

   <int> <int> <int>

1:     1     1    28

2:     2     1     5

3:     2     2     3

4:     1     2     2

5:     1     3     1

6:     3     1     1

[OK] Driver analysis saved to:  /mnt/spareHD_2/nu_287/q2_parallelism/comparison_AK_BC 
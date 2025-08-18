nuclear pi calculation
/work/cyu/nuOXPHOS_genes_final.bed 


# merge BAM files (in the order of the file paths in BAMlist.txt) in a MPILEUP file only retaining nucleotides with BQ >20 and reads with MQ > 20
samtools mpileup -B \
    -f /work/cyu/stickleback_nuclear_only.fa \
    -b bamlist.txt \
    -q 30 \
    -Q 30 \
    -d 5000 \
    | gzip > nuclear.mpileup.gz

gunzip /mnt/spareHD_2/marked_duplicates/nuclear.mpileup.gz

/mnt/spareHD_2/nuclear_marked_duplicates/nuclear.mpileup

perl mpileup2sync.pl \
    --input /mnt/spareHD_2/nuclear_marked_duplicates/crispr.mpileup \
    --output /mnt/spareHD_2/nuclear_marked_duplicates/nuclear.sync \
    --fastq-type sanger \
    --min-qual 20



#ann and create bed files

cd /work/cyu/poolseq/PPalign_output/ann/

grep -P "\tCDS\t" /home/cyu/snpEff/data/Gasterosteus_aculeatus_MT/genes.gff | \
awk '{print $1, $4, $5}' OFS='\t' > pcg.bed

grep -P "\tgene\t" /home/cyu/snpEff/data/Gasterosteus_aculeatus_MT/genes.gff | grep "rRNA" | \
awk '{print $1, $4, $5}' OFS='\t' > rrna.bed

grep -P "\tgene\t" /home/cyu/snpEff/data/Gasterosteus_aculeatus_MT/genes.gff | grep "tRNA" | \
awk '{print $1, $4, $5}' OFS='\t' > trna.bed


grep -P "\tCDS\t" /home/cyu/snpEff/data/Gasterosteus_aculeatus_MT/genes.gff | \
awk '{print $1, $4, $5}' OFS='\t' > cds.bed

grep -P "\tD[-_]loop\t|\tcontrol_region\t" /home/cyu/snpEff/data/Gasterosteus_aculeatus_MT/genes.gff | awk '{print $1, $4, $5}' OFS='\t' > dloop.bed



grep -P "\tCDS\t" /home/cyu/snpEff/data/Gasterosteus_aculeatus_MT/genes.gff | \
grep -w "CYTB" | \
awk '{print $1, $4, $5}' OFS='\t' > complex3.bed

sed -i 's/NC_041244.1/MH205729.1/g' rrna.bed
awk '{print $1, $2, $2+1, $0}' OFS='\t' /work/cyu/poolseq/PPalign_output/mtDNA_bam/fish.sync > fish.sync.bed



sed -i 's/NC_041244.1/MH205729.1/g' rrna.bed
bedtools intersect -a /work/cyu/poolseq/PPalign_output/ann/fish.sync.bed -b /work/cyu/poolseq/PPalign_output/ann/trna.bed > /work/cyu/poolseq/PPalign_output/ann/fish_trna.bed.sync
bedtools intersect -a fish.sync.bed -b ND6.bed > ND6.bed.sync


bedtools intersect -a fish.sync.bed -b complex1.bed > complex1.bed.sync

awk '{print $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, $30, $31, $32, $33}' OFS='\t' complex1.bed.sync > complex1_fixed.sync
awk '{print $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, $30, $31, $32, $33}' OFS='\t' fish_rrna.bed.sync > fish_rrna_fixed.sync





grenedalf diversity \
  --sync-path /work/cyu/poolseq/PPalign_output/diversity/mt_genes_fixed.sync \
  --rename-samples-list /work/cyu/poolseq/PPalign_output/fst/rename.txt \
  --filter-sample-min-count 2 \
  --filter-sample-min-read-depth 4 \
  --window-type regions \
  --window-region-bed /work/cyu/poolseq/PPalign_output/diversity/mt_genes_new.bed \
  --window-average-policy valid-loci \
  --pool-sizes /work/cyu/poolseq/PPalign_output/fst/poolsize.txt \
  --out-dir /work/cyu/poolseq/PPalign_output/diversity \
  --allow-file-overwriting 


  bedtools intersect -a /work/cyu/poolseq/PPalign_output/ann/fish.sync.bed -b /work/cyu/poolseq/PPalign_output/diversity/mt_genes_new.bed > mt_genes.sync
  awk '{print $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, $30, $31, $32, $33}' OFS='\t' mt_genes.sync > mt_genes_fixed.sync

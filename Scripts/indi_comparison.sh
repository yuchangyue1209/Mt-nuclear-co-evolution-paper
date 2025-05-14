#1_extract_mt_reads
#!/bin/bash
# Extract mitochondrial reads from BAM files
cd /mnt/spareHD_2/AK_bam/Threespine_AK_source_lake_BAMs/
MT_CONTIG="NC_041244.1_mitochondion_genome"
OUTDIR="./mito_bams"
mkdir -p "$OUTDIR"

for BAM in *.bam; do
    BASENAME=$(basename "$BAM" .bam)
    OUTBAM="$OUTDIR/${BASENAME}_mt.bam"
    echo "Extracting mtDNA from: $BAM"
    samtools view -b "$BAM" "$MT_CONTIG" > "$OUTBAM"
    samtools index "$OUTBAM"
    echo "✅ Finished: $OUTBAM"
done

echo "🎉 mtDNA extraction complete for all BAMs."

#2_deduplicate_mt_bams
#!/bin/bash
# Mark duplicates in mitochondrial BAMs

export JAVA_HOME=/home/cyu/jdk-17.0.12
export PATH=$JAVA_HOME/bin:$PATH

INPUT_DIR="/mnt/spareHD_2/AK_bam/Threespine_AK_source_lake_BAMs/mito_bams"
OUTPUT_DIR="${INPUT_DIR}_dedup"
mkdir -p "$OUTPUT_DIR"

for BAM_FILE in "$INPUT_DIR"/*_mt.bam; do
    BASENAME=$(basename "$BAM_FILE" .bam)
    OUTPUT_BAM="${OUTPUT_DIR}/${BASENAME}_dedup.bam"
    METRICS_FILE="${OUTPUT_DIR}/${BASENAME}_metrics.txt"
    LOG_FILE="${OUTPUT_DIR}/${BASENAME}_picard.log"

    echo "Deduplicating: $BASENAME"
    picard MarkDuplicates \
        I="$BAM_FILE" \
        O="$OUTPUT_BAM" \
        M="$METRICS_FILE" \
        REMOVE_DUPLICATES=true \
        ASSUME_SORTED=true \
        VALIDATION_STRINGENCY=LENIENT > "$LOG_FILE" 2>&1

    if [[ -f "$OUTPUT_BAM" ]]; then
        samtools index "$OUTPUT_BAM"
        echo "✅ Done: $OUTPUT_BAM"
    else
        echo "❌ Error: $OUTPUT_BAM not created" | tee -a "$OUTPUT_DIR/error.log"
    fi
done

echo "🎉 Deduplication complete."

#3_variant_call_consensus
#!/bin/bash
# Variant calling and consensus generation (haploid)

INPUT_DIR="/mnt/spareHD_2/AK_bam/Threespine_AK_source_lake_BAMs/mito_bams_dedup"
VCF_DIR="/mnt/spareHD_2/AK_bam/vcf_individual_mt"
CONS_DIR="/mnt/spareHD_2/AK_bam/consensus_individual_mt"
REFERENCE="/work/cyu/sequence.fasta"

mkdir -p "$VCF_DIR" "$CONS_DIR"

for BAM in "$INPUT_DIR"/*_dedup.bam; do
    SAMPLE=$(basename "$BAM" _dedup.bam)
    echo "Calling variants: $SAMPLE"

    bcftools mpileup -Ou -f "$REFERENCE" -d 5000 "$BAM" | \
    bcftools call --ploidy 1 -mv -Oz -o "$VCF_DIR/${SAMPLE}_mt.vcf.gz"

    bcftools index "$VCF_DIR/${SAMPLE}_mt.vcf.gz"

    bcftools consensus -f "$REFERENCE" -M N "$VCF_DIR/${SAMPLE}_mt.vcf.gz" > \
        "$CONS_DIR/${SAMPLE}_mt_consensus.fasta"

    awk -v name="$SAMPLE" '/^>/{print ">" name; next} {print}' "$CONS_DIR/${SAMPLE}_mt_consensus.fasta" > \
        "$CONS_DIR/${SAMPLE}_mt_consensus.fasta.tmp" && \
        mv "$CONS_DIR/${SAMPLE}_mt_consensus.fasta.tmp" "$CONS_DIR/${SAMPLE}_mt_consensus.fasta"

    echo "✅ Consensus complete: $SAMPLE"
done

echo "🎉 All consensus FASTAs generated."


#4_rename_chr_header
#!/bin/bash
# Rename mitochondrial contig in BAM headers

INPUT_DIR="/mnt/spareHD_2/AK_bam/Threespine_AK_source_lake_BAMs/mito_bams_dedup"
OUTPUT_DIR="${INPUT_DIR}/renamed"
OLD_CHR="NC_041244.1_mitochondion_genome"
NEW_CHR="MH205729.1"

mkdir -p "$OUTPUT_DIR"

for BAM_FILE in "$INPUT_DIR"/*_dedup.bam; do
    SAMPLE=$(basename "$BAM_FILE" .bam)
    echo "Renaming header for: $SAMPLE"

    samtools view -H "$BAM_FILE" > "$OUTPUT_DIR/${SAMPLE}.header.sam"
    sed "s/${OLD_CHR}/${NEW_CHR}/g" "$OUTPUT_DIR/${SAMPLE}.header.sam" > "$OUTPUT_DIR/${SAMPLE}.newheader.sam"
    samtools reheader "$OUTPUT_DIR/${SAMPLE}.newheader.sam" "$BAM_FILE" > "$OUTPUT_DIR/${SAMPLE}_renamed.bam"
    samtools index "$OUTPUT_DIR/${SAMPLE}_renamed.bam"

    echo "✅ Header updated: $SAMPLE"
done

echo "🎉 All headers renamed and reindexed."



#5_variant_call_from_renamed
#!/bin/bash
# Variant calling and consensus generation from renamed BAMs

INPUT_DIR="/mnt/spareHD_2/AK_bam/Threespine_AK_source_lake_BAMs/mito_bams_dedup/renamed"
VCF_DIR="/mnt/spareHD_2/AK_bam/vcf_individual_mt"
CONS_DIR="/mnt/spareHD_2/AK_bam/consensus_individual_mt"
REFERENCE="/work/cyu/sequence.fasta"

mkdir -p "$VCF_DIR" "$CONS_DIR"

for BAM in "$INPUT_DIR"/*_renamed.bam; do
    SAMPLE=$(basename "$BAM" _renamed.bam)
    echo "Processing renamed BAM: $SAMPLE"

    bcftools mpileup -Ou -f "$REFERENCE" -d 5000 "$BAM" | \
    bcftools call --ploidy 1 -mv -Oz -o "$VCF_DIR/${SAMPLE}_mt.vcf.gz"

    bcftools index "$VCF_DIR/${SAMPLE}_mt.vcf.gz"

    bcftools consensus -f "$REFERENCE" -M N "$VCF_DIR/${SAMPLE}_mt.vcf.gz" > \
        "$CONS_DIR/${SAMPLE}_mt_consensus.fasta"

    awk -v name="$SAMPLE" '/^>/{print ">" name; next} {print}' "$CONS_DIR/${SAMPLE}_mt_consensus.fasta" > \
        "$CONS_DIR/${SAMPLE}_mt_consensus.fasta.tmp" && \
        mv "$CONS_DIR/${SAMPLE}_mt_consensus.fasta.tmp" "$CONS_DIR/${SAMPLE}_mt_consensus.fasta"

    echo "✅ Finished: $SAMPLE"
done

echo "🎉 All renamed BAMs processed for consensus."


#6_align_and_build_tree
#!/bin/bash
# Multiple sequence alignment and phylogenetic tree construction

cd /mnt/spareHD_2/AK_bam/consensus_individual_mt

cat *_consensus.fasta > all_consensus_merged.fasta
cat *_consensus.fasta > all_consensus_merged_allpoolseq.fasta
cat all_consensus_merged_allpoolseq.fasta JS_male_hap1.v20240517.chrM.fa > all_consensus_indi_japan.fasta
# Alignment
mafft --auto all_consensus_merged.fasta > aligned_WT.fasta
mafft --auto all_consensus_merged_allpoolseq.fasta > aligned_WT_merged_allpoolseq.fasta
mafft --auto all_consensus_indi_japan.fasta > aligned_all_consensus_indi_japan.fasta
# Tree inference
iqtree -s aligned_WT.fasta -m MFP -bb 1000 -nt AUTO -pre WT_mtDNA_tree
iqtree -s aligned_WT_merged_allpoolseq.fasta -m MFP -bb 1000 -nt AUTO -pre WT_allpoolseq_mtDNA_tree
iqtree -s aligned_all_consensus_indi_japan.fasta -m MFP -bb 1000 -nt AUTO -pre WT_allpoolseq_japan_mtDNA_tree


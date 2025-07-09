#nuclear data process script
#Step1 trimming
#!/bin/bash

# path
INPUT_DIR="/mnt/spareHD_2/rawdata"
TRIMMED_DIR="/mnt/spareHD_2/trimmed_poolseq"
FASTQC_DIR="/mnt/spareHD_2/fastqc_reports_trimmed"
ADAPTER_REF="/home/cyu/.conda/envs/poolseq_env/opt/bbmap-39.01-1/resources/adapters.fa"
THREADS=48

# make directory
mkdir -p "$TRIMMED_DIR"
mkdir -p "$FASTQC_DIR"

# trimming for
for R1_FILE in "$INPUT_DIR"/*_R1_001.fastq.gz; do
    R2_FILE="${R1_FILE/_R1_/_R2_}"
    PREFIX=$(basename "$R1_FILE" | sed 's/_R1_001.fastq.gz//')

    # output
    OUT_R1="$TRIMMED_DIR/trimmed_R1_${PREFIX}.fastq.gz"
    OUT_R2="$TRIMMED_DIR/trimmed_R2_${PREFIX}.fastq.gz"
    LOG_FILE="$TRIMMED_DIR/${PREFIX}_bbduk_log.txt"

    echo ">>> Trimming $PREFIX ..."
    bbduk.sh \
        in1="$R1_FILE" \
        in2="$R2_FILE" \
        out1="$OUT_R1" \
        out2="$OUT_R2" \
        ref="$ADAPTER_REF" \
        ktrim=rl \
        trimq=20 \
        minlength=25 \
        ftl=10 \
        tossbrokenreads=t \
        threads="$THREADS" > "$LOG_FILE" 2>&1

    echo ">>> Running FastQC on trimmed files for $PREFIX ..."
    fastqc "$OUT_R1" "$OUT_R2" \
        --outdir="$FASTQC_DIR" \
        --threads="$THREADS"
done

echo "✅ All trimming and FastQC complete."



#Step2 mapping
#!/bin/bash

# ---- Configurable Paths ----
TRIMMED_DIR="/mnt/spareHD_2/trimmed_poolseq"
REFERENCE="/work/cyu/stickleback_nuclear_only.fa"
INDEX_PREFIX="/work/cyu/stickleback_nuclear_only"
SAM_DIR="/mnt/spareHD_2/nuclear_mapped_poolseq"
BAM_DIR="/mnt/spareHD_2/nuclear_sorted_bam_poolseq"
THREADS=48


# ---- Create necessary output directories ----
mkdir -p "$SAM_DIR"
mkdir -p "$BAM_DIR"

# ---- Step 1: Build Bowtie2 index if missing ----
if [[ ! -f "${INDEX_PREFIX}.1.bt2" ]]; then
    echo "📌 Building Bowtie2 index..."
    bowtie2-build "$REFERENCE" "$INDEX_PREFIX"
fi

# ---- Step 2: Loop through all paired reads and map ----
for R1_FILE in "$TRIMMED_DIR"/trimmed_R1_*.fastq.gz; do
    BASENAME=$(basename "$R1_FILE" | sed 's/trimmed_R1_//; s/.fastq.gz//')
    R2_FILE="$TRIMMED_DIR/trimmed_R2_${BASENAME}.fastq.gz"

    if [[ -f "$R2_FILE" ]]; then
        SAM_FILE="$SAM_DIR/${BASENAME}.sam"
        RAW_BAM="$BAM_DIR/${BASENAME}.bam"
        SORTED_BAM="$BAM_DIR/${BASENAME}_sorted.bam"
        LOG_FILE="$SAM_DIR/${BASENAME}_bowtie2.log"

        echo "🔄 Mapping $BASENAME..."
        bowtie2 -x "$INDEX_PREFIX" \
            -1 "$R1_FILE" \
            -2 "$R2_FILE" \
            -p "$THREADS" \
            --no-mixed \
            --no-discordant \
            -X 2000 \
            -S "$SAM_FILE" > "$LOG_FILE" 2>&1

        echo "📦 Converting $SAM_FILE to BAM (MAPQ ≥ 20)..."
        samtools view -b -q 20 "$SAM_FILE" > "$RAW_BAM"

        echo "🔃 Sorting $RAW_BAM..."
        samtools sort -o "$SORTED_BAM" "$RAW_BAM"

        echo "🔍 Indexing $SORTED_BAM..."
        samtools index "$SORTED_BAM"

        echo "🧹 Cleaning up intermediate files..."
        rm "$SAM_FILE" "$RAW_BAM"

        echo "✅ Finished processing $BASENAME"
    else
        echo "⚠️  Warning: Missing R2 for $BASENAME" | tee -a "$BAM_DIR/mapping_warnings.log"
    fi
done

echo "🎉 All samples mapped, filtered, sorted, and indexed."



#Step3 dedup
#dedup
export JAVA_HOME=/home/cyu/jdk-17.0.12
export PATH=$JAVA_HOME/bin:$PATH

#!/bin/bash

# Directories
INPUT_DIR="/mnt/spareHD_2/nuclear_sorted_bam_poolseq"
OUTPUT_DIR="/mnt/spareHD_2/nuclear_marked_duplicates"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all sorted BAM files
for BAM_FILE in "$INPUT_DIR"/*_sorted.bam; do
    BASENAME=$(basename "$BAM_FILE" _sorted.bam)
    OUTPUT_BAM="${OUTPUT_DIR}/${BASENAME}_dedup.bam"
    METRICS_FILE="${OUTPUT_DIR}/${BASENAME}_metrics.txt"
    LOG_FILE="${OUTPUT_DIR}/${BASENAME}_picard.log"

    echo "🧹 Marking duplicates for $BASENAME..."

    picard MarkDuplicates \
        I="$BAM_FILE" \
        O="$OUTPUT_BAM" \
        M="$METRICS_FILE" \
        REMOVE_DUPLICATES=true \
        ASSUME_SORTED=true \
        VALIDATION_STRINGENCY=LENIENT > "$LOG_FILE" 2>&1

    echo "✅ Done: $OUTPUT_BAM"
done

echo "🎉 All BAM files deduplicated."


#step4 add groupnames with picard

#!/bin/bash

# Step 4: Add or replace read groups
INPUT_DIR="/mnt/spareHD_2/nuclear_marked_duplicates"
OUTPUT_DIR="/mnt/spareHD_2/nuclear_with_readgroups"
mkdir -p "$OUTPUT_DIR"

for BAM_FILE in "$INPUT_DIR"/*_dedup.bam; do
    BASENAME=$(basename "$BAM_FILE" _dedup.bam)    # e.g. 3_SR_S3
    PREFIX=$(echo "$BASENAME" | cut -d '_' -f1,2)   # e.g. 3_SR

    OUTPUT_BAM="${OUTPUT_DIR}/${PREFIX}_rg.bam"

    echo "📌 Adding RG for $PREFIX..."

    java -Xmx16g -jar /home/cyu/picard/picard.jar AddOrReplaceReadGroups \
        I="$BAM_FILE" \
        O="$OUTPUT_BAM" \
        RGID="$PREFIX" \
        RGLB="$PREFIX" \
        RGPL=ILLUMINA \
        RGPU="${PREFIX}_Unit1" \
        RGSM="$PREFIX" \
        VALIDATION_STRINGENCY=LENIENT

    echo "✅ RG added: $OUTPUT_BAM"
done

echo "🎯 Step 4 done: All BAM files have read groups added."

#!/bin/bash

# Step 5: Index BAM files
INPUT_DIR="/mnt/spareHD_2/nuclear_with_readgroups"

for BAM_FILE in "$INPUT_DIR"/*_rg.bam; do
    echo "🧬 Indexing $(basename "$BAM_FILE")..."
    samtools index "$BAM_FILE"
    echo "✅ Indexed: $(basename "$BAM_FILE")"
done

echo "🎉 Step 5 done: All BAM files indexed."



#Step 6   gtak 4.6.1.0
#!/bin/bash

# path
POOL_INFO="/work/cyu/poolseq/pool_info.txt"
INPUT_DIR="/mnt/spareHD_2/nuclear_with_readgroups"   
OUTPUT_DIR="/mnt/spareHD_2/gatk4_vcf"        
REFERENCE="/work/cyu/stickleback_nuclear_only.fa"                            

mkdir -p "$OUTPUT_DIR"

# except the header
tail -n +2 "$POOL_INFO" | while read -r SAMPLE SIZE; do
    # only numbers 
    CLEAN_SIZE=$(echo "$SIZE" | tr -cd '0-9')
    
    if [ -z "$CLEAN_SIZE" ]; then
        echo "Warning: For sample $SAMPLE, extracted pool size is empty. Skipping."
        continue
    fi

    # path
    BAM_FILE="${INPUT_DIR}/${SAMPLE}_rg.bam"
    if [ ! -f "$BAM_FILE" ]; then
        echo "Warning: BAM file not found for sample $SAMPLE at $BAM_FILE. Skipping."
        continue
    fi

    echo "Processing sample $SAMPLE with pool size $CLEAN_SIZE..."
    
    #  GATK HaplotypeCaller，generate GVCF files for populations
    gatk HaplotypeCaller \
      -R "$REFERENCE" \
      -I "$BAM_FILE" \
      -O "${OUTPUT_DIR}/${SAMPLE}.g.vcf.gz" \
      -ERC GVCF \
      --ploidy "$CLEAN_SIZE" \
      --minimum-mapping-quality 20 \
      -mbq 13 \
      --indel-size-to-eliminate-in-ref-model 12 \
      -G AS_StandardAnnotation \
      -G StandardAnnotation \
      --sample-name "$SAMPLE"

    if [ $? -ne 0 ]; then
        echo "Error processing sample $SAMPLE" >&2
        continue
    fi

    echo "Finished processing sample $SAMPLE"
done

echo "All samples processed."

#genotype 
#!/usr/bin/env bash

REFERENCE="/work/cyu/sequence.fasta"
GVCF_DIR="/work/cyu/poolseq/PPalign_output/gatk4_vcf"
SAMPLE_MAP="${GVCF_DIR}/sample_map_fixed.txt"
GENOMICSDB_WORKSPACE="${GVCF_DIR}/genomicsdb_workspace"
FINAL_VCF="${GVCF_DIR}/all_samples_genotyped.vcf.gz"
INTERVAL="MH205729.1"

echo "=== Step 1: GenomicsDBImport ==="
gatk GenomicsDBImport \
  --sample-name-map "${SAMPLE_MAP}" \
  --genomicsdb-workspace-path "${GENOMICSDB_WORKSPACE}" \
  -L "${INTERVAL}" \
  --genomicsdb-vcf-buffer-size 1048576 \
  --tmp-dir /tmp \
  --reader-threads 4

if [ $? -ne 0 ]; then
  echo "[ERROR] GenomicsDBImport failed!"
  exit 1
fi

echo "=== Step 2: GenotypeGVCFs ==="
gatk GenotypeGVCFs \
  -R "${REFERENCE}" \
  -V "gendb://${GENOMICSDB_WORKSPACE}" \
  -O "${FINAL_VCF}" \
  --tmp-dir /tmp

echo "All done. Final VCF: ${FINAL_VCF}"




1.#create bam.list
ls -1 /mnt/spareHD_2/nuclear_marked_duplicates/*.bam > /mnt/spareHD_2/nuclear_marked_duplicates/bamlist.txt

# merge BAM files (in the order of the file paths in BAMlist.txt) in a MPILEUP file only retaining nucleotides with BQ >20 and reads with MQ > 20
samtools mpileup -B \
    -f /work/cyu/stickleback_nuclear_only.fa \
    -b bamlist.txt \
    -q 30 \
    -Q 30 \
    -d 5000 \
    | gzip > crispr.mpileup.gz

gunzip /mnt/spareHD_2/marked_duplicates/crispr.mpileup.gz



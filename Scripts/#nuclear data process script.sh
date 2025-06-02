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


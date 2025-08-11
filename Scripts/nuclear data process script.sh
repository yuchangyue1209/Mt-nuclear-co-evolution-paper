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






#subset genes bam files
#!/bin/bash

BED="/work/cyu/nuOXPHOS_genes_with_complex_core.bed"
INPUT_DIR="/mnt/spareHD_2/nuclear_with_readgroups"
OUTPUT_DIR="/mnt/spareHD_2/nuclear_with_readgroups_subset"

# 创建输出目录（如果不存在）
mkdir -p "$OUTPUT_DIR"

# 遍历所有 BAM 文件
for BAM_FILE in "$INPUT_DIR"/*_rg.bam; do
    SAMPLE=$(basename "$BAM_FILE" _rg.bam)
    
    echo "✂️ Subsetting $SAMPLE..."

    bedtools intersect -abam "$BAM_FILE" -b "$BED" > "$OUTPUT_DIR/${SAMPLE}_subset.bam"

    samtools index "$OUTPUT_DIR/${SAMPLE}_subset.bam"

    echo "✅ Finished: ${SAMPLE}_subset.bam"
done

echo "🎉 All BAMs subset to OXPHOS regions!"


#use bcftools
#!/bin/bash

REFERENCE="/work/cyu/stickleback_nuclear_only.fa"
VCF_DIR="/mnt/spareHD_2/nuclear_oxphos_vcf/"
SUBSET_DIR="/mnt/spareHD_2/nuclear_with_readgroups_subset/"

mkdir -p "$VCF_DIR"

for f in "$SUBSET_DIR"/*_subset.bam; do
  sample=$(basename "$f" "_subset.bam")

  echo "Calling variants for (diploid mode): $sample"

  # 1. mpileup call variants
  bcftools mpileup -Ou -f "$REFERENCE" -d 5000 \
    --annotate FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP \
    "$f" | \
  bcftools call -mv -Oz -o "$VCF_DIR/${sample}_oxphos_raw.vcf.gz"

  # 2. VCF names
  bcftools reheader -s <(echo "$sample") "$VCF_DIR/${sample}_oxphos_raw.vcf.gz" > \
    "$VCF_DIR/${sample}_oxphos_temp.vcf.gz"

  # 3. Normalization 
  bcftools norm -m -any -f "$REFERENCE" -Oz -o "$VCF_DIR/${sample}_oxphos.vcf.gz" \
    "$VCF_DIR/${sample}_oxphos_temp.vcf.gz"

  rm "$VCF_DIR/${sample}_oxphos_raw.vcf.gz"
  rm "$VCF_DIR/${sample}_oxphos_temp.vcf.gz"

  # 4. Index VCF
  bcftools index "$VCF_DIR/${sample}_oxphos.vcf.gz"
done


#!/bin/bash

VCF_DIR="/mnt/spareHD_2/nuclear_oxphos_vcf/"
CONS_DIR="/mnt/spareHD_2/neclear_consensus_oxphos"
REF="/work/cyu/stickleback_nuclear_only.fa"

mkdir -p "$CONS_DIR"

for f in "$VCF_DIR"/*_oxphos.vcf.gz; do
  sample=$(basename "$f" "_oxphos.vcf.gz")

  echo "Generating diploid consensus sequence for: $sample"

  bcftools consensus -f "$REF" -M N --sample "$sample" "$f" > \
    "$CONS_DIR/${sample}_oxphos_consensus.fasta"
done








1.#create bam.list
ls -1 /mnt/spareHD_2/nuclear_marked_duplicates/*.bam > /mnt/spareHD_2/nuclear_marked_duplicates/bamlist.txt

# merge BAM files (in the order of the file paths in BAMlist.txt) in a MPILEUP file only retaining nucleotides with BQ >20 and reads with MQ > 20
samtools mpileup -B \
    -f /work/cyu/stickleback_nuclear_only.fa \
    -b bamlist.txt \
    -q 30 \
    -Q 30 \
    -d 5000 \
    | gzip > nuclear.mpileup.gz

gunzip /mnt/spareHD_2/nuclear_marked_duplicates/crispr.mpileup.gz /mnt/spareHD_2/nuclear_marked_duplicates/nuclear.mpileup





/mnt/spareHD_2/oxphos_gene_tree/bamlist_nuclear.txt

#!/usr/bin/env bash
set -euo pipefail

REF=/work/cyu/stickleback_nuclear_only.fa.gz
SUB=/mnt/spareHD_2/nuclear_with_readgroups_subset
BEDDIR=/mnt/spareHD_2/oxphos_gene_tree/beds
OUT=/mnt/spareHD_2/oxphos_gene_tree
mkdir -p "$OUT/mpileup"

for B in "$BEDDIR"/*.bed; do
  GENE=$(basename "$B" .bed)
  echo "🧬 mpileup $GENE"

  samtools mpileup \
    -aa \
    -Q 20 \
    -B  \
    -d 500000 \
    -f "$REF" \
    -l "$B" \
    "$SUB"/*_subset.bam \
    > "$OUT/mpileup/${GENE}.mpileup"
done

echo -e "\n🎉 mpileup 文件写入：$OUT/mpileup/"

samtools mpileup -B \
    -f /work/cyu/stickleback_nuclear_only.fa \
    -b bamlist.txt \
    -q 30 \
    -Q 30 \
    -d 5000 \
    | gzip > nuclear.mpileup.gz

gunzip /mnt/spareHD_2/marked_duplicates/nuclear.mpileup.gz
/work/cyu/poolseq/PPalign_output/mtDNA_bam/bamlist.txt

#!/usr/bin/env bash
set -euo pipefail                      # 出错即停，未定义变量即报错

#####################  路径配置  #####################
REF=/work/cyu/stickleback_nuclear_only.fa.gz                # 参考基因组（已 bgzip+faidx）
MT_BAMLIST=/work/cyu/poolseq/PPalign_output/mtDNA_bam/bamlist.txt
NUC_DIR=/mnt/spareHD_2/nuclear_with_readgroups_subset       # *_subset.bam
BEDDIR=/mnt/spareHD_2/oxphos_gene_tree/beds                 # 每基因 .bed
OUT=/mnt/spareHD_2/oxphos_gene_tree                         # 输出目录
THREADS=8                                                   # 当前 samtools 版本不支持 -@，仅留作记录
############################################################

mkdir -p "$OUT"/{mpileup,beds}

#####################  1. 生成核 bamlist，保持与 mtDNA 顺序一致  #####################
BAMLIST_NUC=$OUT/bamlist_nuclear.txt
awk -v d="$NUC_DIR" '
  {
    match($0, /([0-9]+_[A-Z]+)/, a);        # 提取样本 ID（如 10_THE）
    file=d"/"a[1]"_subset.bam";
    if (system("[ -f "file" ]")==0) print file;    # 核 BAM 存在才写
    else {
        printf("⚠️  Warning: %s 不存在，已跳过\n", file) > "/dev/stderr"
    }
  }' "$MT_BAMLIST" > "$BAMLIST_NUC"

echo "📄  核 bamlist 写入 $(wc -l < "$BAMLIST_NUC") 条路径 → $BAMLIST_NUC"

#####################  2. 逐基因 mpileup（参数与示例保持一致）  #####################
for BED in "$BEDDIR"/*.bed; do
    GENE=$(basename "$BED" .bed)
    echo "🧬  mpileup  $GENE"

    samtools mpileup \
        -B \
        -f "$REF" \
        -b "$BAMLIST_NUC" \
        -q 30 \
        -Q 30 \
        -d 5000 \
        -l "$BED" \
      | gzip > "$OUT/mpileup/${GENE}.mpileup.gz"
done

echo -e "\n🎉  mpileup.gz 文件全部写入：$OUT/mpileup/"

cd /mnt/spareHD_2/oxphos_gene_tree/mpileup
gunzip *.mpileup.gz           # 串行
# 或并行
ls *.mpileup.gz | xargs -n1 -P8 gunzip





MP_DIR=/mnt/spareHD_2/oxphos_gene_tree/mpileup
SYNC_DIR=/mnt/spareHD_2/oxphos_gene_tree/sync
mkdir -p "$SYNC_DIR"

for MP in "$MP_DIR"/*.mpileup; do
  GENE=$(basename "$MP" .mpileup)
  perl mpileup2sync.pl \
       --input "$MP" \
       --output "$SYNC_DIR/${GENE}.sync" \
       --fastq-type sanger \
       --min-qual 30
done


#!/usr/bin/env bash
set -euo pipefail

# -----------  路径 ------------
CF_DIR=/mnt/spareHD_2/oxphos_gene_tree/counts   # 你的 .cf 文件目录
TREE_DIR=/mnt/spareHD_2/oxphos_gene_tree/trees  # 输出目录
mkdir -p "$TREE_DIR"

THREADS=AUTO      # AUTO = IQ-TREE 自动检测；也可写 16、32 等
BOOT=1000         # Ultrafast bootstrap 轮数
ALRT=1000         # SH-aLRT 轮数

# -----------  循环跑 -----------
for CF in "$CF_DIR"/*.cf; do
  GENE=$(basename "$CF" .cf)
  echo "🌳  IQ-Tree  $GENE  ($(date +%T))"

  iqtree2 -s "$CF" \
          -m GTR+P \
          -nt "$THREADS" \
          -B "$BOOT" \
          -alrt "$ALRT" \
          -pre "$TREE_DIR/$GENE" \
          -quiet
done

echo -e "\n✅  所有基因树完成！树文件在：$TREE_DIR/*.treefile"

#构建master树 用nuclear oxpho 82genetree
java -Xmx8g -jar astral.5.7.8.jar \
  -i /mnt/spareHD_2/oxphos_gene_tree/genes_astral.tre \
  -o /mnt/spareHD_2/oxphos_gene_tree/species_astral.tre \
  2> /mnt/spareHD_2/oxphos_gene_tree/astral.log
astral.5.7.8.jar 
cat /mnt/spareHD_2/oxphos_gene_tree/species_astral.tre
(10_THE,(16_AMO,(((9_SWA,(11_JOE,12_BEA)0.49:0.052704063490568435)0.86:0.12852164098344246,((13_MUC,(22_FRED,24_PACH)1:0.7981770630522537)0.74:0.09167667908034977,((27_LB,(17_SAY,(28_CH,(25_RS,26_SC)1:0.4418327522790392)1:1.0091776080159818)0.68:0.16080153569716776)0.98:0.22450967855654755,((4_SL,(8_WK,7_WT)0.92:0.15588829972579915)0.44:0.02734415112952233,((5_TL,3_SR)0.64:0.06916681982473395,(2_LG,(1_FG,6_WB)0.98:0.20763936477824518)0.51:0.03889609052423365)0.49:0.04632609721427273)0.83:0.12022862571625097)0.5:0.03645182279797124)0.5:0.03708014865417043)0.44:0.022098788242902498,((14_PYE,19_ROB)0.66:0.07594115323538288,(23_LAW,((18_GOS,21_ECHO)0.97:0.2709576295105863,20_BOOT)1:0.8748478121554343)1:0.5341491814983622)0.38:0.007488021780721204)0.49:0.03425810992000614):0.0);



#!/usr/bin/env bash
set -euo pipefail

MASTER_TRE="/mnt/spareHD_2/oxphos_gene_tree/species_astral.tre"
CF_DIR="/mnt/spareHD_2/oxphos_gene_tree/counts"
OUT_DIR="/mnt/spareHD_2/oxphos_gene_tree/pruned_trees"
THREADS=8

mkdir -p "$OUT_DIR"

# 提取 master tree 物种名（去掉分支长度和支持值）
tr '(),:;' '\n' < "$MASTER_TRE" | grep -vE '^$|^[0-9.]+$' | sort -u > /tmp/master.tips

for CF in "$CF_DIR"/*.cf; do
    GENE=$(basename "$CF" .cf)
    echo "🌿 处理 $GENE ..."

    # 从 .cf 文件第2行提取物种名
    awk 'NR==2 && $1=="CHROM"{for(i=3;i<=NF;i++) print $i}' "$CF" | sort -u > /tmp/cf.taxa

    # 取交集
    comm -12 /tmp/master.tips /tmp/cf.taxa > /tmp/keep.taxa

    NUM_TAXA=$(wc -l < /tmp/keep.taxa)
    if (( NUM_TAXA < 4 )); then
        echo "⚠️  $GENE 物种数太少 ($NUM_TAXA)，跳过"
        continue
    fi

    # 用 IQ-TREE 固定 topology 重新估分支长度
    iqtree2 -s "$CF" \
            -m GTR+P \
            -g "$MASTER_TRE" \
            -nt "$THREADS" \
            --safe \
            -pre "$OUT_DIR/${GENE}_fixed" \
            -quiet \
    || { echo "❌ $GENE 运行失败，跳过"; continue; }

    echo "✅ $GENE 完成"
done

echo -e "\n🎯 所有基因树已剪枝并重新估分支长度，结果在 $OUT_DIR"

cd /mnt/spareHD_2/oxphos_gene_tree/pruned_trees

for f in $(ls *_fixed.treefile | sort); do
    gene=${f%%_*}   # 文件名前缀为 gene 名
    echo "[$gene]"  # 这是注释行
    cat "$f"
done > /mnt/spareHD_2/oxphos_gene_tree/genes_astral_named.tre

treefile = "/users/changyueyu/desktop/genes_astral_named.tre"

outputfile = "/users/changyueyu/desktop/out.RDS"
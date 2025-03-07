A. #Trim raw_data with bbduk 

#!/bin/bash

# Set input and output directories and parameters
INPUT_DIR="/work/cyu/poolseq/raw_data"
OUTPUT_DIR="/work/cyu/poolseq/PPalign_output/trimmed"
ADAPTER_REF="/home/cyu/adapters.fa"
THREADS=48

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all paired R1 and R2 files
for R1_FILE in "$INPUT_DIR"/*_R1_001.fastq; do
    # Get the corresponding R2 file
    R2_FILE="${R1_FILE/_R1_/_R2_}"
    
    # Extract the file prefix (e.g., 16_AMO_S16)
    PREFIX=$(basename "$R1_FILE" | sed 's/_R1_001.fastq//')
    
    # Set output file paths
    OUT_R1="$OUTPUT_DIR/trimmed_R1_${PREFIX}.fastq"
    OUT_R2="$OUTPUT_DIR/trimmed_R2_${PREFIX}.fastq"
    LOG_FILE="$OUTPUT_DIR/${PREFIX}_bbduk_log.txt"
    
    # Run bbduk.sh
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
done

2. Map trimmed data to reference mt genome with bowtie2

#!/bin/bash
# Set input and output paths
TRIMMED_DIR="/work/cyu/poolseq/PPalign_output/trimmed"
MAPPED_DIR="/work/cyu/poolseq/PPalign_output//mapped"
REFERENCE="/work/cyu/chrM_index"

# Create output directory (if it does not exist)
mkdir -p "$MAPPED_DIR"

# Loop through all R1 files and find the corresponding R2 files
for R1_FILE in "$TRIMMED_DIR"/trimmed_R1_*.fastq; do
    # Get the base name of the file (remove the path and prefix)
    BASENAME=$(basename "$R1_FILE" | sed 's/trimmed_R1_//; s/.fastq//')
    R2_FILE="$TRIMMED_DIR/trimmed_R2_$BASENAME.fastq"
    
    # Check if the R2 file exists
    if [[ -f "$R2_FILE" ]]; then
        # Output file paths
        SAM_FILE="$MAPPED_DIR/${BASENAME}.sam"
        LOG_FILE="$MAPPED_DIR/${BASENAME}_bowtie2.log"
        
        # Perform mapping and save logs
        echo "Mapping $BASENAME..."
        bowtie2 -x "$REFERENCE" \
            -1 "$R1_FILE" \
            -2 "$R2_FILE" \
            -p 48 \
            --very-sensitive-local \
            --no-mixed \
            --no-discordant \
            -X 2000 \
            -S "$SAM_FILE" > "$LOG_FILE" 2>&1
    else
        echo "Warning: R2 file for $BASENAME not found, skipping." | tee -a "$MAPPED_DIR/mapping_warnings.log"
    fi
done

echo "Mapping completed."

3. Convert mapped sam files to bam and index
#!/bin/bash

# Define input and output directories
INPUT_DIR="/work/cyu/poolseq/PPalign_output/mapped"
DEPTH_OUTPUT_DIR="/work/cyu/poolseq/PPalign_output/depth"

# Create output directory if it doesn't exist
mkdir -p "$DEPTH_OUTPUT_DIR"

# Process all sorted BAM files
for BAM_FILE in "${INPUT_DIR}"/*_sorted.bam; do
    # Extract the basename (remove directory and extension)
    BASENAME=$(basename "$BAM_FILE" _sorted.bam)
    
    # Define output file path
    DEPTH_FILE="${DEPTH_OUTPUT_DIR}/${BASENAME}_depth.txt"
    
    echo "Calculating depth for $BAM_FILE..."
    
    # Compute depth using samtools
    samtools depth "$BAM_FILE" > "$DEPTH_FILE"
    
    echo "Depth file created: $DEPTH_FILE"
done

echo "All sorted BAM files have been processed for depth calculation."


4. Check the coverage/depth

#!/bin/bash

# Input and output directories
INPUT_DIR="/work/cyu/poolseq/PPalign_output/mapped"
OUTPUT_DIR="/work/cyu/poolseq/PPalign_output/mapped"

# Coverage output directory
COVERAGE_DIR="/work/cyu/poolseq/PPalign_output/coverage"
mkdir -p "$COVERAGE_DIR"

# Process all sorted BAM files
for BAM_FILE in "${OUTPUT_DIR}"/*_sorted.bam; do
    # Get the basename (removing directory and extension)
    BASENAME=$(basename "$BAM_FILE" _sorted.bam)
    
    # Coverage output file
    COVERAGE_FILE="${COVERAGE_DIR}/${BASENAME}_coverage.txt"
    
    echo "Calculating depth for $BAM_FILE..."
    
    # Calculate depth
    samtools depth "$BAM_FILE" > "$COVERAGE_FILE"
    
    echo "Coverage file created: $COVERAGE_FILE"
done

echo "All sorted BAM files have been processed for depth calculation."
6. Piscard add group name
#!/bin/bash

# Directories
INPUT_DIR="/work/cyu/poolseq/PPalign_output/mapped"
OUTPUT_DIR="/work/cyu/poolseq/PPalign_output/mapped_with_rg"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all marked BAM files (after dedup BAM 文件)
for BAM_FILE in "${INPUT_DIR}"/*_sorted.bam; do
    # Get the base name without extension
    BASENAME=$(basename "$BAM_FILE" _sorted.bam)
    
    # Output BAM file path
    OUTPUT_BAM="${OUTPUT_DIR}/${BASENAME}_rg.bam"
    
    # Run Picard AddOrReplaceReadGroups to add read groups
    echo "Processing $BAM_FILE..."
    java -Xmx16g -jar /home/cyu/picard/picard.jar AddOrReplaceReadGroups \
        -I "$BAM_FILE" \
        -O "$OUTPUT_BAM" \
        -RGID "$BASENAME" \
        -RGLB "$BASENAME" \
        -RGPL ILLUMINA \
        -RGPU "${BASENAME}_Unit1" \
        -RGSM "$BASENAME"
    
    echo "Finished processing $BAM_FILE: Read groups added."
done

echo "All BAM files have been processed."


7. Index group name added files

#!/bin/bash

# Input directory
INPUT_DIR="/work/cyu/poolseq/PPalign_output/mapped_with_rg"  

# Loop through all _rg.bam files
for BAM_FILE in "${INPUT_DIR}"/*_rg.bam; do
    # Get the sample name (remove the extension)
    SAMPLE_NAME=$(basename "$BAM_FILE" _rg.bam)
    
    # Index the BAM file
    echo "Indexing BAM file for ${SAMPLE_NAME}..."
    samtools index "$BAM_FILE"  # Create the BAM file index
    
    echo "Finished indexing ${SAMPLE_NAME}"
done

echo "All BAM files have been indexed."


8.gtak 4.6.1.0
#!/bin/bash

# path
POOL_INFO="/work/cyu/poolseq/pool_info.txt"
INPUT_DIR="/work/cyu/poolseq/PPalign_output/mapped_with_rg"   
OUTPUT_DIR="/work/cyu/poolseq/PPalign_output/gatk4_vcf"        
REFERENCE="/work/cyu/sequence.fasta"                            

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


9. # run bcftools ---2nd pipeline
mkdir -p /work/cyu/poolseq/PPalign_output/mtDNA_bam  

MTDNA_ID="MH205729.1"  # right name

for f in /work/cyu/poolseq/PPalign_output/mapped/*_sorted.bam
do
  sample=$(basename "$f" "_sorted.bam")
  
  echo "Extracting mtDNA from: $sample"

  # BAM 
  samtools view -b -q 20 -o /work/cyu/poolseq/PPalign_output/mtDNA_bam/${sample}_mtDNA.bam "$f" $MTDNA_ID

  # index bam
  samtools index /work/cyu/poolseq/PPalign_output/mtDNA_bam/${sample}_mtDNA.bam
done



mkdir -p /work/cyu/poolseq/PPalign_output/vcf

for f in /work/cyu/poolseq/PPalign_output/mtDNA_bam/*_mtDNA.bam
do
  sample=$(basename "$f" "_mtDNA.bam")

  echo "Calling variants for (haploid mode): $sample"
  # 1. mpileup call variants
  bcftools mpileup -Ou -f "$REFERENCE" -d 5000 \
    --annotate FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP \
    "$f" | \
  bcftools call --ploidy 1 -mv -Oz -o "$VCF_DIR/${sample}_mtDNA_raw.vcf.gz"
  # 2. VCF names
  bcftools reheader -s <(echo "$sample") /work/cyu/poolseq/PPalign_output/vcf/${sample}_mtDNA_raw.vcf.gz > \
  /work/cyu/poolseq/PPalign_output/vcf/${sample}_mtDNA_temp.vcf.gz

  # 3. Normalization 
  bcftools norm -m -any -f /work/cyu/sequence.fasta -Oz -o /work/cyu/poolseq/PPalign_output/vcf/${sample}_mtDNA.vcf.gz \
  /work/cyu/poolseq/PPalign_output/vcf/${sample}_mtDNA_temp.vcf.gz

  rm /work/cyu/poolseq/PPalign_output/vcf/${sample}_mtDNA_raw.vcf.gz
  rm /work/cyu/poolseq/PPalign_output/vcf/${sample}_mtDNA_temp.vcf.gz

  # 4. Index VCF
  bcftools index /work/cyu/poolseq/PPalign_output/vcf/${sample}_mtDNA.vcf.gz
done

mkdir -p /work/cyu/poolseq/PPalign_output/direct_consensus

REF_MT="/work/cyu/sequence.fasta"

for f in /work/cyu/poolseq/PPalign_output/vcf/*_mtDNA.vcf.gz
do
  sample=$(basename "$f" "_mtDNA.vcf.gz")

  echo "Generating haploid consensus sequence for: $sample"

  bcftools consensus -f "$REF_MT" -M N --sample "$sample" "$f" > \
  /work/cyu/poolseq/PPalign_output/direct_consensus/${sample}_mtDNA_consensus.fasta
done


for f in /work/cyu/poolseq/PPalign_output/direct_consensus/*_mtDNA_consensus.fasta
do
  sample=$(basename "$f" "_mtDNA_consensus.fasta")
  pop_name=$(echo "$sample" | cut -d'_' -f2)  

  echo "Renaming header for: $sample → $pop_name"

  awk -v name="$pop_name" '/^>/{print ">" name; next} {print}' "$f" > "${f}.tmp" && mv "${f}.tmp" "$f"
done

10. Intersection

#!/usr/bin/env bash

# ============== Configuration ==============
BCF_VCF_DIR="/work/cyu/poolseq/PPalign_output/vcf"             # Directory containing bcftools VCF files
GATK_VCF_DIR="/work/cyu/poolseq/PPalign_output/gatk4_vcf/single_sample_vcf"  # Directory containing GATK VCF files
OUT_DIR="/work/cyu/poolseq/PPalign_output/compare_results"     # Output directory for intersection results

mkdir -p "$OUT_DIR"

# ============== Execution ==============
# Iterate over *_mtDNA.vcf.gz files in the bcftools directory
for bcftools_vcf in "$BCF_VCF_DIR"/*_mtDNA.vcf.gz
do
    # If no matching files are found, skip
    [ -e "$bcftools_vcf" ] || { echo "No bcftools VCF found in $BCF_VCF_DIR"; exit 0; }

    # Extract sample name (remove _mtDNA.vcf.gz)
    fname=$(basename "$bcftools_vcf")   # e.g., 10_THE_mtDNA.vcf.gz
    sample="${fname%_mtDNA.vcf.gz}"     # e.g., 10_THE

    # Construct corresponding GATK VCF path
    gatk_vcf="${GATK_VCF_DIR}/${sample}_final.vcf.gz"

    # If the corresponding GATK file is not found, skip
    if [ ! -f "$gatk_vcf" ]; then
        echo "[WARNING] GATK VCF not found for sample $sample at $gatk_vcf, skipping..."
        continue
    fi

    # Define output directory: e.g., /.../compare_results/10_THE_isec
    outdir="${OUT_DIR}/${sample}_isec"
    mkdir -p "$outdir"

    echo "Running isec for sample: $sample"
    echo "  bcftools VCF: $bcftools_vcf"
    echo "  GATK VCF:     $gatk_vcf"
    echo "  Output Dir:   $outdir"

    # Execute bcftools isec 
    bcftools isec \
      "$bcftools_vcf" \
      "$gatk_vcf" \
      -p "$outdir" \
      -Oz

    echo "======================================================"
done

echo "All isec done. Check results in $OUT_DIR"

11. Filter and qc
#!/usr/bin/env bash
set -euo pipefail

# ============== Configuration ==============
COMPARE_OUT_DIR="/work/cyu/poolseq/PPalign_output/compare_results"
OVERLAP_DIR="/work/cyu/poolseq/PPalign_output/overlap.vcf"
mkdir -p "$OVERLAP_DIR"

# Filtering parameters (for single samples, using --minDP and --maxDP)
MIN_DP=10
MAX_DP=5000
MAX_MISSING=1
MINQ=30

# ============== Execution ==============
# Iterate over all *_isec directories in the compare results directory
for sample_dir in "$COMPARE_OUT_DIR"/*_isec; do
    # Extract sample name, e.g., 10_THE_isec -> 10_THE
    sample=$(basename "$sample_dir" _isec)
    
    # Define the intersection file path
    intersection_file="${sample_dir}/0002.vcf.gz"
    
    if [ -f "$intersection_file" ]; then
        echo "Processing sample: $sample"
        
        # Filter the intersection file and output it to OVERLAP_DIR with prefix <sample>_filtered
        vcftools --gzvcf "$intersection_file" \
          --minDP "$MIN_DP" \
          --maxDP "$MAX_DP" \
          --max-missing "$MAX_MISSING" \
          --minQ "$MINQ" \
          --recode --recode-INFO-all \
          --out "${OVERLAP_DIR}/${sample}_filtered"
        
        # The filtered file will be named ${sample}_filtered.recode.vcf
        filtered_vcf="${OVERLAP_DIR}/${sample}_filtered.recode.vcf"
        
        if [ -f "$filtered_vcf" ]; then
            # Compress the filtered VCF and rename it to <sample>_overlap.vcf.gz
            bgzip -c "$filtered_vcf" > "${OVERLAP_DIR}/${sample}_overlap.vcf.gz"
            bcftools index -f "${OVERLAP_DIR}/${sample}_overlap.vcf.gz"
            echo "  => Overlap VCF saved as: ${OVERLAP_DIR}/${sample}_overlap.vcf.gz"
            # Optionally delete intermediate files
            rm "$filtered_vcf"
            rm "${OVERLAP_DIR}/${sample}_filtered.log" 2>/dev/null || true
        else
            echo "Warning: Filtering did not produce an output for sample $sample."
        fi
    else
        echo "Warning: No intersection file (0002.vcf.gz) found for sample $sample"
    fi
    
    echo "======================================================"
done

echo "All overlap VCF files have been saved in: $OVERLAP_DIR"


12. header
#!/usr/bin/env bash
set -euo pipefail

# Define the directory containing consensus files
CONSENSUS_DIR="/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus"

# Iterate over all *_consensus.fasta files in the directory
for f in "$CONSENSUS_DIR"/*_consensus.fasta; do
    [ -e "$f" ] || { echo "No consensus files found in $CONSENSUS_DIR"; exit 0; }
    
    # Extract the sample name, e.g., "10_THE_consensus.fasta" -> "10_THE"
    sample=$(basename "$f" _consensus.fasta)
    # Extract the population name, using "_" as the delimiter and taking the second field, e.g., "10_THE" -> "THE"
    pop_name=$(echo "$sample" | cut -d'_' -f2)
    
    echo "Renaming header for sample: $sample → $pop_name"
    
    # Use awk to replace the FASTA header, only modifying the first line (starting with ">") to ">$pop_name"
    awk -v name="$pop_name" '/^>/{print ">" name; next} {print}' "$f" > "${f}.tmp" && mv "${f}.tmp" "$f"
done

echo "All consensus headers have been renamed."

# Concatenate all consensus FASTA files into a single file
cat /work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/*.fasta > /work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/all_mt_overlap.fasta

13. align
mafft --auto all_mt_overlap.fasta > aligned_mt_overlap.fasta

14. tree
iqtree -s aligned_mt_overlap.fasta -m GTR+G -nt AUTO

15. trans to nex FORMAT 
from Bio import SeqIO

def fasta_to_nexus(input_fasta, output_nexus):
    """
    Convert a FASTA file to NEXUS format for PopART, BEAST, MrBayes, and PAUP.
    Ensures all sequences have the same length.
    """

    # Read the FASTA file
    sequences = list(SeqIO.parse(input_fasta, "fasta"))
    
    # Ensure all sequences have the same length
    seq_lengths = {len(seq.seq) for seq in sequences}
    if len(seq_lengths) > 1:
        raise ValueError("Error: Sequences have different lengths. Alignment is required.")

    n_tax = len(sequences)  # Number of taxa
    n_char = len(sequences[0].seq)  # Sequence length

    # Write to NEXUS file
    with open(output_nexus, "w") as nexus:
        # NEXUS header
        nexus.write("#NEXUS\n\n")

        # TAXA block (defining samples)
        nexus.write("BEGIN TAXA;\n")
        nexus.write(f"  DIMENSIONS NTAX={n_tax};\n")
        nexus.write("  TAXLABELS\n")
        for seq in sequences:
            nexus.write(f"    {seq.id}\n")  # Sample name
        nexus.write("  ;\nEND;\n\n")

        # CHARACTERS block (defining sequences)
        nexus.write("BEGIN CHARACTERS;\n")
        nexus.write(f"  DIMENSIONS NCHAR={n_char};\n")
        nexus.write("  FORMAT DATATYPE=DNA MISSING=? GAP=-;\n")
        nexus.write("  MATRIX\n")
        for seq in sequences:
            nexus.write(f"    {seq.id} {str(seq.seq)}\n")  # Sample name + sequence
        nexus.write("  ;\nEND;\n")

# Set input & output file paths
input_fasta = "/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/aligned_mt_overlap.fasta"
output_nexus = "/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/aligned_mt_overlap.nex"

# Run conversion
fasta_to_nexus(input_fasta, output_nexus)

print(f"✅ Nexus file created: {output_nexus}")



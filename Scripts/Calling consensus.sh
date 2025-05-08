1. #Trim raw_data with bbduk 38.9

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

for r1 in /work/cyu/poolseq/PPalign_output/trimmed/trimmed_R1_*.fastq; do
    r2="${r1/R1/R2}"
    fastqc "$r1" "$r2" -o /work/cyu/poolseq/PPalign_output/quality_after_trim/
done


2. Map trimmed data to reference mt genome with bowtie2.5.4

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

3. Convert mapped sam files to bam files and index them
#!/bin/bash

# path
SAM_DIR="/work/cyu/poolseq/PPalign_output/mapped"
BAM_DIR="/work/cyu/poolseq/PPalign_output/mapped"

# sam
for SAM_FILE in "$SAM_DIR"/*.sam; do
    # name
    BASENAME=$(basename "$SAM_FILE" .sam)
    BAM_FILE="$BAM_DIR/${BASENAME}.bam"
    SORTED_BAM_FILE="$BAM_DIR/${BASENAME}_sorted.bam"
    # MAPQ ≥ 20 
    if [[ -f "$SAM_FILE" ]]; then
        echo "Converting $SAM_FILE to BAM with MAPQ ≥ 20..."
        samtools view -b -q 20 "$SAM_FILE" > "$BAM_FILE"

        echo "Sorting $BAM_FILE..."
        samtools sort -o "$SORTED_BAM_FILE" "$BAM_FILE"


        echo "Finished processing $BASENAME."
    else
        echo "Warning: No SAM files found in $SAM_DIR."
    fi
done

echo "Conversion and sorting completed."



4. Check the mapping quality: depth for each site
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


5. Piscard mark and de duplicates java17 I didn't run since data type--- mt is special！But for nuclear poolseq analysis, this step is requied
#!/bin/bash

# Directories
INPUT_DIR="/work/cyu/poolseq/PPalign_output/mapped"
OUTPUT_DIR="/work/cyu/poolseq/PPalign_output/marked_duplicates"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all sorted BAM files
for BAM_FILE in "${INPUT_DIR}"/*_sorted.bam; do
    # Get the base name without extension
    BASENAME=$(basename "$BAM_FILE" _sorted.bam)
    
    # Output BAM file path for the marked duplicates
    OUTPUT_BAM="${OUTPUT_DIR}/${BASENAME}_marked_duplicates.bam"
    
    # Run Picard MarkDuplicates
    echo "Processing $BAM_FILE..."
    picard MarkDuplicates \
        I="$BAM_FILE" \
        O="$OUTPUT_BAM" \
        M="${OUTPUT_DIR}/${BASENAME}_metrics.txt" \
        REMOVE_DUPLICATES=true
    
    echo "Finished processing $BAM_FILE: Duplicates marked."
done

echo "All BAM files have been processed"



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


9. # Run bcftools ---this is 2nd pipeline for variants calling 
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



10. After running these two standard pipeline, intersect the results
_mtDNA_run2.vcf.gz
#!/usr/bin/env bash

BCF_VCF_DIR="/work/cyu/poolseq/PPalign_output/vcf"             # bcftools vcf
GATK_VCF_DIR="/work/cyu/poolseq/PPalign_output/gatk4_vcf/single_sample_vcf"  # GATK vcf
OUT_DIR="/work/cyu/poolseq/PPalign_output/compare_results"     # output path

mkdir -p "$OUT_DIR"

#  bcftools *_mtDNA.vcf.gz 
for bcftools_vcf in "$BCF_VCF_DIR"/*_mtDNA.vcf.gz
do
    # exit if no vcf
    [ -e "$bcftools_vcf" ] || { echo "No bcftools VCF found in $BCF_VCF_DIR"; exit 0; }

    # extract sample name ( _mtDNA.vcf.gz)
    fname=$(basename "$bcftools_vcf")   #  10_THE_mtDNA.vcf.gz
    sample="${fname%_mtDNA.vcf.gz}"     #  10_THE

    # build path
    gatk_vcf="${GATK_VCF_DIR}/${sample}_final.vcf.gz"

    # exit if no gatk vcf
    if [ ! -f "$gatk_vcf" ]; then
        echo "[WARNING] GATK VCF not found for sample $sample at $gatk_vcf, skipping..."
        continue
    fi

    # output path： e.g. /.../compare_results/10_THE_isec
    outdir="${OUT_DIR}/${sample}_isec"
    mkdir -p "$outdir"

    echo "Running isec for sample: $sample"
    echo "  bcftools VCF: $bcftools_vcf"
    echo "  GATK VCF:     $gatk_vcf"
    echo "  Output Dir:   $outdir"

    # bcftools isec 
    bcftools isec \
      "$bcftools_vcf" \
      "$gatk_vcf" \
      -p "$outdir" \
      -Oz

    echo "======================================================"
done

echo "All isec done. Check results in $OUT_DIR"



11. Filter and QC
#!/usr/bin/env bash
set -euo pipefail

COMPARE_OUT_DIR="/work/cyu/poolseq/PPalign_output/compare_results"
OVERLAP_DIR="/work/cyu/poolseq/PPalign_output/overlap.vcf"
mkdir -p "$OVERLAP_DIR"

# filter parameters --minDP &--maxDP
MIN_DP=10
MAX_DP=5000
MAX_MISSING=1
MINQ=30

# for all compare  *_isec 
for sample_dir in "$COMPARE_OUT_DIR"/*_isec; do
    # extract sample name 10_THE_isec -> 10_THE
    sample=$(basename "$sample_dir" _isec)
    
    # intersection path defined
    intersection_file="${sample_dir}/0002.vcf.gz"
    
    if [ -f "$intersection_file" ]; then
        echo "Processing sample: $sample"
        
        # filter intersection files to output path *_filtered
        vcftools --gzvcf "$intersection_file" \
          --minDP "$MIN_DP" \
          --maxDP "$MAX_DP" \
          --max-missing "$MAX_MISSING" \
          --minQ "$MINQ" \
          --recode --recode-INFO-all \
          --out "${OVERLAP_DIR}/${sample}_filtered"
        
        # vcftools filtered file: ${sample}_filtered.recode.vcf
        filtered_vcf="${OVERLAP_DIR}/${sample}_filtered.recode.vcf"
        
        if [ -f "$filtered_vcf" ]; then
            # rename <sample>_overlap.vcf.gz
            bgzip -c "$filtered_vcf" > "${OVERLAP_DIR}/${sample}_overlap.vcf.gz"
            bcftools index -f "${OVERLAP_DIR}/${sample}_overlap.vcf.gz"
            echo "  => Overlap VCF saved as: ${OVERLAP_DIR}/${sample}_overlap.vcf.gz"
            # delete 
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



#indel check：no indels
bcftools merge --force-samples /work/cyu/poolseq/PPalign_output/overlap.vcf/*.vcf.gz -Oz -o /work/cyu/poolseq/PPalign_output/indel_regions/all_samples_merged.vcf.gz
tabix -p vcf /work/cyu/poolseq/PPalign_output/indel_regions/all_samples_merged.vcf.gz
(poolseq_env) cyu@stickleback:~$ bcftools view -v indels -Oz -o /work/cyu/poolseq/PPalign_output/indel_regions/all_samples_indels.vcf.gz \
/work/cyu/poolseq/PPalign_output/indel_regions/all_samples_merged.vcf.gz
(poolseq_env) cyu@stickleback:~$ bcftools query -l /work/cyu/poolseq/PPalign_output/indel_regions/all_samples_indels.vcf.gz



12. Call consensus for intersection snps from 2 pipelines
#!/usr/bin/env bash
set -euo pipefail

# path
OVERLAP_DIR="/work/cyu/poolseq/PPalign_output/overlap.vcf"
REFERENCE="/work/cyu/sequence.fasta"
CONSENSUS_DIR="${OVERLAP_DIR}/consensus"
mkdir -p "$CONSENSUS_DIR"

# for all *_overlap.vcf.gz 
for vcf in "$OVERLAP_DIR"/*_overlap.vcf.gz; do
    [ -e "$vcf" ] || { echo "No overlap VCF files found in $OVERLAP_DIR"; exit 0; }
    sample=$(basename "$vcf" _overlap.vcf.gz)
    echo "Generating consensus for sample: $sample"
  # bcftools consensus calling
    bcftools consensus -f "$REFERENCE" "$vcf" > "${CONSENSUS_DIR}/${sample}_consensus.fasta"
    
    echo "  => Consensus sequence saved to ${CONSENSUS_DIR}/${sample}_consensus.fasta"
done

echo "All consensus sequences have been generated in: $CONSENSUS_DIR"




12. Change header and combine the fasta files 
#!/usr/bin/env bash
set -euo pipefail

CONSENSUS_DIR="/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus"

for f in "$CONSENSUS_DIR"/*_consensus.fasta; do
    [ -e "$f" ] || { echo "No consensus files found in $CONSENSUS_DIR"; exit 0; }
    
    # extract name "10_THE_consensus.fasta" -> "10_THE"
    sample=$(basename "$f" _consensus.fasta)
    # "10_THE" -> "THE"
    pop_name=$(echo "$sample" | cut -d'_' -f2)
    
    echo "Renaming header for sample: $sample → $pop_name"
    
    # replace FASTA header to ">$pop_name"
    awk -v name="$pop_name" '/^>/{print ">" name; next} {print}' "$f" > "${f}.tmp" && mv "${f}.tmp" "$f"
done

echo "All consensus headers have been renamed."

#Merge the files into one fasta
cat /work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/*.fasta > /work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/all_mt_overlap.fasta

13. Realignment
mafft --auto all_mt_overlap.fasta > aligned_mt_overlap.fasta
muscle -in all_mt_overlap.fasta -out mc_aligned_mt_overlap.fasta

14. IQtree 
iqtree -s aligned_mt_overlap.fasta -m GTR+G -bb 1000 -alrt 1000 -nt AUTO
iqtree -s mc_aligned_mt_overlap.fasta -m GTR+G -bb 1000 -alrt 1000 -nt AUTO

15. Covert nexus format for popart
#nexus
from Bio import SeqIO

def fasta_to_nexus(input_fasta, output_nexus):
    """
    Convert a FASTA file to NEXUS format for PopART, BEAST, MrBayes, and PAUP.
    Ensures all sequences have the same length.
    """

    # fasta
    sequences = list(SeqIO.parse(input_fasta, "fasta"))
    
    # same length
    seq_lengths = {len(seq.seq) for seq in sequences}
    if len(seq_lengths) > 1:
        raise ValueError("Error: Sequences have different lengths. Alignment is required.")

    n_tax = len(sequences)  
    n_char = len(sequences[0].seq)  

    #  NEXUS 
    with open(output_nexus, "w") as nexus:
        # NEXUS head
        nexus.write("#NEXUS\n\n")

        # TAXA 
        nexus.write("BEGIN TAXA;\n")
        nexus.write(f"  DIMENSIONS NTAX={n_tax};\n")
        nexus.write("  TAXLABELS\n")
        for seq in sequences:
            nexus.write(f"    {seq.id}\n")  # sample name
        nexus.write("  ;\nEND;\n\n")

        # CHARACTERS 
        nexus.write("BEGIN CHARACTERS;\n")
        nexus.write(f"  DIMENSIONS NCHAR={n_char};\n")
        nexus.write("  FORMAT DATATYPE=DNA MISSING=? GAP=-;\n")
        nexus.write("  MATRIX\n")
        for seq in sequences:
            nexus.write(f"    {seq.id} {str(seq.seq)}\n")  
        nexus.write("  ;\nEND;\n")

# path
input_fasta = "/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/aligned_mt_overlap.fasta"
output_nexus = "/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/aligned_mt_overlap.nex"


fasta_to_nexus(input_fasta, output_nexus)

print(f"✅ Nexus file created: {output_nexus}")


16. Stat for SNPs in 4 steps
#!/usr/bin/env bash
set -euo pipefail

#path
BCF_VCF_DIR="/work/cyu/poolseq/PPalign_output/vcf"             # bcftools 
GATK_VCF_DIR="/work/cyu/poolseq/PPalign_output/gatk4_vcf/single_sample_vcf"  # GATK 
OUT_DIR="/work/cyu/poolseq/PPalign_output/compare_results"             # output intersection
FILTER_DIR="/work/cyu/poolseq/PPalign_output/overlap.vcf"              # after filtering

# output files for stat
OUTPUT_FILE="${OUT_DIR}/snp_stats.txt"

# header
echo -e "Sample\tGATK_SNP_count\tBCF_SNP_count\tOverlap_SNP_count\tAfterfilter_SNP_count" > "$OUTPUT_FILE"

# for all *_mtDNA_run2.vcf.gz
for bcf_vcf in "$BCF_VCF_DIR"/*_mtDNA.vcf.gz; do
    # extract name "10_THE_mtDNA_run2.vcf.gz" → "10_THE"
    sample=$(basename "$bcf_vcf" _mtDNA.vcf.gz)
    
    # stat bcftools SNP number（withour header）
    bcf_count=$(zcat "$bcf_vcf" | grep -v "^#" | wc -l)
    
    #  GATK VCF path "<sample>_final.vcf.gz"
    gatk_vcf="${GATK_VCF_DIR}/${sample}_final.vcf.gz"
    if [ ! -f "$gatk_vcf" ]; then
        echo "Warning: GATK VCF not found for sample $sample at $gatk_vcf, skipping..."
        continue
    fi
    gatk_count=$(zcat "$gatk_vcf" | grep -v "^#" | wc -l)
    
    # output）
    sample_isec_dir="${OUT_DIR}/${sample}_isec"
    mkdir -p "$sample_isec_dir"
    
    # bcftools isec  (0002.vcf.gz)
    if [ ! -f "${sample_isec_dir}/0002.vcf.gz" ]; then
        bcftools isec \
          "$bcf_vcf" \
          "$gatk_vcf" \
          -p "$sample_isec_dir" \
          -Oz
    fi
    
    #  SNP number
    if [ -f "${sample_isec_dir}/0002.vcf.gz" ]; then
        overlap_count=$(zcat "${sample_isec_dir}/0002.vcf.gz" | grep -v "^#" | wc -l)
    else
        overlap_count=0
    fi
    
    # FILTER_DIR  VCF  SNP number
    # path "<sample>_overlap.vcf.gz"
    filtered_vcf="${FILTER_DIR}/${sample}_overlap.vcf.gz"
    if [ -f "$filtered_vcf" ]; then
        filtered_count=$(zcat "$filtered_vcf" | grep -v "^#" | wc -l)
    else
        filtered_count=0
    fi

    # output
    echo -e "${sample}\t${gatk_count}\t${bcf_count}\t${overlap_count}\t${filtered_count}" >> "$OUTPUT_FILE"
done

echo "All sample comparisons done. Results have been saved in: $OUTPUT_FILE"



17. Check numts contamination
/work/cyu/poolparty/stickleback_v5_assembly_no_chrM.fa

makeblastdb -in stickleback_v5_assembly_no_chrM.fa -dbtype nucl -out nuclear_genome_no_chrM
mkdir -p /work/cyu/poolseq/PPalign_output/blast_results

for f in /work/cyu/poolseq/PPalign_output/direct_consensus/*_mtDNA_consensus.fasta
do
  sample=$(basename "$f" "_mtDNA_consensus.fasta")
  
  echo "Running BLAST for: $sample"

  blastn -query "$f" -db nuclear_genome_no_chrM -outfmt 6 \
  -out /work/cyu/poolseq/PPalign_output/blast_results/${sample}_blast.txt
done
/work/cyu/poolseq/PPalign_output/blast_results/




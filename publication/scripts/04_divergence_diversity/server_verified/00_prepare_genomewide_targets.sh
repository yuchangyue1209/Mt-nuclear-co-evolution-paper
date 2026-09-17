#!/usr/bin/env bash
set -euo pipefail

CLASSIFICATION="/path/to/workspace/Kuster2026_genomewide_reanalysis/stickleback_mapping/stickleback_gene_classification.final.tsv"
GFF3="/path/to/workspace/oxphos_from_ref_no_biomart/04_gff/stickleback_with_contact.gff3"
REF="/path/to/workspace/stickleback_nuclear_only.fa"

OUT="/path/to/data/genomewide_codeml_kuster/00_targets"
mkdir -p "$OUT"

GENES="$OUT/selected_genes.txt"
TRANSCRIPTS="$OUT/selected_transcripts.txt"
TX2GENE="$OUT/selected_transcript_to_gene.tsv"

CANON_GFF3="$OUT/stickleback_20426_canonical.gff3"
CANON_GTF="$OUT/stickleback_20426_canonical.gtf"

CDS_BY_GENE="$OUT/stickleback_20426_canonical_cds_by_gene.bed"
CDS_MERGED="$OUT/stickleback_20426_canonical_cds_merged.bed"

echo "[00] Extracting selected genes and transcripts"

awk -F'\t' '
NR>1 {
    print $1
}' "$CLASSIFICATION" |
sort -u > "$GENES"

awk -F'\t' '
NR>1 {
    print $2
}' "$CLASSIFICATION" |
sort -u > "$TRANSCRIPTS"

awk -F'\t' '
BEGIN {OFS="\t"}
NR>1 {
    print $2,$1
}' "$CLASSIFICATION" |
sort -u > "$TX2GENE"

echo "Selected genes:       $(wc -l < "$GENES")"
echo "Selected transcripts: $(wc -l < "$TRANSCRIPTS")"

if [[ $(wc -l < "$GENES") -ne 20426 ]]; then
    echo "ERROR: expected 20426 genes" >&2
    exit 1
fi

if [[ $(wc -l < "$TRANSCRIPTS") -ne 20426 ]]; then
    echo "ERROR: expected 20426 representative transcripts" >&2
    exit 1
fi

echo "[00] Filtering GFF3"

awk -F'\t' \
    -v genes_file="$GENES" \
    -v tx_file="$TRANSCRIPTS" '
BEGIN {
    OFS="\t"

    while ((getline x < genes_file) > 0)
        keep_gene[x]=1
    close(genes_file)

    while ((getline x < tx_file) > 0)
        keep_tx[x]=1
    close(tx_file)

    print "##gff-version 3"
}

$0 ~ /^#/ {
    next
}

{
    feature=$3
    attributes=$9
    id=""
    parent=""

    n=split(attributes,a,";")
    for (i=1; i<=n; i++) {
        if (a[i] ~ /^ID=/) {
            id=a[i]
            sub(/^ID=/,"",id)
        }
        else if (a[i] ~ /^Parent=/) {
            parent=a[i]
            sub(/^Parent=/,"",parent)
        }
    }

    keep=0

    if (feature=="gene" && id in keep_gene) {
        keep=1
    }
    else if ((feature=="mRNA" || feature=="transcript") && (id in keep_tx)) {
        keep=1
    }
    else if (parent!="") {
        m=split(parent,p,",")
        for (j=1; j<=m; j++) {
            if (p[j] in keep_tx) {
                keep=1
                break
            }
        }
    }

    if (keep)
        print
}
' "$GFF3" > "$CANON_GFF3"

echo "[00] Converting canonical GFF3 to GTF"

gffread \
    -T \
    -o "$CANON_GTF" \
    "$CANON_GFF3"

echo "[00] Creating gene-labelled CDS BED"

awk -F'\t' \
    -v tx2gene_file="$TX2GENE" '
BEGIN {
    OFS="\t"

    while ((getline x < tx2gene_file) > 0) {
        split(x,a,"\t")
        tx2gene[a[1]]=a[2]
    }
    close(tx2gene_file)
}

$0 !~ /^#/ && $3=="CDS" {
    parent=""

    n=split($9,a,";")
    for (i=1; i<=n; i++) {
        if (a[i] ~ /^Parent=/) {
            parent=a[i]
            sub(/^Parent=/,"",parent)
        }
    }

    parent_count=split(parent,p,",")
    for (j=1; j<=parent_count; j++) {
        tx=p[j]

        if (tx in tx2gene) {
            gene=tx2gene[tx]

            # BED: 0-based start, half-open end
            print $1,$4-1,$5,gene,0,$7,tx
        }
    }
}
' "$CANON_GFF3" > "$CDS_BY_GENE"

echo "[00] Sorting CDS BED"

if [[ ! -f "${REF}.fai" ]]; then
    samtools faidx "$REF"
fi

bedtools sort \
    -faidx "${REF}.fai" \
    -i "$CDS_BY_GENE" \
> "$OUT/stickleback_20426_canonical_cds_by_gene.sorted.bed"

mv \
  "$OUT/stickleback_20426_canonical_cds_by_gene.sorted.bed" \
  "$CDS_BY_GENE"

echo "[00] Creating merged CDS target intervals"

cut -f1-3 "$CDS_BY_GENE" |
bedtools merge \
> "$CDS_MERGED"

echo "[00] Calculating CDS lengths"

awk -F'\t' '
BEGIN {OFS="\t"}
{
    length_by_gene[$4]+=$3-$2
    transcript[$4]=$7
}
END {
    print "gene_id","transcript_id","cds_length"

    for (gene in length_by_gene)
        print gene,transcript[gene],length_by_gene[gene]
}
' "$CDS_BY_GENE" |
sort -t$'\t' -k1,1 \
> "$OUT/stickleback_20426_canonical_cds_lengths.tsv"

echo "[00] QA"

printf "Classification genes\t"
wc -l < "$GENES"

printf "Classification transcripts\t"
wc -l < "$TRANSCRIPTS"

printf "GFF3 gene features\t"
awk -F'\t' '$0!~/^#/ && $3=="gene" {n++} END {print n+0}' \
  "$CANON_GFF3"

printf "GFF3 mRNA features\t"
awk -F'\t' '
$0!~/^#/ && ($3=="mRNA" || $3=="transcript") {n++}
END {print n+0}
' "$CANON_GFF3"

printf "Genes with CDS\t"
awk -F'\t' '{print $4}' "$CDS_BY_GENE" |
sort -u |
wc -l

printf "Transcripts with CDS\t"
awk -F'\t' '{print $7}' "$CDS_BY_GENE" |
sort -u |
wc -l

printf "CDS BED rows\t"
wc -l < "$CDS_BY_GENE"

printf "Merged CDS intervals\t"
wc -l < "$CDS_MERGED"

printf "CDS lengths not divisible by 3\t"
awk -F'\t' '
NR>1 && $3%3!=0 {n++}
END {print n+0}
' "$OUT/stickleback_20426_canonical_cds_lengths.tsv"

echo "[00] Classification counts"

awk -F'\t' '
NR>1 {count[$4]++}
END {
    for (class in count)
        print class,count[class]
}
' "$CLASSIFICATION" |
sort

echo "[00] Outputs"
echo "Canonical GFF3: $CANON_GFF3"
echo "Canonical GTF:  $CANON_GTF"
echo "CDS by gene:    $CDS_BY_GENE"
echo "Merged CDS BED: $CDS_MERGED"

#!/usr/bin/env bash
set -u

GENE=${1:?Gene ID required}

ROOT="/path/to/data/genomewide_codeml_kuster"
PHY="$ROOT/06_codeml/phylip/${GENE}.phy"

TREE="/path/to/workspace/Kuster2026_genomewide_reanalysis/genomewide_erc/00_inputs/erc_master_chr21_27_unrooted_paml.tree"

OUT="$ROOT/07_codeml_genomewide"
DIR="$OUT/results/$GENE"

COMPLETE="$DIR/.complete"
RESULT="$DIR/result.tsv"
LOG="$DIR/codeml.stdout.log"
MLC="$DIR/mlc"

mkdir -p "$DIR"

if [[ -s "$COMPLETE" &&
      -s "$RESULT" &&
      -s "$MLC" ]] &&
   grep -q 'lnL(ntime' "$MLC" &&
   grep -q 'omega (dN/dS)' "$MLC"; then

    echo "[skip] $GENE"
    exit 0
fi

if [[ ! -s "$PHY" ]]; then
    echo "[fail] $GENE missing PHYLIP" >&2
    exit 2
fi

if [[ ! -s "$TREE" ]]; then
    echo "[fail] $GENE missing tree" >&2
    exit 2
fi

rm -f \
    "$COMPLETE" \
    "$RESULT" \
    "$DIR/result.tsv.tmp" \
    "$DIR/alignment.phy" \
    "$DIR/tree.nwk" \
    "$DIR/codeml.ctl" \
    "$DIR/mlc" \
    "$DIR/rst" \
    "$DIR/rst1" \
    "$DIR/2NG.dN" \
    "$DIR/2NG.dS" \
    "$DIR/2NG.t" \
    "$DIR/lnf" \
    "$DIR/rub"

cp "$PHY" "$DIR/alignment.phy"
cp "$TREE" "$DIR/tree.nwk"

cat > "$DIR/codeml.ctl" <<'CTL'
seqfile = alignment.phy
treefile = tree.nwk
outfile = mlc

noisy = 0
verbose = 0
runmode = 0

seqtype = 1
CodonFreq = 2
clock = 0
aaDist = 0
model = 0
NSsites = 0
icode = 0

fix_kappa = 0
kappa = 2

fix_omega = 0
omega = 0.2

fix_alpha = 1
alpha = 0
Malpha = 0
ncatG = 8

getSE = 0
RateAncestor = 0
Small_Diff = 5e-7
cleandata = 1
method = 0
CTL

START_EPOCH=$(date +%s)

(
    cd "$DIR"
    codeml codeml.ctl
) > "$LOG" 2>&1

STATUS=$?
END_EPOCH=$(date +%s)
ELAPSED=$((END_EPOCH - START_EPOCH))

SUCCESS="no"

if [[ "$STATUS" -eq 0 &&
      -s "$MLC" ]] &&
   grep -q 'lnL(ntime' "$MLC" &&
   grep -q 'omega (dN/dS)' "$MLC"; then
    SUCCESS="yes"
fi

if [[ "$SUCCESS" == "yes" ]]; then
    LNL=$(
        awk '
        /lnL\(ntime/ {
            for (i=1; i<=NF; i++) {
                if ($i ~ /^-[0-9]/ || $i ~ /^[0-9]/) {
                    value=$i
                }
            }
        }
        END {print value}
        ' "$MLC"
    )

    KAPPA=$(
        awk -F'= *' '
        /kappa \(ts\/tv\)/ {
            print $2
            exit
        }
        ' "$MLC"
    )

    OMEGA=$(
        awk -F'= *' '
        /omega \(dN\/dS\)/ {
            print $2
            exit
        }
        ' "$MLC"
    )

    TREE_DN=$(
        awk -F': *' '
        /tree length for dN:/ {
            print $2
            exit
        }
        ' "$MLC"
    )

    TREE_DS=$(
        awk -F': *' '
        /tree length for dS:/ {
            print $2
            exit
        }
        ' "$MLC"
    )

    {
        printf "gene_id\texit_status\telapsed_seconds\tlnL\tkappa\tomega\ttree_length_dN\ttree_length_dS\n"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$GENE" \
            "$STATUS" \
            "$ELAPSED" \
            "${LNL:-NA}" \
            "${KAPPA:-NA}" \
            "${OMEGA:-NA}" \
            "${TREE_DN:-NA}" \
            "${TREE_DS:-NA}"
    } > "$DIR/result.tsv.tmp"

    mv "$DIR/result.tsv.tmp" "$RESULT"

    {
        echo "gene=$GENE"
        echo "completed=$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
        echo "elapsed_seconds=$ELAPSED"
    } > "$COMPLETE"

    # Inputs are available centrally and can be recreated.
    rm -f \
        "$DIR/alignment.phy" \
        "$DIR/tree.nwk" \
        "$DIR/codeml.ctl" \
        "$DIR/rst" \
        "$DIR/rst1" \
        "$DIR/2NG.dN" \
        "$DIR/2NG.dS" \
        "$DIR/2NG.t" \
        "$DIR/lnf" \
        "$DIR/rub"

    echo "[done] $GENE seconds=$ELAPSED omega=${OMEGA:-NA}"
    exit 0
fi

{
    printf "gene_id\texit_status\telapsed_seconds\tstatus\n"
    printf "%s\t%s\t%s\tfailed\n" \
        "$GENE" "$STATUS" "$ELAPSED"
} > "$DIR/failed.tsv"

echo "[fail] $GENE exit=$STATUS seconds=$ELAPSED" >&2
exit 1

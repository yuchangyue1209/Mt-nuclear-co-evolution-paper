#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Convert an aligned FASTA file to NEXUS format for use in PopART,
including a BEGIN TRAITS block for tagging populations by Environment and Region.
"""

from Bio import SeqIO

def fasta_to_nexus_with_traits(input_fasta, output_nexus):
    """
    1) Read aligned FASTA and convert it into a NEXUS file (BEGIN DATA/CHARACTERS block).
    2) Add a BEGIN TRAITS block with metadata annotations (Environment, Region).
    """

    # Define population metadata
    # Keys = sequence IDs from FASTA; Values = (Environment, Region)
    # Environment: marine, freshwater, recent_colonized
    # Region: Alaska, BC
    pop_map = {
        # ----------- Alaska -----------
        "FG":  ("freshwater",       "Alaska"),
        "LG":  ("freshwater",       "Alaska"),
        "SR":  ("freshwater",       "Alaska"),  # Assuming SR is the marine population
        "SL":  ("freshwater",       "Alaska"),
        "TL":  ("freshwater",       "Alaska"),
        "WB":  ("freshwater",       "Alaska"),
        "WT":  ("freshwater",       "Alaska"),
        "WK":  ("freshwater",       "Alaska"),
        "RS":  ("marine",           "Alaska"),
        "SC":  ("recent_colonized", "Alaska"),
        "LB":  ("freshwater",       "Alaska"),
        "CH":  ("recent_colonized", "Alaska"), 

        # ----------- British Columbia -----------
        "SWA":   ("freshwater",       "BC"),
        "THE":   ("freshwater",       "BC"),
        "JOE":   ("freshwater",       "BC"),
        "BEA":   ("freshwater",       "BC"),
        "MUC":   ("freshwater",       "BC"),
        "PYE":   ("freshwater",       "BC"),
        "AMO":   ("freshwater",       "BC"),
        "SAY":   ("marine",           "BC"),
        "GOS":   ("freshwater",       "BC"),
        "ROB":   ("freshwater",       "BC"),
        "BOOT":  ("freshwater",       "BC"),
        "ECHO":  ("freshwater",       "BC"),
        "FRED":  ("recent_colonized", "BC"),
        "LAW":   ("freshwater",       "BC"),
        "PACH":  ("recent_colonized", "BC"),
    }

    # Read aligned FASTA file
    sequences = list(SeqIO.parse(input_fasta, "fasta"))
    if not sequences:
        raise ValueError("No sequences found in the input FASTA.")

    # Ensure all sequences are aligned (equal length)
    seq_lengths = {len(seq.seq) for seq in sequences}
    if len(seq_lengths) > 1:
        raise ValueError("Sequences differ in length. Alignment is required.")

    n_tax = len(sequences)
    n_char = len(sequences[0].seq)

    with open(output_nexus, "w") as nexus:
        # NEXUS header
        nexus.write("#NEXUS\n\n")

        # TAXA block
        nexus.write("BEGIN TAXA;\n")
        nexus.write(f"  DIMENSIONS NTAX={n_tax};\n")
        nexus.write("  TAXLABELS\n")
        for seq in sequences:
            nexus.write(f"    {seq.id}\n")
        nexus.write("  ;\nEND;\n\n")

        # CHARACTERS block (DNA sequences)
        nexus.write("BEGIN CHARACTERS;\n")
        nexus.write(f"  DIMENSIONS NCHAR={n_char};\n")
        nexus.write("  FORMAT DATATYPE=DNA MISSING=? GAP=-;\n")
        nexus.write("  MATRIX\n")
        for seq in sequences:
            nexus.write(f"    {seq.id} {str(seq.seq)}\n")
        nexus.write("  ;\nEND;\n\n")

        # TRAITS block (environment and region)
        nexus.write("BEGIN TRAITS;\n")
        nexus.write("  DIMENSIONS NTRAITS=2;\n")
        nexus.write("  FORMAT labels=yes missing=? separator=Space;\n")
        nexus.write("  TraitLabels Environment Region;\n")
        nexus.write("  MATRIX\n")
        for seq in sequences:
            seq_id = seq.id
            if seq_id not in pop_map:
                env, reg = ("freshwater", "BC")  # Default fallback
            else:
                env, reg = pop_map[seq_id]
            nexus.write(f"    {seq_id} {env} {reg}\n")
        nexus.write("  ;\nEND;\n")

    print(f"✅ NEXUS file with TRAITS block created: {output_nexus}")


# ---------- Example usage ----------
input_fasta = "/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/aligned_mt_overlap.fasta"
output_nexus = "/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/aligned_mt_overlap_trait.nex"

fasta_to_nexus_with_traits(input_fasta, output_nexus)

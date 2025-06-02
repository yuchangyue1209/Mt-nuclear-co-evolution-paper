# dn/ds calculation
#Merge the files into one fasta
cat /work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/*.fasta > /work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/all_mt_overlap.fasta

13. Realignment
mafft --auto all_mt_overlap.fasta > aligned_mt_overlap.fasta
muscle -in all_mt_overlap.fasta -out mc_aligned_mt_overlap.fasta

14. IQtree 
iqtree -s aligned_mt_overlap.fasta -m GTR+G -bb 1000 -alrt 1000 -nt AUTO
iqtree -s mc_aligned_mt_overlap.fasta -m GTR+G -bb 1000 -alrt 1000 -nt AUTO


/home/cyu/snpEff/data/Gasterosteus_aculeatus_MT/genes.gff
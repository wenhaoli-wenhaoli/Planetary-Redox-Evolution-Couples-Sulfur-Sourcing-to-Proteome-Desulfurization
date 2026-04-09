#===============================================================================
# Build phylogenetic tree of life using GToTree and IQ-TREE
#===============================================================================

# 1. Run GToTree to extract 16 ribosomal marker genes, concatenate and trim
# Input: genome_list.txt (698 representative genomes, including 7 plastid and 12 eukaryotic genomes)
# Output: concatenated.aa.trim.aln (2,176 amino acid positions)

time GToTree -A genome_list.txt -H Universal_Hug_et_al.hmm -G 0.5 -o gtotree_output  -j 100

# 2. Run IQ-TREE with partitioned analysis 

time iqtree -s gtotree_output/Aligned_SCGs.faa -Q gtotree_output/Partitions.txt -pre iqtree_out -m MFP+MERGE -B 1000 -alrt 1000 --nmax 2000 --boot-trees --runs 3  -T AUTO

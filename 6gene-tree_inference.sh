#1===============================================================================
# Build phylogenetic trees for 226 key COG gene families (Single batch)
#1===============================================================================

# For each COG family in cog_fa directory

mkdir -p cog_trees

for f in cog_fa/*.fa; do
    base=$(basename "$f" .fa)
    echo "===> Processing COG $base at $(date)"
    
    # Step 1: Multiple sequence alignment with MAFFT
    mafft --maxiterate 1000 --localpair "$f" > "${base}.alig" || continue
    
    # Step 2: Create COG-specific directory
    mkdir -p "cog_trees/${base}"
    
    # Step 3: Trim poorly aligned regions with TrimAl
    trimal -in "${base}.alig" -out "cog_trees/${base}/${base}.trim" -automated1 || continue
    rm "${base}.alig"
    
    # Step 4: Model selection and tree inference with IQ-TREE
    time iqtree  -s "cog_trees/${base}/${base}.trim" -m MFP -mrate E,I,G,I+G,R \
        -madd C10,C20,C30,C40,C50,C60,EX2,EX3,EHO,UL2,UL3,EX_EHO,LG4M,LG4X,CF4,LG+C10,LG+C20,LG+C30,LG+C40,LG+C50,LG+C60 \
        --nmax 2000 -B 1000 --boot-trees -T AUTO
    echo "===> Completed COG $base at $(date)"
done

echo "All COG tree building completed. Results in cog_trees/"

#2===============================================================================
# Build phylogenetic trees for sulfonate-sulfite interconversion enzymes
#2===============================================================================
# For each enzyme family identified by custom HMM profiles from the species tree

mkdir -p sulfonate_trees

for i in *.fa; do
    base=$(basename "$i" .fa)
    echo "===> Processing sulfonate enzyme $base at $(date)"
    
    # Step 1: Multiple sequence alignment with MAFFT
    mafft --maxiterate 1000 --localpair "$i" > "${base}.alig"
    
    # Step 2: Create enzyme-specific directory
    mkdir -p "sulfonate_trees/${base}"
    
    # Step 3: Trim poorly aligned regions with TrimAl
    trimal -in "${base}.alig" -out "sulfonate_trees/${base}/${base}.trim" -automated1
    
    # Step 4: Model selection and tree inference with IQ-TREE
    time iqtree \
        -s "sulfonate_trees/${base}/${base}.trim" -m MFP -mrate E,I,G,I+G,R \
        -madd C10,C20,C30,C40,C50,C60,EX2,EX3,EHO,UL2,UL3,EX_EHO,LG4M,LG4X,CF4,LG+C10,LG+C20,LG+C30,LG+C40,LG+C50,LG+C60 \
        --nmax 2000 -B 1000 --boot-trees -T AUTO
    
    echo "===> Completed sulfonate enzyme $base at $(date)"
done

echo "All sulfonate enzyme tree building completed. Results in sulfonate_trees/"
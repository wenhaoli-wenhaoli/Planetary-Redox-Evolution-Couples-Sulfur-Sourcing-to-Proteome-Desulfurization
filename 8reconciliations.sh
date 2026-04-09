#!/bin/bash
#===============================================================================
# Reconciliations: 226 COG gene families + sulfonate metabolism enzymes
# with time-calibrated species tree using AleRax
#===============================================================================

# =============================================================================
# PART 1: 226 COG gene families reconciliation
# =============================================================================

# 1. Prepare family files for AleRax
mkdir -p families

tail -n +2 selected_features.tsv | tr -d '\r' | while read -r cog; do
    cog=$(echo "$cog" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    [ -z "$cog" ] && continue
    {
        echo "[FAMILIES]"
        echo "- $cog"
        echo "gene_tree = data/cog_tree/$cog/${cog}.trim.ufboot"
        echo "mapping = data/cog_mappings/$cog.mapping"
    } > "families/${cog}.family"
done

echo "Generated $(ls families/*.family 2>/dev/null | wc -l) COG family files"

# 2. Create output directories for three molecular clock models
mkdir -p output_reconciliation/cog/cir output_reconciliation/cog/ln output_reconciliation/cog/ugam

# 3. Run AleRax for all COG families with three molecular clock models
# --gene-tree-samples 1000: use 1000 bootstrap trees

# CIR model
nohup bash -c '
for f in families/*.family; do
    [ -f "$f" ] && cog=$(basename "$f" .family)
    mkdir -p "output_reconciliation/cog/cir/$cog"
    echo "Processing COG $cog with CIR model at $(date)"
    time alerax -f "$f" \
        -s data/timetree/cir1_sample.chronogram \
        -p "output_reconciliation/cog/cir/$cog/" \
        --gene-tree-samples 1000
done
' > alerax_cog_cir.log 2>&1 &

# LN model
nohup bash -c '
for f in families/*.family; do
    [ -f "$f" ] && cog=$(basename "$f" .family)
    mkdir -p "output_reconciliation/cog/ln/$cog"
    echo "Processing COG $cog with LN model at $(date)"
    time alerax -f "$f" \
        -s data/timetree/ln1_sample.chronogram \
        -p "output_reconciliation/cog/ln/$cog/" \
        --gene-tree-samples 1000
done
' > alerax_cog_ln.log 2>&1 &

# UGAM model
nohup bash -c '
for f in families/*.family; do
    [ -f "$f" ] && cog=$(basename "$f" .family)
    mkdir -p "output_reconciliation/cog/ugam/$cog"
    echo "Processing COG $cog with UGAM model at $(date)"
    time alerax -f "$f" \
        -s data/timetree/ugam1_sample.chronogram \
        -p "output_reconciliation/cog/ugam/$cog/" \
        --gene-tree-samples 1000
done
' > alerax_cog_ugam.log 2>&1 &

# =============================================================================
# PART 2: Sulfonate metabolism enzymes reconciliation (11 enzymes)
# =============================================================================

# 1. Prepare family files for sulfonate enzymes
mkdir -p families_sulfonate

# List of 11 sulfonate-sulfite interconversion enzymes
SULFONATE_ENZYMES="sqdB comA cs smoC sqoD xsc cuyA ssuD suyB tauD hpsG"

for enzyme in $SULFONATE_ENZYMES; do
    {
        echo "[FAMILIES]"
        echo "- $enzyme"
        echo "gene_tree = data/sulfonate_tree/$enzyme/${enzyme}.trim.ufboot"
        echo "mapping = data/sulfonate_mappings/$enzyme.mapping"
    } > "families_sulfonate/${enzyme}.family"
done

echo "Generated $(ls families_sulfonate/*.family 2>/dev/null | wc -l) sulfonate enzyme family files"

# 2. Create output directories for sulfonate enzymes
mkdir -p output_reconciliation/sulfonate/cir output_reconciliation/sulfonate/ln output_reconciliation/sulfonate/ugam

# 3. Run AleRax for sulfonate enzymes with CIR model
nohup bash -c '
for f in families_sulfonate/*.family; do
    [ -f "$f" ] && enzyme=$(basename "$f" .family)
    mkdir -p "output_reconciliation/sulfonate/cir/$enzyme"
    echo "Processing sulfonate enzyme $enzyme with CIR model at $(date)"
    time alerax -f "$f" \
        -s data/timetree/cir1_sample.chronogram \
        -p "output_reconciliation/sulfonate/cir/$enzyme/" \
        --gene-tree-samples 1000
done
' > alerax_sulfonate_cir.log 2>&1 &

# 4. Run AleRax for sulfonate enzymes with LN model
nohup bash -c '
for f in families_sulfonate/*.family; do
    [ -f "$f" ] && enzyme=$(basename "$f" .family)
    mkdir -p "output_reconciliation/sulfonate/ln/$enzyme"
    echo "Processing sulfonate enzyme $enzyme with LN model at $(date)"
    time alerax -f "$f" \
        -s data/timetree/ln1_sample.chronogram \
        -p "output_reconciliation/sulfonate/ln/$enzyme/" \
        --gene-tree-samples 1000
done
' > alerax_sulfonate_ln.log 2>&1 &

# 5. Run AleRax for sulfonate enzymes with UGAM model
nohup bash -c '
for f in families_sulfonate/*.family; do
    [ -f "$f" ] && enzyme=$(basename "$f" .family)
    mkdir -p "output_reconciliation/sulfonate/ugam/$enzyme"
    echo "Processing sulfonate enzyme $enzyme with UGAM model at $(date)"
    time alerax -f "$f" \
        -s data/timetree/ugam1_sample.chronogram \
        -p "output_reconciliation/sulfonate/ugam/$enzyme/" \
        --gene-tree-samples 1000
done
' > alerax_sulfonate_ugam.log 2>&1 &

# =============================================================================
# Clean up unnecessary files to save disk space
# =============================================================================

# Wait for all jobs to complete (optional)
# wait

# Remove large intermediate files
rm -rf output_reconciliation/cog/*/*/reconciliations/origins
rm -rf output_reconciliation/sulfonate/*/*/reconciliations/origins

find output_reconciliation/cog/*/*/reconciliations/all -type f -name "*_perspecies_eventcount_*" -print0 | xargs -0 rm -f
find output_reconciliation/cog/*/*/reconciliations/all -type f -name "*_transfers_*" -print0 | xargs -0 rm -f
find output_reconciliation/sulfonate/*/*/reconciliations/all -type f -name "*_perspecies_eventcount_*" -print0 | xargs -0 rm -f
find output_reconciliation/sulfonate/*/*/reconciliations/all -type f -name "*_transfers_*" -print0 | xargs -0 rm -f

echo "================================================================================"
echo "All reconciliations completed!"
echo "COG families: output_reconciliation/cog/"
echo "Sulfonate enzymes: output_reconciliation/sulfonate/"
echo "================================================================================"
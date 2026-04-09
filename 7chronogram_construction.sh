#===============================================================================
# Molecular clock analysis with PhyloBayes v4.1
#===============================================================================
# Calibration points:
# 
# 1. Root (LUCA divergence: Bacteria vs Archaea)
#    - Uniform distribution: 4.4 - 3.5 Ga
#    - Based on earliest evidence for liquid water and oldest microfossils
#
# 2. MRCA of TACK superphylum and Euryarchaeota (methanogens)
#    - Minimum age: 3.42 Ga
#    - Filamentous microfossils from Barberton greenstone belt, South Africa
#
# 3. MRCA of thylakoid-containing cyanobacteria
#    - Calibration: 1.75 Ga (Navifusa majensis fossil)
#
# 4. MRCA of heterocyst-forming cyanobacteria
#    - Calibration: 1.0 Ga (Polysphaeroides filiformis fossil)
#
# 5. Divergence of eukaryotes
#    - Minimum age: 1.64 Ga (Tappania plana fossil)
#
# 6. Divergence of Chromatiaceae (purple sulfur bacteria)
#    - Minimum age: 1.64 Ga (Okenane biomarker from Barney Creek Formation)
#===============================================================================

# Create output directories for three molecular clock models
mkdir -p ln cir ugam

#===============================================================================
# Model 1: Autocorrelated Lognormal clock (LN)
#===============================================================================
# Two independent chains run for >50,000 cycles

nohup time ../data/pb -d alig.phy -T rm.treefile -r micro.outgroup -cal calib -x 5 80000 -catfix C20 -ln ln/ln1 > ln1.log 2>&1 &
nohup time ../data/pb -d alig.phy -T rm.treefile -r micro.outgroup -cal calib -x 5 80000 -catfix C20 -ln ln/ln2 > ln2.log 2>&1 &


#===============================================================================
# Model 2: Autocorrelated Cox-Ingersoll-Ross clock (CIR)
#===============================================================================
# Two independent chains run for >50,000 cycles

nohup time ../data/pb -d alig.phy -T rm.treefile -r micro.outgroup -cal calib -x 5 80000 -catfix C20 -cir cir/cir1 > cir1.log 2>&1 &
nohup time ../data/pb -d alig.phy -T rm.treefile -r micro.outgroup -cal calib -x 5 80000 -catfix C20 -cir cir/cir2 > cir2.log 2>&1 &


#===============================================================================
# Model 3: Uncorrelated Gamma multipliers clock (UGAM)
#===============================================================================
# Two independent chains run for >50,000 cycles

nohup time ../data/pb -d alig.phy -T rm.treefile -r micro.outgroup -cal calib -x 5 80000 -catfix C20 -ugam ugam/ugam1 > ugam1.log 2>&1 &
nohup time ../data/pb -d alig.phy -T rm.treefile -r micro.outgroup -cal calib -x 5 80000 -catfix C20 -ugam ugam/ugam2 > ugam2.log 2>&1 &


#===============================================================================
# Check convergence using tracecomp
#===============================================================================
# Convergence criteria: max discrepancy < 0.3 and min effective size > 50

# Check LN model convergence
../data/tracecomp -x 10000 ln/ln1 ln/ln2 > ln_convergence.txt

# Check CIR model convergence
../data/tracecomp -x 10000 cir/cir1 cir/cir2 > cir_convergence.txt

# Check UGAM model convergence
../data/tracecomp -x 10000 ugam/ugam1 ugam/ugam2 > ugam_convergence.txt

#===============================================================================
# Generate chronograms using readdiv (discard first 25% as burn-in)
#===============================================================================

# LN model chronogram
../data/readdiv -x 10000 ln/ln1 > ln_chronogram.tre

# CIR model chronogram
../data/readdiv -x 10000 cir/cir1 > cir_chronogram.tre

# UGAM model chronogram
../data/readdiv -x 10000 ugam/ugam1 > ugam_chronogram.tre

echo "All molecular clock analyses completed!"
echo "Chronograms saved as: ln_chronogram.tre, cir_chronogram.tre, ugam_chronogram.tre"
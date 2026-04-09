#!/bin/bash
# run eggNOG-mapper for order-level representative genomes

INPUT_DIR="data/processed/order/"
OUTPUT_DIR="data/processed/eggnog_output/"
EGGNOG_DB="/path/to/eggNOG_db"
THREADS=64

# Create output directory
mkdir -p $OUTPUT_DIR
cat ${INPUT_DIR}*.faa > ${OUTPUT_DIR}/order.fa
time emapper.py -i ${OUTPUT_DIR}/order.fa -o ${OUTPUT_DIR}/all_annotations --cpu $THREADS -m diamond --data_dir $EGGNOG_DB --query_cover 50 --evalue 1e-7
echo "eggNOG-mapper completed. Results saved to ${OUTPUT_DIR}/all_annotations.emapper.annotations"
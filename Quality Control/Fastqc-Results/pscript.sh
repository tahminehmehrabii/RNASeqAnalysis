# CONTENTS OF pscript.sh
#!/bin/bash
DATA_DIR=~/RNASeq_PROJECT/DATA
OUTPUT_DIR=~/RNASeq_PROJECT/fastqc-results

for file in $DATA_DIR/*.fastq.gz
do
    echo "Running FastQC on $file"
    fastqc $file -o $OUTPUT_DIR
done

echo "FastQC analysis complete. Results are in the $OUTPUT_DIR directory."
# BASENAME=$(sed "${SGE_TASK_ID}q;d" ./data/annotations/basenames.txt)

# handle arguments
ARGS=$(awk -F "," -v task_id="${SGE_TASK_ID}" 'NR==task_id+1 { print $2, $4 }' ./data/annotations/sample_annotations_beads.csv)
BASENAME=$(cut -d " " -f 1 <<< $ARGS)
DIRECTORY=$(cut -d " " -f 2 <<< $ARGS)
mkdir -p ./data/genomecov/"${DIRECTORY}"

# index BAM
# samtools index ../../"${DIRECTORY}"/data/alignments/"${BASENAME}".bam

# genome coverage
bedtools genomecov -ibam ../../"${DIRECTORY}"/data/alignments/"${BASENAME}".bam -dz -split > ./data/genomecov/"${DIRECTORY}"/"${BASENAME}".full.txt

# only keep rRNAs
awk '$1 == "I" && $2 > 15044838 && $2 < 15088346 {print $0}' ./data/genomecov/"${DIRECTORY}"/"${BASENAME}".full.txt
rm -f ./data/genomecov/"${DIRECTORY}"/"${BASENAME}".full.txt

# compress
gzip -f ./data/genomecov/"${DIRECTORY}"/"${BASENAME}".txt

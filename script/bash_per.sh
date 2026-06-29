#sample="aKG"
#round="19"

sample=$1
round=$2
SAMPLE="${sample}${round}"
echo ">> Start analysis for ${SAMPLE}"

# Resolve repository root based on this script location.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"



# 1. Create directories for results and logs.
echo "1. Create directories for results and logs"
mkdir -p ${DIR}/results/${SAMPLE}/QC
mkdir -p ${DIR}/log
echo "1. Done ------------------------------------------------------------"
echo ""


# 2. Define input fastq files.
echo "2. Define input fastq files"
fq_in_R1=`ls ${DIR}/data/${sample}/${SAMPLE}_L001/${SAMPLE}*_R1_001.fastq.gz`
fq_in_R2=`ls ${DIR}/data/${sample}/${SAMPLE}_L001/${SAMPLE}*_R2_001.fastq.gz`
echo "-Input fastq files:"
echo "--R1: ${fq_in_R1}" 
echo "--R2: ${fq_in_R2}"
echo "2. Done ------------------------------------------------------------"
echo ""


# 3. Quality check of raw fastq files.
echo "3. Quality check of raw fastq files (FastQC)"
fastqc -q ${fq_in_R1} ${fq_in_R2} \
	-o ${DIR}/results/${SAMPLE}/QC \
	> ${DIR}/log/${SAMPLE}_RawFastQC.txt
echo "3. Done ------------------------------------------------------------"
echo ""


# 4. Adaptor trimming.
echo "4. Adaptor trimming"
cutadapt \
	-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
	-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
	-q 30,30 \
	-m 30 \
	-o ${DIR}/results/${SAMPLE}/${SAMPLE}_R1_trimmed.fastq.gz \
	-p ${DIR}/results/${SAMPLE}/${SAMPLE}_R2_trimmed.fastq.gz \
	${fq_in_R1} \
	${fq_in_R2} \
	> ${DIR}/log/${SAMPLE}_AdapterTrimming.txt
echo "4. Done ------------------------------------------------------------"
echo ""


# 5. Merge R1,R2 reads to make DNA Fragment.
echo "5. Merge R1,R2 reads to make DNA Fragment"
pear \
	-f ${DIR}/results/${SAMPLE}/${SAMPLE}_R1_trimmed.fastq.gz \
	-r ${DIR}/results/${SAMPLE}/${SAMPLE}_R2_trimmed.fastq.gz \
	-o ${DIR}/results/${SAMPLE}/${SAMPLE}_PEAR \
	-v 10 \
	> ${DIR}/log/${SAMPLE}_PEAR.txt
echo "5. Done ------------------------------------------------------------"
echo ""


# 6. Quality check of merged fastq files.
echo "6. Quality check of merged fastq files (FastQC)"
# make fastqc silent.
fastqc -q ${DIR}/results/${SAMPLE}/${SAMPLE}_PEAR.assembled.fastq \
	-o ${DIR}/results/${SAMPLE}/QC \
	> ${DIR}/log/${SAMPLE}_MergedFastQC.txt
echo "6. Done ------------------------------------------------------------"
echo ""


# Define DNAzyme sequence.
# 'Enzyme' sequence: 5'- TAATACGACTCACTATAGGG - N40 - ATCTGACGGTAACGCTATAGTGTCACCTAAATAGC -3’
#                        ***TACGACTCACTATAGGG - N40 - ATCTGACGGTAACGCTATAGTGTCACCTAAAT***
## Remove first (5`) and last (3`) 3 bases for buffering. Since cutadapt can remove a few more/less bases, we will allow 3 bases buffering.
## Define flanking sequence.

seq_5="TACGACTCACTATAGGG"
seq_3="ATCTGACGGTAACGCTATAGTGTCACCTAAAT"
seq_3_rc="ATTTAGGTGACACTATAGCGTTACCGTCAGAT" # Reverse complemented

# 7. Select DNA fragments with expected flanking sequences.
echo "7. Select DNA fragments with expected flanking sequences"

# 7.1 Select 5` flanking sequence.
echo "7.1 Select 5' flanking sequence"
grep --no-group-separator -B1 -A2 ${seq_5} ${DIR}/results/${SAMPLE}/${SAMPLE}_PEAR.assembled.fastq \
	> ${DIR}/results/${SAMPLE}/${SAMPLE}_PEAR.assembled.5to3.fastq
	
# 7.2 Select 3` flanking sequence.
echo "7.2 Select 3' flanking sequence"
grep --no-group-separator -B1 -A2 ${seq_3_rc} ${DIR}/results/${SAMPLE}/${SAMPLE}_PEAR.assembled.fastq \
        > ${DIR}/results/${SAMPLE}/${SAMPLE}_PEAR.assembled.3to5.fastq

# 7.3 Reversecomplement 3to5 sequences.
echo "7.3 Reversecomplement 3to5 sequences"
fastx_reverse_complement \
	-i ${DIR}/results/${SAMPLE}/${SAMPLE}_PEAR.assembled.3to5.fastq \
	-o ${DIR}/results/${SAMPLE}/${SAMPLE}_PEAR.assembled.3to5.RC.fastq

# 7.4 Concat sequences in the same 5-3 direction.
echo "7.4 Concat sequences in the same 5-3 direction"
cat \
	${DIR}/results/${SAMPLE}/${SAMPLE}_PEAR.assembled.5to3.fastq \
	${DIR}/results/${SAMPLE}/${SAMPLE}_PEAR.assembled.3to5.RC.fastq \
	> ${DIR}/results/${SAMPLE}/${SAMPLE}_PEAR.assembled.5to3.merged.fastq
echo "7. Done ------------------------------------------------------------"
echo ""


# 8. Select DNA fragments with expected 5`/3` flanking sequences.
echo "8. Select DNA fragments with expected 5'/3' flanking sequences"
grep "TACGACTCACTATAGGG" ${DIR}/results/${SAMPLE}/${SAMPLE}_PEAR.assembled.5to3.merged.fastq |\
	grep "ATCTGACGGTAACGCTATAGTGTCACCTAAAT" \
	> ${DIR}/results/${SAMPLE}/${SAMPLE}_PEAR.assembled.5to3.merged.final.txt
echo "8. Done ------------------------------------------------------------"
echo ""


# 9. Count N40 sequence motif/position frequency and make 2D heatmap.
echo "9. Count N40 sequence motif/position frequency and make 2D heatmap"
echo "9.1 Use all N40 sequences (w/o deduplication)."
python ${SCRIPT_DIR}/MotifFreq_byRounds_v1.py \
	${DIR}/results \
	${sample} \
	${round} \
	${seq_5} \
	${seq_3} \
	"UseAll"
	
echo "9.2 Use unique N40 sequences (w/ deduplication)."
python ${SCRIPT_DIR}/MotifFreq_byRounds_v1.py \
	${DIR}/results \
	${sample} \
	${round} \
	${seq_5} \
	${seq_3} \
	"Deduplicate"
echo "9. Done ------------------------------------------------------------"
echo ""
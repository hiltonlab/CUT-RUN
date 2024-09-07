#!/bin/bash

# Function to display a message and exit the script
function exit_with_message {
    echo "$1"
    exit 1
}

# Check if the correct number of arguments are provided
if [ "$#" -ne 8 ]; then
    echo "Usage: $0 <input_R1> <input_R2> <ref_genome> <ecoli_genome> <output_prefix> <bigwig_output> <output_dir> <fastqc_report_dir>"
    exit 1
fi

# Assign command-line arguments to variables
INPUT_R1="$1"
INPUT_R2="$2"
REF_GENOME="$3"
ECOILI_GENOME="$4"
OUTPUT_PREFIX="$5"
BIGWIG_OUTPUT="$6"
OUTPUT_DIR="$7"
FASTQC_REPORT_DIR="$8"

# Check if input files exist
[ ! -f "$INPUT_R1" ] && exit_with_message "Input R1 file does not exist: $INPUT_R1"
[ ! -f "$INPUT_R2" ] && exit_with_message "Input R2 file does not exist: $INPUT_R2"
[ ! -f "$REF_GENOME" ] && exit_with_message "Reference genome file does not exist: $REF_GENOME"
[ ! -f "$ECOILI_GENOME" ] && exit_with_message "E. coli genome file does not exist: $ECOILI_GENOME"

# Create the output directories if they do not exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$FASTQC_REPORT_DIR"

# Define file paths
TRIMMED_R1_OUT="$OUTPUT_DIR/$(basename ${INPUT_R1%_001.fastq.gz}_trimmed.fastq.gz)"
TRIMMED_R2_OUT="$OUTPUT_DIR/$(basename ${INPUT_R2%_001.fastq.gz}_trimmed.fastq.gz)"
SAM_HG38="$OUTPUT_DIR/$(basename ${OUTPUT_PREFIX}_hg38.sam)"
SAM_ECOILI="$OUTPUT_DIR/$(basename ${OUTPUT_PREFIX}_ecoli.sam)"
BAM_HG38="$OUTPUT_DIR/$(basename ${SAM_HG38%.sam}.bam)"
BAM_ECOILI="$OUTPUT_DIR/$(basename ${SAM_ECOILI%.sam}.bam)"
SORTED_BAM_HG38="$OUTPUT_DIR/$(basename ${BAM_HG38%.bam}_sorted.bam)"
SORTED_BAM_ECOILI="$OUTPUT_DIR/$(basename ${BAM_ECOILI%.bam}_sorted.bam)"
STATS_HG38="$OUTPUT_DIR/$(basename ${OUTPUT_PREFIX}_hg38_sorted.txt)"
STATS_ECOILI="$OUTPUT_DIR/$(basename ${OUTPUT_PREFIX}_ecoli_sorted.txt)"
BIGWIG_OUTPUT="$OUTPUT_DIR/$(basename ${BIGWIG_OUTPUT}.bigWig)"

# Trimming reads
echo "Trimming reads..."
fastp -i "$INPUT_R1" -o "$TRIMMED_R1_OUT" -I "$INPUT_R2" -O "$TRIMMED_R2_OUT" -5 20 -3 20 -q 20 -l 25 --cut_mean_quality 15

# Quality control of trimmed reads
echo "Running FastQC on trimmed reads..."
fastqc "$TRIMMED_R1_OUT" "$TRIMMED_R2_OUT" --outdir="$FASTQC_REPORT_DIR"

# Align reads to the reference genomes
echo "Aligning reads to hg38..."
bwa mem "$REF_GENOME" -t 32 "$TRIMMED_R1_OUT" "$TRIMMED_R2_OUT" > "$SAM_HG38"

echo "Aligning reads to E. coli..."
bwa mem "$ECOILI_GENOME" -t 32 "$TRIMMED_R1_OUT" "$TRIMMED_R2_OUT" > "$SAM_ECOILI"

# Convert SAM to BAM, sort, and index BAM files
echo "Converting SAM to BAM and sorting..."
samtools view -bS "$SAM_HG38" > "$BAM_HG38"
samtools view -bS "$SAM_ECOILI" > "$BAM_ECOILI"
samtools sort "$BAM_HG38" > "$SORTED_BAM_HG38"
samtools sort "$BAM_ECOILI" > "$SORTED_BAM_ECOILI"
samtools index "$SORTED_BAM_HG38"
samtools index "$SORTED_BAM_ECOILI"

# Generate stats for BAM files
echo "Generating BAM statistics..."
bamtools stats -in "$SORTED_BAM_ECOILI" > "$STATS_ECOILI"
bamtools stats -in "$SORTED_BAM_HG38" > "$STATS_HG38"

# Count reads in BAM files
sample_count=$(samtools view -c "$SORTED_BAM_HG38")
ecoli_count=$(samtools view -c "$SORTED_BAM_ECOILI")

# Calculate normalization scale factor
normalization=$(echo "scale=2; ($ecoli_count / $sample_count) * 100" | bc)
scale_factor=$(echo "scale=2; 1 / ($normalization / 100)" | bc)

# Create BigWig file with normalization
echo "Creating BigWig file..."
bamCoverage -b "$SORTED_BAM_HG38" -o "$BIGWIG_OUTPUT" --normalizeUsing RPKM --scaleFactor "$scale_factor" --binSize 20 --smoothLength 200

echo "Script completed successfully."


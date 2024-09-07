#!/bin/bash

# Function to display a message and exit the script
function exit_with_message {
    echo "$1"
    exit 1
}

# Check if the correct number of arguments are provided
if [ "$#" -ne 7 ]; then
    echo "Usage: $0 <sample_prefixes> <ref_genome> <ecoli_genome> <output_dir> <fastqc_report_dir> <bin_size> <smooth_length>"
    exit 1
fi

# Assign command-line arguments to variables
SAMPLE_PREFIXES="$1"  # Space-separated list of sample prefixes
REF_GENOME="$2"
ECOILI_GENOME="$3"
OUTPUT_DIR="$4"
FASTQC_REPORT_DIR="$5"
BIN_SIZE="$6"
SMOOTH_LENGTH="$7"

# Convert the space-separated list into an array
IFS=' ' read -r -a SAMPLE_PREFIX_ARRAY <<< "$SAMPLE_PREFIXES"

# Check if input files exist
[ ! -f "$REF_GENOME" ] && exit_with_message "Reference genome file does not exist: $REF_GENOME"
[ ! -f "$ECOILI_GENOME" ] && exit_with_message "E. coli genome file does not exist: $ECOILI_GENOME"

# Create the output directories if they do not exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$FASTQC_REPORT_DIR"

# Loop over each sample prefix
for SAMPLE_PREFIX in "${SAMPLE_PREFIX_ARRAY[@]}"; do
    # Define sample-specific output directory
    SAMPLE_OUTPUT_DIR="$OUTPUT_DIR/${SAMPLE_PREFIX}"
    SAMPLE_FASTQC_REPORT_DIR="$FASTQC_REPORT_DIR/${SAMPLE_PREFIX}"

    # Create the sample-specific output directories if they do not exist
    mkdir -p "$SAMPLE_OUTPUT_DIR"
    mkdir -p "$SAMPLE_FASTQC_REPORT_DIR"

    # Define file paths based on sample prefix
    INPUT_R1="${SAMPLE_PREFIX}_R1_001.fastq.gz"
    INPUT_R2="${SAMPLE_PREFIX}_R2_001.fastq.gz"
    TRIMMED_R1_OUT="$SAMPLE_OUTPUT_DIR/$(basename ${INPUT_R1%_001.fastq.gz}_trimmed.fastq.gz)"
    TRIMMED_R2_OUT="$SAMPLE_OUTPUT_DIR/$(basename ${INPUT_R2%_001.fastq.gz}_trimmed.fastq.gz)"
    SAM_HG38="$SAMPLE_OUTPUT_DIR/$(basename ${SAMPLE_PREFIX}_hg38.sam)"
    SAM_ECOILI="$SAMPLE_OUTPUT_DIR/$(basename ${SAMPLE_PREFIX}_ecoli.sam)"
    BAM_HG38="$SAMPLE_OUTPUT_DIR/$(basename ${SAM_HG38%.sam}.bam)"
    BAM_ECOILI="$SAMPLE_OUTPUT_DIR/$(basename ${SAM_ECOILI%.sam}.bam)"
    SORTED_BAM_HG38="$SAMPLE_OUTPUT_DIR/$(basename ${BAM_HG38%.bam}_sorted.bam)"
    SORTED_BAM_ECOILI="$SAMPLE_OUTPUT_DIR/$(basename ${BAM_ECOILI%.bam}_sorted.bam)"
    STATS_HG38="$SAMPLE_OUTPUT_DIR/$(basename ${SAMPLE_PREFIX}_hg38_sorted.txt)"
    STATS_ECOILI="$SAMPLE_OUTPUT_DIR/$(basename ${SAMPLE_PREFIX}_ecoli_sorted.txt)"
    BIGWIG_OUTPUT="$SAMPLE_OUTPUT_DIR/$(basename ${SAMPLE_PREFIX}.bigWig)"

    # Debugging: Print file paths
    echo "Processing sample: $SAMPLE_PREFIX"
    echo "Input R1: $INPUT_R1"
    echo "Input R2: $INPUT_R2"
    echo "Trimmed R1 Output: $TRIMMED_R1_OUT"
    echo "Trimmed R2 Output: $TRIMMED_R2_OUT"
    echo "SAM (hg38): $SAM_HG38"
    echo "SAM (E. coli): $SAM_ECOILI"
    echo "BAM (hg38): $BAM_HG38"
    echo "BAM (E. coli): $BAM_ECOILI"
    echo "Sorted BAM (hg38): $SORTED_BAM_HG38"
    echo "Sorted BAM (E. coli): $SORTED_BAM_ECOILI"
    echo "BigWig Output: $BIGWIG_OUTPUT"

    # Check if input files exist
    [ ! -f "$INPUT_R1" ] && exit_with_message "Input R1 file does not exist: $INPUT_R1"
    [ ! -f "$INPUT_R2" ] && exit_with_message "Input R2 file does not exist: $INPUT_R2"

    # Trimming reads
    echo "Trimming reads for sample: $SAMPLE_PREFIX..."
    fastp -i "$INPUT_R1" -o "$TRIMMED_R1_OUT" -I "$INPUT_R2" -O "$TRIMMED_R2_OUT" -5 20 -3 20 -q 20 -l 25 --cut_mean_quality 15

    # Quality control of trimmed reads
    echo "Running FastQC on trimmed reads for sample: $SAMPLE_PREFIX..."
    fastqc -t 64 "$TRIMMED_R1_OUT" "$TRIMMED_R2_OUT" --outdir="$SAMPLE_FASTQC_REPORT_DIR"

    # Align reads to the reference genomes
    echo "Aligning reads to hg38 for sample: $SAMPLE_PREFIX..."
    bwa mem "$REF_GENOME" -t 64 "$TRIMMED_R1_OUT" "$TRIMMED_R2_OUT" > "$SAM_HG38"

    echo "Aligning reads to E. coli for sample: $SAMPLE_PREFIX..."
    bwa mem "$ECOILI_GENOME" -t 64 "$TRIMMED_R1_OUT" "$TRIMMED_R2_OUT" > "$SAM_ECOILI"

    # Convert SAM to BAM, sort, and index BAM files
    echo "Converting SAM to BAM and sorting for sample: $SAMPLE_PREFIX..."
    samtools view -bS "$SAM_HG38" > "$BAM_HG38"
    samtools view -bS "$SAM_ECOILI" > "$BAM_ECOILI"
    samtools sort -@ 64 "$BAM_HG38" > "$SORTED_BAM_HG38"
    samtools sort -@ 64 "$BAM_ECOILI" > "$SORTED_BAM_ECOILI"
    samtools index "$SORTED_BAM_HG38"
    samtools index "$SORTED_BAM_ECOILI"

    # Generate stats for BAM files
    echo "Generating BAM statistics for sample: $SAMPLE_PREFIX..."
    bamtools stats -in "$SORTED_BAM_ECOILI" > "$STATS_ECOILI"
    bamtools stats -in "$SORTED_BAM_HG38" > "$STATS_HG38"

    # Count reads in BAM files
    sample_count=$(samtools view -c "$SORTED_BAM_HG38")
    ecoli_count=$(samtools view -c "$SORTED_BAM_ECOILI")

    # Calculate normalization scale factor
    normalization=$(echo "scale=2; ($ecoli_count / $sample_count) * 100" | bc)
    scale_factor=$(echo "scale=2; 1 / ($normalization / 100)" | bc)

    # Create BigWig file with normalization
    echo "Creating BigWig file for sample: $SAMPLE_PREFIX..."
    bamCoverage -b "$SORTED_BAM_HG38" -o "$BIGWIG_OUTPUT" --numberOfProcessors 64 --normalizeUsing RPKM --scaleFactor "$scale_factor" --binSize "$BIN_SIZE" --smoothLength "$SMOOTH_LENGTH"

    echo "Processing completed for sample: $SAMPLE_PREFIX"
done

echo "Script completed successfully."


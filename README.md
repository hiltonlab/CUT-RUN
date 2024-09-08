# CUT&RUN
CUT&RUN is the leading technology for mapping genomic enrichment of chromatin targets, such as transcription factors and histone PTMs. This script can be used for CUT&amp;RUN analysis starting from raw fastq to normalized bigwigs (ecoli spike-ins) for visualization in IGV. 

# Authors
  * Selvalakshmi Selvaraj Anand (https://github.com/Selvalakshmi27)
  * Guy Bedford
# Prerequisites
Download and index ```hg38``` and ```ecoli``` reference genomes before starting analysis.
  * fastqc (https://github.com/s-andrews/FastQC)
  * fastp (https://github.com/OpenGene/fastp)
  * bwa (https://github.com/lh3/bwa)
  * samtools (https://github.com/samtools/samtools)
  * bamtools (https://github.com/pezmaster31/bamtools)
  * deepTools (https://deeptools.readthedocs.io/en/latest/)

# Usage
For one sample at a time use ```cut_run.sh``` 


```
./cut_run.sh <input_R1> <input_R2> <hg38_ref> <ecoli_ref> <output_prefix> <bigwig_prefix> <output_dir> <fastqc_dir>
```

```input_R1```  - Raw fastq.gz R1

```input_R2```  - Raw fastq.gz R2

```hg38_ref```  - hg38 reference genome directory with the index files

```ecoli_ref``` - ecoli reference genome directory with the index files

```output_prefix``` - Prefix for the output files

```bigwig_prefix``` - Prefix for the output bigwig files

```output_dir``` - Output directory

```fastqc_dir``` - Output directory for fastqc results

Example Usage 

```
./cut_run.sh sample1_R1.fastq.gz sample1_R2.fastq.gz hg38_genome/hg38.fa ecoli_genome/ecoli.fa sample1 sample1 /path/to/output /path/to/fastqc_output
```

For multiple samples at the same time use ```cut_run_for.sh```

```
./cut_run_for.sh <sample_prefix> <hg38_ref> <ecoli_ref> <output_dir> <fastqc_report_dir> <bin_size> <smooth_length>
```

```sample_prefix```  - Prefix of the samples to run

```hg38_ref```  - hg38 reference genome directory with the index files

```ecoli_ref``` - ecoli reference genome directory with the index files

```output_dir``` - Output directory

```fastqc_report_dir``` - Output directory for fastqc results

```bin_size``` - Specify bin size for normalization of the bigwig files

```smooth_length``` - Specify smooth length for the bigwig files 

Example Usage

```
./cut_run_for.sh "sample1 sample2 sample3" hg38_genome/hg38.fa ecoli_genome/ecoli.fa /path/to/results /path/to/fastqc_reports 20 200
```


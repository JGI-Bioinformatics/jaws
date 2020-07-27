# DAP seq pipeline

## Summary
This pipeline is intended for analysis of DAP-Seq libraries.

The main steps are:
- adapter trimming and quality filtering with BBTools
- alignment to ref genome with bowtie2
- generate bigWig with deeptools
- peak-calling with macs2
- motif calling with meme-suite
- assigning peaks to genes

Negative control samples are optional (but recommended)

All samples submitted together in a single inputs.json file are treated as one batch:
- they should all align to the same genome
- the negative control sample(s) (if provided) are applied to all experimental samples
- all negative control samples are merged and treated as a single bam file

#### Required variables are:
```
adapters:           [file path]             # illumina adapters sequences in fasta format
genome_fasta:       [file path]             # reference genome seq in fasta format
bt2index_file1:     [file path]             # one of the bowtie2 reference index files (the rest are inferred)
effgsize:           [int]                   # effective/mappable genome size
genes_gff:          [file path]             # sorted gene annotation file in gff format
bgmodel:            [file path]             # background model for motif calling
outdir:             [file path]             # where to copy the final results
amplified:          [string]                # how many amplification cycles for input libraries?
maxfrags:           [int]                   # subsample fastq to at most this many fragments
find_motifs:        [boolean]               # if set to false, motif calling is skipped
ctl_raw_fastqs:     [array of file paths]   # absolute paths to control fastq files (.gz OK)
expt_raw_fastqs:    [array of file paths]   # absolute paths to experimental fastq files (.gz OK)
library_names_map:  [raw_fastq: name]       # library names for each fastq (for tracking)
sample_names_map:   [raw_fastq: name]       # sample names for each fastq (descriptive names)
```
Note: `ctl_raw_fastqs` _must_ be specified, but it can be an empty array: `[]`

## Tools required (these are loaded in the docker images)

### python3.6+ environment
- samtools
- deeptools
- bbtools (bbduk and reformat.sh)
- bedtools
- bowtie2
- biopython
- matplotlib

### python2.7 environment
- meme-suite
- macs2

## Detailed task info
### trimAlign
**this step is run for both control and experimental samples**
- subsample fastq file
- trim Illumina adapters using bbduk (BBTools)
- align to reference genome using bowtie2
- calculate median fragment size
- generate BigWig file using deeptools

### mergeBams
**this step is run for control samples only**
- merge all control bam files into a single bam using samtools

### findPeaks
**this and all subsequent steps are run for experimental samples only** \
**if zero peaks pass filtering: motifInputs, findMotifs and assignGenes are skipped**

- find peaks using macs2 and with the control bam for background (if provided)
- filter for peaks that are ≥5-fold above background

### motifInputs
- python script to put the peak regions into a fasta file for meme input

*generates two files:*

1. for the entire peak regions
2. only the summits +/- 30bp

### findMotifs
- meme-suite to find motifs
- fimo to map motifs back to sequences

*generates 4 motifs, each in its own subdirectory:*

- entire peak region (regular mode)
- entire peak region (palindromic mode)
- summits +/- 30bp (regular mode)
- summits +/- 30bp (palindromic mode)

### assignGenes
- bedtools to assign peaks to genes

*assignment is done by the following steps:*

1. trim all peaks by 1/3 the avg. fragment size on each end
2. remove any peaks that have been trimmed to <1bp
3. trimmed peaks that are entirely overlapping a gene are written to `<NAME>_overlap_genes.bed`
4. remaining intergenic peaks are assigned to genes that are directly adjacent and oriented away from the peak in `<NAME>_assigned_genes.bed`

### dapStats
- calculate FRIP score using samtools and bedtools
- plot peak strength histogram and peak significance (p-score) vs summit height
- make a table with basic DAP stats

### copyOutput
- copy the output to the specified `<outdir>` in a subdirectory `<NAME>`: \
    `<NAME>` = `<sample_name>`+`_`+`<library_name>`

1. `<NAME>_trim_stats.txt`            bbduk stats
2. `<NAME>_align_stats.txt`           bowtie2 stats
3. `<NAME>_macs2_stats.txt`           macs2 stats
4. `<NAME>_dap_stats.tsv`             summary stats for library
5. `<NAME>_reference.txt`             reference genome files path
6. `<NAME>.bam` and `<NAME>.bai`      aligned reads and index
7. `<NAME>_dedup_norm.bw`:            coverage (BigWig)
8. `<NAME>_peaks.narrowPeak`:         macs2 peaks output
9. `<NAME<_peaks_filt.narrowPeak`:    peaks filtered for  ≥5-fold above background
10. `<NAME>_overlap_genes.bed`:       peaks within genes
11. `<NAME>_assigned_genes.bed`:      peaks upstream of genes
12. `<NAME>_peak_plot.png`:           histogram and scatterplot of peaks
13.  `<NAME>_meme_*` (4 directories): 4 different motifs and their locations in peaks \
    - `<NAME>_meme_peaks`\
    - `<NAME>_meme_peaks_pal`\
    - `<NAME>_meme_summits`\
    - `<NAME>_meme_summits_pal`
14. `<NAME>_trim_stats.txt`:           output from bbduk
15. `<NAME>_align_stats.txt`:          output from bowtie2
16. `<NAME>_macs2_stats.txt`:         output from macs2
17. `<NAME>_dap_stats.tsv`:           number of aligned reads, peaks, and FRIP score

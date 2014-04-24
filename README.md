BreaKmer
========

A method to identify genomic structural variation in target regions/genes from reference-aligned high-throughput sequence data. It uses a “kmer” strategy to assemble misaligned sequence reads for predicting insertions, deletions, inversions, tandem duplications, and translocations at base-pair resolution.

Installation
----------

Prior to installation the following are required for installation:
- [Biopython](http://biopython.org/wiki/Main_Page)
- [Pysam](https://code.google.com/p/pysam/)

Download the python scripts and run the command:
```
python setup.py install
```
Use appropriate commands for installing locally:
```
python setup.py install --user
```
Once installed, BreaKmer can be run with the following command:
```
python breakmer.py <options> <path to config file>
```

Requirements
---------

### Programs
- BLAT standalone and server binaries ([blat, gfServer, gfClient](http://hgdownload.cse.ucsc.edu/admin/exe/)).
  - Re-alignment to reference sequence.
- [Cutadapt](https://code.google.com/p/cutadapt/)
  - Trims adapter sequence from aligned reads.
- [Jellyfish](http://www.cbcb.umd.edu/software/jellyfish/) 
  - Generating kmers.

Paths to these binaries need to be specified in the configuration file input to breakmer.py

### Data
- BreaKmer requires sequence reads that have been aligned to a reference sequence in a binary alignment format (BAM). The alignment program must soft-clip or trim reads that can be partially aligned to the reference sequence. The partially aligned reads and unmapped reads with mapped mates (paired-end data) are used to build contigs with potential SV. 
- There are a number of aligners that soft-clip partially aligned sequences, bwa and Bowtie are two well-known tools:
  - [bwa](http://bio-bwa.sourceforge.net/)
  - [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml)

Configuration file
------------

- A configuration file with all the options is available, breakmer.cfg.
- The options in the configuration file determine where the output files are directed and stored as well as where key input files are located:
```
analysis_name=<sample_id>
targets_bed_file=<path to bed file containing locations of target regions>
sample_bam_file=<path to sample bam file>
normal_bam_file=<path to normal bam file>
analysis_dir=<path to analysis directory>
reference_data_dir=<path to where reference files will/are stored> 
other_regions_file=<path to bed file containing targeted unannotated cluster regions> 
cutadapt=<path to cutadapt binary> 
cutadapt_config_file=<path to cutadapt configuration file> 
jellyfish=<path to Jellyfish binary> 
blat=<path to BLAT binaries, blat, gfServer, gfClient>
reference_fasta=<path to whole genome reference fasta file>
gene_annotation_file=<path to ucsc_hg19_refgene.txt>
repeat_mask_file=<path to ucsc_hg19_rmsk.bed>
kmer_size=15
```

Input file formats
-----------

- targets_bed_file = tab-delimited file with columns: chr, start, end, region_name, feature_name
```
15      45003745        45003811        B2M     exon
15      45007621        45007922        B2M     exon
15      45008527        45008540        B2M     exon
```
- other_regions_file = tab-delimited file with columns: chr, start, end, region_name
   - This file is intended to cover regions that are not annotated in the annotation file. These are useful for cluster regions that are not well annotated in the annotation files.
```
14   22090057        23021075        TRA
7    141998851       142510972       TRB
```
- cutadapt_config_file = each row corresponds to a parameter for cutadapt (see cutadapt.cfg example file or cutadapt documentation)
  - The file provided is intended for data generated using the paired-end Illumina TruSeq library.
  - Many of the Illumina library sequences have been annotated [elsewhere](https://wikis.utexas.edu/display/GSAF/Illumina+-+all+flavors).


BreaKmer parameters
-------------
| Parameter | Description | Default |
|---------- | ----------- | ------- |
| -l, --log_level    | Logging level | Debug |
| -a, --keep_repeat_regions | Keep indels in repeat regions. Requires a repeat mask bed file. | False |
| -p, --preset_ref_data | Preset all the reference dta for all the target regions before running analysis. | False |
| -s, --indel_size | Indel size filter | 15 |
| -c, --trl_sr_thresh | Assembled read support threshold for translocations | 2 |
| -d, --indel_sr_thresh | Assembled read support threshold for indels | 5 | 
| -r, --rearr_sr_thresh | Assembled read support threshold for inversions and tandem duplications | 3 |
| -g, --gene_list | File containing a list of target region names to consider for analysis. Names must match targets_bed_file region names. | None |
| -k, --keep_intron_vars | Keep indels or rearrangements with breakpoints in intron regions | False |
| -v, --var_filter | Variant types to report (all, indel, trl, rearrangement) | all |
| -m, --rearr_min_seg_len | Threshold for minimum segment length to be rearranged | 30 |
| -n, --trl_min_seg_len | Threshold for minimum length of a translocation segment | 25 |
| -t, --align_thresh | Threshold for minimum read alignment for assembly | .90 |
| -z, --no_output_header | Suppress headers on output files. | False |

- kmer_size = option to change the length of the kmer size used (default = 15 bp).

Output files and formats
-----------
- While the program is running a log file analysis_dir/log.txt will continually be updated with information regarding the status of the analysis. 

### Logging

### Target data
  - For each region analyzed, a directory with the name of the region, as specified in the targets bed file, is created in a 'targets' directory (/<analysis_dir/>/targets).
  - Each target directory contains 'data', 'contigs', and 'kmers' directories. 
    - data - Contains the extracted reference-aligned reads from this target region as well as the kmers created from these extracted reads.
    - contigs - Contains a directory for each contig that was created, which contains the contig sequence in fasta format, and the reads used to assemble the contig in fastq format. The BLAT results are also stored with formatted output if a SV was called.
    - kmers - Contains the file with all the sample-only kmer sequences and how many reads in which it was contained as well as a file with the kmers and read ids that were used in assembling each contig.

### Final output 
  - When the program completes, the final output files are directed into a directory labeled 'output' within the specified analysis directory. 
  - Summary file
    - A summary file labeled /<analysis_name/>_summary.out contains columns: Target name, number of contigs assembled, total number of variants detected, number of indels, number of inversions and tandem duplications, number of translocations, and a list of translocation gene partners.
  - Output for each SV type are put in respective tab-delimited files, labeled /<analysis_name/>_/<indel,trl,inv_rearrangement,td_rearrangement/>_svs.out
    - The columns are:
       - Gene/Target names
       - Target genomic breakpoints (chr:pos)
       - Re-alignment CIGAR string
       - Number of re-alignment mismatches
       - Strand(s) of contig when realigned
       - Percentage re-alignment overlaps with repeat regions, length of aligned BLAT segments
       - Type of variation detected
       - Number of reads that cover the assembled contig at the inferred breakpoint
       - Number of kmers used to assemble the contig
       - Number of discordantly-mapped paired-end reads that support the event (Not applicable to indels).
       - Contig ID
       - Contig sequence. 


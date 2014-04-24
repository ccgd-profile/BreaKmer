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
  - Each "target" region can contain multiple subregions that are annotated by feature type (i.e., intron, exon). These feature types are used in the filtering steps with certain parameters. Note that the minimum coordinate and maximum coordinate +/- 200 bp are used as the boundaries for each region (e.g., 45003745-200,45008540+200 for B2M).
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
- reference_fasta = genome reference fasta formatted file containing all the chromosome reference sequences that were used to initially align the data.
- gene_annotation = Annotation file containing the location of reference genes. These can be downloaded from UCSC Genome Browser, (i.e., [hg19 refGene table](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/))
```
#bin    name    chrom   strand  txStart txEnd   cdsStart        cdsEnd  exonCount       exonStarts      exonEnds        score   name2   cdsStartStat    cdsEndStat      exonFrames
0       NM_032291       chr1    +       66999824        67210768        67000041        67208778        25      66999824,67091529,67098752,67101626,67105459,67108492,67109226,67126195,67
133212,67136677,67137626,67138963,67142686,67145360,67147551,67154830,67155872,67161116,67184976,67194946,67199430,67205017,67206340,67206954,67208755, 67000051,67091593,67098777,6710169
8,67105516,67108547,67109402,67126207,67133224,67136702,67137678,67139049,67142779,67145435,67148052,67154958,67155999,67161176,67185088,67195102,67199563,67205220,67206405,67207119,6721
0768,   0       SGIP1   cmpl    cmpl    0,1,2,0,0,0,1,0,0,0,1,2,1,1,1,1,0,1,1,2,2,0,2,1,1,
1       NM_032785       chr1    -       48998526        50489626        48999844        50489468        14      48998526,49000561,49005313,49052675,49056504,49100164,49119008,49128823,49
332862,49511255,49711441,50162984,50317067,50489434,    48999965,49000588,49005410,49052838,49056657,49100276,49119123,49128913,49332902,49511472,49711536,50163109,50317190,50489626,  0AGBL4    cmpl    cmpl    2,2,1,0,0,2,1,1,0,2,0,1,1,0,
1       NM_018090       chr1    +       16767166        16786584        16767256        16785385        8       16767166,16770126,16774364,16774554,16775587,16778332,16782312,16785336, 16767348,16770227,16774469,16774636,16775696,16778510,16782388,16786584, 0       NECAP2  cmpl    cmpl    0,2,1,1,2,0,1,2,
1       NM_052998       chr1    +       33546713        33585995        33547850        33585783        12      33546713,33546988,33547201,33547778,33549554,33557650,33558882,33560148,33
562307,33563667,33583502,33585644,      33546895,33547109,33547413,33547955,33549728,33557823,33559017,33560314,33562470,33563780,33583717,33585995,    0       ADC     cmpl    cmpl    -1
,-1,-1,0,0,0,2,2,0,1,0,2,
...
```
- repeat_mask_file = A BED formatted file containing repeat masked regions. These can be found from [UCSC Genome Table Browser](http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=370921603_HSUaLVPi7dEbtDqCy5W1ANaqC7Fz) with group = Repeats, track = RepeatMasker, table = rmsk
```
chr1    16777160        16777470        AluSp   2147    +
chr1    25165800        25166089        AluY    2626    -
chr1    33553606        33554646        L2b     626     +
chr1    50330063        50332153        L1PA10  12545   +
chr1    58720067        58720973        L1PA2   8050    -
chr1    75496180        75498100        L1MB7   10586   +
chr1    83886030        83886750        ERVL-E-int      980     -
chr1    100662895       100663391       L2a     1422    -
chr1    117440426       117440514       L1ME1   532     +
chr1    117440494       117441457       L1ME1   4025    +
...
```

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
  - While the program is running a log file will (\<analysis\_directory\>/log.txt) continually be updated with information regarding the status of the analysis.

### Target data
  - For each region analyzed, a directory with the name of the region, as specified in the targets bed file, is created in a 'targets' directory (\<analysis\_dir\>/targets).
  - Each target directory contains 'data', 'contigs', and 'kmers' directories. 
    - data - Contains the extracted reference-aligned reads from this target region as well as the kmers created from these extracted reads.
    - contigs - Contains a directory for each contig that was created, which contains the contig sequence in fasta format, and the reads used to assemble the contig in fastq format. The BLAT results are also stored with formatted output if a SV was called.
    - kmers - Contains the file with all the sample-only kmer sequences and how many reads in which it was contained as well as a file with the kmers and read ids that were used in assembling each contig.

### Final output 
  - When the program completes, the final output files are directed into a directory labeled 'output' within the specified analysis directory. 
  - Summary file
    - A summary file labeled \<analysis_name\>\_summary.out contains columns: Target name, number of contigs assembled, total number of variants detected, number of indels, number of inversions and tandem duplications, number of translocations, and a list of translocation gene partners.
  - Output for each SV type are put in respective tab-delimited files, labeled \<analysis_name\>\_\<indel,trl,inv\_rearrangement,td\_rearrangement\>\_svs.out
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
  - Each target gene in which a SV was detected has a separate output directory (\<analysis\_dir\>/output/\<target\_name\>) containing formatted output specific to the target and the related reference-aligned sequence reads for the contigs that contain the structural variants detected in BAM format.

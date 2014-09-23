BreaKmer
========

A method to identify genomic structural variation in target regions/genes from reference-aligned high-throughput sequence data. It uses a “kmer” strategy to assemble misaligned sequence reads for predicting insertions, deletions, inversions, tandem duplications, and translocations at base-pair resolution.

Installation
----------

The following are required for installation:
- [Python 2.7](https://www.python.org/download/releases/2.7)

Download the python scripts and run the command:
```
python setup.py install
```
Use appropriate commands for installing locally:
```
python setup.py install --user
```

Using the setup.py script for installation should setup the required python module dependencies for appropriate usage. If the install script is not used, the two modules need to be downloaded and installed if not already:
- [Biopython 1.62](http://biopython.org/wiki/Main_Page)
- [Pysam 0.6](https://code.google.com/p/pysam/)

BreaKmer has currently been installed and tested on:
- 64bit linux using CentOS release 5.5 and python2.7.2 
- 64bit linux using Ubuntu and python2.7.6

Usage
---------

List the available command line parameters.
```
python <PATH_TO_BREAKMER_DIR>/breakmer.py -h
```

Analyze all the target genes specified in the targets bed file.
```
python <path to breakmer scripts>/breakmer.py <options> <path to breakmer configuration file>
```

Analyze a subset of genes specified in a file.
```
python <path to breakmer scripts>/breakmer.py -g <file containing list of target genes to analyze> <path to breakmer configuration file>
```

Requirements
---------

### Programs
- BLAT standalone and server binaries ([blat, gfServer, gfClient, faToTwoBit](http://hgdownload.cse.ucsc.edu/admin/exe/)).
  - Re-alignment to reference sequence.
  - Versions tested :
    - standalone BLAT v35x1
    - gfServer v35x1
    - gfClient v35x1
- [Cutadapt](https://code.google.com/p/cutadapt/)
  - Trims adapter sequence from aligned reads.
  - v1.5 tested
- [Jellyfish](http://www.cbcb.umd.edu/software/jellyfish)
  - Generating kmers.
  - [v1.1.11](http://www.cbcb.umd.edu/software/jellyfish/jellyfish-1.1.11.tar.gz) tested
  - [v2.1.3](http://www.genome.umd.edu/jellyfish.html) untested

When these programs are installed, the paths to the binaries can be either be specified in the BreaKmer configuration file or put in the path (e.g. export PATH=$PATH:/path/to/binary).

### Sequence data
- Sample bam file
  - BreaKmer requires sequence reads that have been aligned to a reference sequence in a binary alignment format (BAM). The alignment program must soft-clip or trim reads that can be partially aligned to the reference sequence (e.g., bwa, bowtie, novoalign, mosaik). The partially aligned reads and unmapped reads with mapped mates (paired-end data) are used to build contigs with potential SV. This bam file is required as input in the configuration file as the "sample_bam_file".
    - The BAM files need to be sorted and indexed, with the indexed files in the same directory as the BAM file.
    - There are a number of aligners that soft-clip partially aligned sequences, bwa and Bowtie are two well-known tools:
      - [bwa](http://bio-bwa.sourceforge.net/)
      - [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml)
    - It is also useful to mark duplicate reads prior to using BreaKmer, as these will inflate read counts for identified variants.
      - [Picard MarkDuplicates](http://picard.sourceforge.net/command-line-overview.shtml#MarkDuplicates)
- Normal bam file
  - A matched normal with similarly targeted sequencing data as the sample can be input to help filter germline events. The sequences in the normal bam file are processed for each target in a similar manner to the sample sequences and used to further filter our kmers. This is an optional input in the configuration file as "normal_bam_file".

### Reference data
- Reference fasta file 
  - BreaKmer makes use of reference sequence data throughout the program. The genomic reference sequence used to align the short sequence reads is required as an input in the configuration file as the "reference_fasta".
  - Format: single fasta file with chromosome number/id as names (i.e. '>1', '>2', '>3')
  - Hg19 fasta files can be downloaded from UCSC Genome Browser(http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/)
  - This file should be placed in a writeable directory. A 2bit file will be generated from this file to start the blat server.
- Reference genome annotation 
  - This file is used by BreaKmer to annotate the locations of the breakpoints identified and the genes involved. This is required as input in the configuration file as the "gene_annotation_file".
  - Format: tab delimited file containing a row for each RefSeq transcript with multiple columns describing the coding coordinates of the transcript.
  - Hg19 RefSeq annotation file can be downloaded from UCSC Genome Browser [Table Browser](http://genome.ucsc.edu/cgi-bin/hgTables?command=start) using assembly: Feb. 2009 (GRCh37/hg19), group: Genes and Gene Predictions, track: RefSeq Genes, table: refGene
- Repeat mask bed file 
  - This file is used by BreaKmer to determine whether breakpoints and identified variants lie within simple and low-complexity repeat regions that have been annotated. This is an optional input in the configuration file as "repeat_mask_file".
  - Format: tab delimited file containing the coordinates for the various repeat regions.
  - Hg19 repeat regions bed file can be downloaded from UCSC Genome Browser [Table Browser](http://genome.ucsc.edu/cgi-bin/hgTables?command=start) with group = Repeats, track = RepeatMasker, table = rmsk, output format: BED
- Alternate reference assembly sequences 
  - This is an optional input in the configuration file as the "alternate_fastas". These are used to supplement the required reference fasta sequence and helpful for removing indel variants that exist in the reference sequence but not either of the alternate assembly sequences.
  - Format: single fasta file with chromosome number/id as names (i.e. '>1', '>2', '>3')
  - There are two alternate assemblies available for human, [CHM1_1.1](http://www.ncbi.nlm.nih.gov/assembly/GCF_000306695.2/) and [HuRef](http://www.ncbi.nlm.nih.gov/assembly/GCA_000002125.2).

### Other files
- Target regions bed file
  - The target regions sequenced need to be specified in a bed file for BreaKmer to know the coordinates of the regions sequenced and which are to be analyzed. This is a required input in the configuration file as the "targets_bed_file".
  - Format: a tab delimited file containing the chromosome, start, end, HUGO gene name, and exon/intron feature of a tile region of the genome. The tiled regions should not be overlapping.

Configuration file
------------

- A template configuration file with all the options is available, breakmer.cfg.
- Below lists and describes the parameters that can be used in the configuration file. Do not keep commented text (i.e., #Required parameters) in the configuration file.
- Note that the paths to the six required program binaries (Cutadapt, Jellyfish, blat, gfServer, and gfClient) can be set in the configuration file
  or these binaries can be included in the users path (e.g., for linux users: export PATH=$PATH:/path/to/binary).
- Use full paths (e.g., /home/bob and not ~/)
```
# Required parameters
analysis_name=<sample_id, string value that all the output files will contain>
targets_bed_file=<path to bed file containing locations of target regions>
sample_bam_file=<path to sorted and indexed sample bam file>
analysis_dir=<path to analysis directory, location where analysis and final output files will be located>
reference_data_dir=<path to where reference files will/are stored for each of the targeted genes analyzed> 
cutadapt=<path to cutadapt binary v1.5, i.e. /usr/bin/cutadapt-1.5/bin/cutadapt> 
cutadapt_config_file=<path to cutadapt configuration file> 
jellyfish=<path to Jellyfish binary v1.1.11 required, i.e. /usr/bin/jellyfish>
blat=<path to blat binary, i.e. /usr/bin/blat>
gfserver=<path to gfServer binary, i.e. /usr/bin/gfServer>
gfclient=<path to gfClient binary, i.e. /usr/bin/gfClient>
fatotwobit=<path to faToTwoBit binary, i.e. /usr/bin/faToTwoBit>
reference_fasta=<path to whole genome reference fasta file, one file with all records>
gene_annotation_file=<path to gene annotation file, e.g., ucsc_hg19_refgene.txt>
kmer_size=15

# Optional parameters
other_regions_file=<path to bed file containing coordinates for targeted unannotated cluster regions if they exist, such as IGH, IGK> 
normal_bam_file=<path to normal bam file, can be used to filter germline events with matched-normal sample>
alternate_fastas=<comma delimited list of the paths to alternate fasta files, such as HuRef or CHM1>
repeat_mask_file=<path to ucsc_hg19_rmsk.bed> # Only used when available, useful for helping filtering events in simple repeat regions.
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
| -a, --keep_repeat_regions | Keep indels in repeat regions. No repeat mask bed file required if set. | False |
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

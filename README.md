BreaKmer
========

A method to identify genomic structural variation in target regions/genes from reference-aligned high-throughput sequence data. It uses a “kmer” strategy to assemble misaligned sequence reads for predicting insertions, deletions, inversions, tandem duplications, and translocations at base-pair resolution.

http://a-bioinformatician.github.io/BreaKmer

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
- 64bit linux using Ubuntu release 14.04 and python2.7.6

Usage
---------

List the available command line parameters.
```
python <PATH_TO_BREAKMER_DIR>/breakmer.py -h
```

Analyze all the target genes specified in the targets bed file.
```
python <path to breakmer scripts>/breakmer.py run <options> -c <path to breakmer configuration file>
```

Prepare the reference data files before starting the analysis.
```
python <path to breakmer scripts>/breakmer.py prepare_reference_data <options> -c <path to breakmer configuration file>
```

Start the blat server on a specific host and port.
```
python <path to breakmer scripts>/breakmer.py start_blat_server <options> -c <path to breakmer configuration file>
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
  - [v2.1.3](http://www.genome.umd.edu/jellyfish.html) tested

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
  - This file is used by BreaKmer to annotate the locations of the breakpoints identified and the genes involved. This is required as input in the configuration file as the "gene_annotation_file" only when annotation is desired for the output files (currently not implemented) or when the --generate_image option is set and the corresponding genes and transcripts of the variant are to be visualized.
  - Format: GTF formatted file (see below).
  - Gencode transcript annotation file can be downloaded from [Gencode Release 19 (GRCh37.p13)](http://www.gencodegenes.org/releases/19.html).

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
jellyfish=<path to Jellyfish binary, i.e. /usr/bin/jellyfish>
blat=<path to blat binary, i.e. /usr/bin/blat>
gfserver=<path to gfServer binary, i.e. /usr/bin/gfServer>
gfclient=<path to gfClient binary, i.e. /usr/bin/gfClient>
fatotwobit=<path to faToTwoBit binary, i.e. /usr/bin/faToTwoBit>
reference_fasta=<path to whole genome reference fasta file, one file with all records>
gene_annotation_file=<path to gene annotation file>
kmer_size=15

# Optional parameters
other_regions_file=<path to bed file containing coordinates for targeted unannotated cluster regions if they exist, such as IGH, IGK> 
normal_bam_file=<path to normal bam file, can be used to filter germline events with matched-normal sample>
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
  - This shold be a single file containing all the sequences.
  - Note that the chromosome IDs should match the chromosome IDs of the alignment file (i.e., using chr1 vs. 1).
- gene_annotation = Annotation file containing the genomic coordinates for annotated gene transcripts and their corresponding exon regions. These can be downloaded from Gencode, (i.e., [Gencode Release 19 (GRCh37.p13)](http://www.gencodegenes.org/releases/19.html))
```
column-number content values/format
1 chromosome name chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,M}
2 annotation source {ENSEMBL,HAVANA}
3 feature-type  {gene,transcript,exon,CDS,UTR,start_codon,stop_codon,Selenocysteine}
4 genomic start location  integer-value (1-based)
5 genomic end location  integer-value
6 score (not used)  .
7 genomic strand  {+,-}
8 genomic phase (for CDS features)  {0,1,2,.}
9 additional information as key-value pairs see below

chr21   HAVANA  transcript      10862622        10863067        .       +       .       gene_id "ENSG00000169861"; transcript_id "ENST00000302092"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "IGHV1OR15-5"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "IGHV1OR15-5-001"; level 2; havana_gene "OTTHUMG00000074130"; havana_transcript "OTTHUMT00000157419";
chr21   HAVANA  exon    10862622        10862667        .       +       .       gene_id "ENSG00000169861"; transcript_id "ENST00000302092"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "IGHV1OR15-5"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "IGHV1OR15-5-001"; level 2; havana_gene "OTTHUMG00000074130"; havana_transcript "OTTHUMT00000157419";
chr21   HAVANA  CDS     10862622        10862667        .       +       0       gene_id "ENSG00000169861"; transcript_id "ENST00000302092"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "IGHV1OR15-5"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "IGHV1OR15-5-001"; level 2; havana_gene "OTTHUMG00000074130"; havana_transcript "OTTHUMT00000157419";
chr21   HAVANA  start_codon     10862622        10862624        .       +       0       gene_id "ENSG00000169861"; transcript_id "ENST00000302092"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "IGHV1OR15-5"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "IGHV1OR15-5-001"; level 2; havana_gene "OTTHUMG00000074130"; havana_transcript "OTTHUMT00000157419";
chr21   HAVANA  exon    10862751        10863067        .       +       .       gene_id "ENSG00000169861"; transcript_id "ENST00000302092"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "IGHV1OR15-5"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "IGHV1OR15-5-001"; level 2; havana_gene "OTTHUMG00000074130"; havana_transcript "OTTHUMT00000157419";
chr21   HAVANA  CDS     10862751        10863064        .       +       2       gene_id "ENSG00000169861"; transcript_id "ENST00000302092"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "IGHV1OR15-5"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "IGHV1OR15-5-001"; level 2; havana_gene "OTTHUMG00000074130"; havana_transcript "OTTHUMT00000157419";
chr21   HAVANA  stop_codon      10863065        10863067        .       +       0       gene_id "ENSG00000169861"; transcript_id "ENST00000302092"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "IGHV1OR15-5"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "IGHV1OR15-5-001"; level 2; havana_gene "OTTHUMG00000074130"; havana_transcript "OTTHUMT00000157419";
chr21   HAVANA  UTR     10863065        10863067        .       +       .       gene_id "ENSG00000169861"; transcript_id "ENST00000302092"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "IGHV1OR15-5"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "IGHV1OR15-5-001"; level 2; havana_gene "OTTHUMG00000074130"; havana_transcript "OTTHUMT00000157419";
...
```

BreaKmer parameters
-------------
| Program Mode | Parameter | Description | Default |
| -------- |---------- | ----------- | ------- |
| run | | Perform structural variant detection analysis | NA |
| | -h, --help | Show the help message and exit | NA |
| | --log_level | The level of logging detail to record while program is running. | DEBUG |
| | --indel_size | The minimum indel size (in base pairs) to keep and output. | 15 |
| | --trl_sr_thresh | Split read support threshold for translocation events (i.e., the minimum number of split reads required to keep a translocation event). | 2 |
| | --indel_sr_thresh | Split read support threshold for indels | 5 |
| | --rearr_sr_thresh | Split read support threshold for general rearrangements, inversions, tandem duplications, or other unclassified rearrangements. | 2 |
| | --rearr_min_seg_len | Minimum length (in base pairs) of a realignment portion of an assembled contig sequence to consider for filtering a rearrangement event (inversions, tandem duplications, other unclassified rearrangements). | 30 |
| | --trl_min_seg_len | Minimum length (in base pairs) of a realignment portion of an assembled contig sequence to consider for filtering a translocation event. | 25 |
| | --align_thresh | Threshold for minimum read alignment during assembly. | 0.9 |
| | --no_output_header | Suppress column headers in output files. | False |
| | --discread_only_thresh | The number of discordant read pairs in a cluster to output without evidence from a split read event. | 2 |
| | --generate_image | Generate a plot containing the assembled reads used to create a contig sequence with a called variant along with genes and transcripts associated with the event - requires gene_annotation_file in configuration file to be set with GTF file containing transcript annotations. | False |
| | --hostname | The hostname for the blat server. | localhost |
| | -g, --gene_list | File containing the names of specific regions to analyze, these target names must match the names in the bed file defining the target regions to analyze. | None |
| | -f, --filter_list | A file containing specific variant events to filter out. | None |
| | -n, --nprocessors | The number of processors to use for analysis. | 1 |
| | -s, --start_blat_sever | Option to indicate that a blat server needs to be started before analysis begins. A random part number and the localhost will be used if neither is specified. | False |
| | -k, --keep_blat_server | Option to keep the blat server running after the analysis completes. | False |
| | -p, --port_number | The port number for the blat server to either start on or already is running on. | None |
| | -c, --config | The configuration filename that contains additional parameters. | None |
| | | | |
| start_blat_server | | Start the blat server prior to performing the analysis | NA |
| | -h, --help | Show the help message and exit. | NA |
| | --hostname | The hostname for the blat server. | localhost |
| | -p, --port_number | The port number for the blat server to either start on or already is running on. | None |
| | -c, --config | The configuration filename that contains additional parameters. | None |
| | | | |
| prepare_reference_data | | Extract the reference sequence from all the input target regions and store them in files to access during analysis. | NA |
| | -h, --help | Show the help message and exit | NA |
| | -g, --gene_list | File containing the names of specific regions to analyze, these target names must match the names in the bed file defining the target regions to analyze. | None |
| | -c, --config | The configuration filename that contains additional parameters. | None |
| | -n, --nprocessors | The number of processors to use for analysis. | 1 |

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
       - Target_Name - Target region name
       - SV_type - Rearrangement type (rearrangement, indel)
       - SV_subtype - Rearrangement subtype (trl, inversion, tandem duplication)
       - Description - Description of the indel (I5 = Insertion of 5 base pairs, D5 = deletion of 5 base pairs)
       - All_genomic_breakpoints - All the genomic breakpoints for the structural variant event.
       - Target_genomic_breakpoints - Breakpoints only within the target region.
       - Split_read_counts - Number of assembled reads that cover the assembled contig at the inferred breakpoint.
       - Discordant_read_counts - Number of discordantly-mapped paired-end reads that support the event (Not applicable to indels).
       - Read_depth_at_genomic_breakpoints - Depth of coverage at the inferred breakpoints.
       - Align_cigar - Cigar formatted string to indicate realignment to the reference sequence.
       - Strands - Strands the contig sequence realigned to on the reference sequence.
       - Total_mismatches - Number of base pairs that were mismatched in the realignment to the reference sequence.
       - Realignment_uniqueness - Uniqueness metric for the realignment.
       - Contig_ID - Contig ID
       - Contig_length - Number of base pairs in the contig sequence.
       - Contig_sequence - Contig sequence.
       - Filtered - Boolean field to indicate if the variant passed basic filtering criteria.
       - Filtered_reason - Message for filter reason if Filtered is True.
       - Filter_values - Additional information from the realignment.
  - Each target gene in which a SV was detected has a separate output directory (\<analysis\_dir\>/output/\<target\_name\>) containing formatted output specific to the target and the related reference-aligned sequence reads for the contigs that contain the structural variants detected in BAM format.

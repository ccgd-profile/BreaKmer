BreaKmer
========

A method to identify structural variation from target-captured NGS data.

Installation
============

Download the python scripts and run the command:
python setup.py install

Use appropriate commands for installing locally:
python setup.py install --user

Once installed, BreaKmer can be run with hte following command:
python breakmer.py <options> <path to config file>

Note that BreaKmer requires the following to operate properly:
- BLAT standalone and server (gfServer, gfClient) (http://hgdownload.cse.ucsc.edu/admin/exe/)
- Cutadapt - used for trimming adapter sequence from aligned reads (https://code.google.com/p/cutadapt/)
- Jellyfish - used for generating kmers (http://www.cbcb.umd.edu/software/jellyfish/)

Paths to these binaries need to be specified in the configuration file input to breakmer.py

Configuration file
==================

- A configuration file with all the options is available, kmer_region.cfg.
- The options in the configuration file determine where the output files are directed and stored as well as where key input files are located.


Input file formats
==================

- targets_bed_file = tab-delimited file with columns <chr,start,end,region_name,feature_name>
   - example row = 1       2489165 2489273 TNFRSF14        exon
- other_regions_file = tab-delimited file with columns <chr,start,end,region_name>
   - This file is intended to cover regions that are not annotated in the annotation file.
- cutadapt_config_file = each row corresponds to a parameter for cutadapt (see cutadapt.cfg example file or cutadapt documentation)
 

#! /usr/bin/local/python

import os
import logging
import shutil
from Bio import SeqIO
import subprocess
import random
import time
from pysam import *
import multiprocessing
from itertools import izip, islice, repeat, izip_longest


################################################################
# Utility functions
###############################################################

def calc_contig_complexity(seq, N=3, w=6) :
  cmers = [] 
  for i in range(len(seq)) :
    s = max(0,i-w)
    e = min(len(seq),i+w)
    cmer = count_nmers(seq[s:e], N)
    n = len(cmer)
    nmod = float(n) / float(e-s)
    cmers.append(round(nmod,2))
  cmers_mean = sum(cmers) / len(cmers)
  return cmers_mean, cmers

def count_nmers(seq, N) :
  nmers = {}
  total_possible = len(seq) - 2
  for i in range(len(seq) - (N - 1)) :
    mer = str(seq[i:i+N]).upper()
    if mer not in nmers : nmers[mer] = 0 
    nmers[mer] += 1
  return nmers
#-----------------------------------------------------------
def is_number(s):
  try:
    float(s)
    return True
  except ValueError:
    pass
 
  try:
    import unicodedata
    unicodedata.numeric(s)
    return True
  except (TypeError, ValueError):
    pass
  return False
#-----------------------------------------------------------

#-----------------------------------------------------------
def filter_by_feature(brkpts, query_region, keep_intron_vars) :
  in_filter = False
  span_filter = False
  if not keep_intron_vars :
    in_vals, span_vals = check_intervals(brkpts, query_region)
    if in_vals[0] : 
      if 'exon' not in in_vals[1] : 
        in_filter = True 
    else :
      in_filter = True
    if span_vals[0] :
      if 'exon' not in span_vals[1] :
        span_filter = True
    else :
      span_filter = True
  return in_filter, span_filter
#-----------------------------------------------------------

#-----------------------------------------------------------
def check_intervals(breakpts, query_region ) :
  in_values = [False, [], []]
  span_values = [False, [], []]
  bp_in_interval = False
  contains_interval = False
  in_features = []
  in_interval = None
  for bp in breakpts :
    for interval in query_region[4] :
      if (int(bp) >= (interval[1]-20)) and (int(bp) <= (interval[2]+20)) : 
        in_values[1].append(interval[4])
        in_values[0] = True
        in_values[2].append(interval)
      if (interval[2] <= max(breakpts)) and (interval[1] >= min(breakpts)) :
        span_values[0] = True
        span_values[1].append(interval[4])
        span_values[2].append(interval)
  return in_values, span_values
#-----------------------------------------------------------

#-----------------------------------------------------------
def setup_logger(param_opts, name) :
  output_path = os.path.abspath(os.path.normpath(param_opts['analysis_dir']))
  if not os.path.exists(output_path) : os.makedirs(output_path)

  logger = logging.getLogger(name)
  logger.setLevel(logging.DEBUG)

  # FileHandler
  fh = logging.FileHandler(os.path.join(output_path,'log.txt'),mode='w')
  fh.setLevel(logging.DEBUG)
  # ConsoleHandler
  ch = logging.StreamHandler()
  ch.setLevel(logging.ERROR)

  formatter = logging.Formatter(fmt='%(asctime)s - %(name)s - %(levelname)s - %(message)s',datefmt='%m/%d/%Y %I:%M:%S %p')
  fh.setFormatter(formatter)
  ch.setFormatter(formatter)

  logger.addHandler(fh)
  logger.addHandler(ch)
#-----------------------------------------------------------

#-----------------------------------------------------------
def check_repeat_regions(coords, repeat_locs) :
  start, end = coords
  seg_len = float(end-start)
  in_repeat = False
  rep_overlap = 0.0
  rep_coords = []
  filter_reps_edges = [False, False]
  for rloc in repeat_locs :
    rchr, rbp1, rbp2, rname = rloc
    if (rbp1 >= start and rbp1 <= end) or (rbp2 >= start and rbp2 <= end) or (rbp1 <= start and rbp2 >= end):
      in_repeat = True
      rep_overlap += float(min(rbp2,end)-max(rbp1,start))
      rep_coords.append((rbp1,rbp2))
      # Simple or low complexity seq repeat for filtering
      if rname.find(")n") > -1 or rname.find("_rich") > -1 : 
        if (rbp1<=start and rbp2>=start) : filter_reps_edges[0] = True
        elif (rbp1<=end and rbp2>=end) : filter_reps_edges[1] = True
#      if rep_overlap >= seg_len :
#        break
  roverlap = round( (float(min(rep_overlap, seg_len)) / float(seg_len))*100, 2)
  #break
  return in_repeat, roverlap, rep_coords, filter_reps_edges
#-----------------------------------------------------------

#-----------------------------------------------------------
def get_marker_fn(fn) :
  return os.path.join(os.path.split(fn)[0],"."+os.path.basename(fn))
#-----------------------------------------------------------

#-----------------------------------------------------------
def run_jellyfish(fa_fn, jellyfish, kmer_size) :
  logger = logging.getLogger('root')
  file_path = os.path.split(fa_fn)[0]
  file_base = os.path.basename(fa_fn)
  dump_fn = os.path.join(file_path,file_base + "_" + str(kmer_size) + "mers_dump")
  dump_marker_fn = get_marker_fn(dump_fn)
  if not os.path.isfile(dump_marker_fn) :
    count_fn = os.path.join(file_path, file_base + "_" + str(kmer_size) + "mers_counts")
    logger.info('Running %s on file %s to determine kmers'%(jellyfish,fa_fn)) 
    cmd = '%s count -m %d -s %d -t %d -o %s %s'%(jellyfish,kmer_size,100000000,8,count_fn,fa_fn)
    logger.info('Jellyfish counts system command %s'%cmd)
    p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
    output, errors = p.communicate()
    logger.info('Jellyfish count output %s'%output)
    logger.info('Jellyfish count errors %s'%errors)
    cmd = '%s dump -c -o %s %s'%(jellyfish,dump_fn,count_fn+"_0")
    logger.info('Jellyfish dump system command %s'%cmd)
    p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
    output, errors = p.communicate()
    logger.info('Jellyfish dump output %s'%output)
    logger.info('Jellyfish dump errors %s'%errors)
    cmd = 'touch %s'%dump_marker_fn
    p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
    output, errors = p.communicate()  
    logger.info('Completed jellyfish dump %s, touching marker file %s'%(dump_fn,dump_marker_fn))
  else :
    logger.info('Jellfish already run and kmers already generated for target.')
  return dump_fn
#-----------------------------------------------------------

#-----------------------------------------------------------
def setup_ref_data(setup_params) :
   genes = setup_params[0]
   rep_mask, ref_fa, ref_path, jfish_path, kmer_size = setup_params[1]

   logger = logging.getLogger('root')

   for gene in genes :
     chr, bp1, bp2, name, intvs = gene
     gene_ref_path = os.path.join(ref_path,name)
     if rep_mask : 
       logger.info('Extracting repeat mask regions for target gene %s.'%name)
       setup_rmask(gene, gene_ref_path, rep_mask)
     
     logger.info('Extracting refseq sequence for %s'%name)
     f_fn = extract_refseq_fa(gene, gene_ref_path, ref_fa, "forward")
     r_fn = extract_refseq_fa(gene, gene_ref_path, ref_fa, "reverse")
     run_jellyfish(f_fn, jfish_path, kmer_size)
     run_jellyfish(r_fn, jfish_path, kmer_size)
#-----------------------------------------------------------

#-----------------------------------------------------------
def get_fastq_reads(fn, sv_reads) :
  read_len = 0
  filtered_fq_fn = fn.split(".fastq")[0] + "_filtered.fastq"
  filt_fq = open(filtered_fq_fn, 'w')
  fq_recs = {}
#  f = open(fn,'r')
#  fq_recs = list(SeqIO.parse(f,'fastq'))
  for header,seq,qual in FastqFile(fn) : 
    qname_split = header.lstrip("@").split("_")
    indel_only = qname_split[-1]
    qname = "_".join(qname_split[0:len(qname_split)-1])
    if qname in sv_reads : 
      oseq, sc_seqs, clip_coords, indel_meta = sv_reads[qname]
      cleaned_seq = seq 
      old_seq = oseq.seq
      add = True
      if str(cleaned_seq) != str(old_seq) and sc_seqs :
        sc_clips = sc_seqs['clipped']
        idx = old_seq.find(cleaned_seq)
        trimmed_seq = ''
        if idx == 0 : 
          trimmed_seq = old_seq[len(cleaned_seq):len(old_seq)]
        else : trimmed_seq = old_seq[0:idx]    
        sc_lens = 0
        for sc_seq in sc_clips : 
          sc_lens += len(sc_seq)
          if trimmed_seq.find(sc_seq) > -1 : 
            add = False
        if len(cleaned_seq) == (len(old_seq) - sc_lens) :
          for sc_seq in sc_clips : 
            if cleaned_seq.find(sc_seq) == -1 :
              # Don't add, just trimmed clipped portion.
              add = False
#    else : print qname, 'not in sv reads'
    if add :
      filt_fq.write(header + "\n" + seq + "\n+\n" + qual + "\n") 
      fr = fq_read(header, seq, qual, indel_meta)
      read_len = max(read_len, len(fr.seq)) 
      seq = fr.seq
      if seq not in fq_recs : 
        fq_recs[seq] = []
      fq_recs[seq].append(fr)
  filt_fq.close()
  return filtered_fq_fn, fq_recs, read_len
#-----------------------------------------------------------

#-----------------------------------------------------------
def get_fastq_reads_old(fn, sv_reads) :
  read_len = 0
  fq_recs = {}
  f = open(fn,'r')
#  fq_recs = list(SeqIO.parse(f,'fastq'))
  for header,seq,qual in FastqFile(fn) : 
    qname = header.lstrip("@")
    if qname in sv_reads : 
      oseq, sc_seqs, clip_coords = sv_reads[qname]
      cleaned_seq = seq 
      old_seq = oseq.seq
      add = True
      if str(cleaned_seq) != str(old_seq) and sc_seqs :
        idx = old_seq.find(cleaned_seq)
        trimmed_seq = ''
        if idx == 0 : 
          trimmed_seq = old_seq[len(cleaned_seq):len(old_seq)]
        else : trimmed_seq = old_seq[0:idx]    
        sc_lens = 0
        for sc_seq in sc_seqs : 
          sc_lens += len(sc_seq)
          if trimmed_seq.find(sc_seq) > -1 : 
            add = False
        if len(cleaned_seq) == (len(old_seq) - sc_lens) :
          for sc_seq in sc_seqs : 
            if cleaned_seq.find(sc_seq) == -1 :
              # Don't add, just trimmed clipped portion.
              add = False
#    else : print qname, 'not in sv reads'
    if add :
      fr = fq_read(header, seq, qual)
      read_len = max(read_len, len(fr.seq))   
      fq_recs[fr.id] = fr
  return fq_recs, read_len
#-----------------------------------------------------------

#-----------------------------------------------------------
def load_kmers(fns,kmers) :
  fns = fns.split(",")
  for fn in fns :
    f = open(fn,'rU')
    for line in f.readlines() :
      line = line.strip()
      mer, count = line.split()
      if not mer in kmers : kmers[mer] = 0
      kmers[mer] += int(count)
  return kmers
#-----------------------------------------------------------

#-----------------------------------------------------------
# Store all repeats in a dict by chrom
#-----------------------------------------------------------
def setup_rmask_all(rmask_fn) :
  logger = logging.getLogger('root')

  rmask = {}
  f = open(rmask_fn,'rU')
  flines = f.readlines()
  for line in flines :
    line = line.strip()
    rchr,rbp1,rbp2,rname = line.split("\t")[0:4]
    rchr = rchr.replace('chr','')
    if rchr not in rmask : rmask[rchr] = []
    rmask[rchr].append((rchr,int(rbp1),int(rbp2),rname))
  return rmask
#-----------------------------------------------------------

#-----------------------------------------------------------
# Stores repeats for a specific gene
#-----------------------------------------------------------
def setup_rmask(gene_coords, ref_path, rmask_fn) :
  logger = logging.getLogger('root')
  chrom, s, e, name, intervals = gene_coords
  mask_out_fn = os.path.join(ref_path,name+'_rep_mask.bed')
  marker_fn = get_marker_fn(mask_out_fn)

  rmask = []
  if not os.path.isfile(marker_fn) :
    fout = open(mask_out_fn,'w')
    f = open(rmask_fn,'rU')
    flines = f.readlines()
    for line in flines :
      line = line.strip()
      rchr,rbp1,rbp2,rname = line.split("\t")[0:4]
      rchr = rchr.replace('chr','')
      if rchr == chrom :
        if int(rbp1) >= int(s) and int(rbp2) <= int(e) : 
          fout.write("\t".join([str(x) for x in [rchr,int(rbp1),int(rbp2),rname]])+"\n")
          rmask.append((rchr,int(rbp1),int(rbp2),rname))
    f.close()
    fout.close()

    cmd = 'touch %s'%marker_fn
    p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
    output, errors = p.communicate()  
    logger.info('Completed writing repeat mask file %s, touching marker file %s'%(mask_out_fn,marker_fn))
  else :
    rep_f = open(mask_out_fn,'rU')
    rep_flines = rep_f.readlines()
    for line in rep_flines :
      line = line.strip()
      rchr,rbp1,rbp2,rname = line.split()
      rmask.append((rchr,int(rbp1),int(rbp2),rname))
  return rmask
#-----------------------------------------------------------

#-----------------------------------------------------------
def extract_refseq_fa(gene_coords, ref_path, ref_fa, direction) :
  logger = logging.getLogger('root')

  chrom, s, e, name, intervals = gene_coords
  fa_fn = os.path.join(ref_path,name+'_'+direction+'_refseq.fa')
  marker_fn = get_marker_fn(fa_fn) 

  if not os.path.isfile(marker_fn) :
    ref_d = SeqIO.to_dict(SeqIO.parse(ref_fa, 'fasta'))
    seq_str = ''
    seq = ref_d[chrom].seq[(s-200):(e+200)]
    if direction == "reverse" :
      seq_str = str(seq.reverse_complement())
    else :
      seq_str = str(seq)
    fa = open(fa_fn,'w')
    fa.write(">"+name+"\n"+seq_str+"\n")
    fa.close()
    cmd = 'touch %s'%marker_fn
    p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
    output, errors = p.communicate()
    logger.info('Completed writing refseq fasta file %s, touching marker file %s'%(fa_fn, marker_fn))
  else :
    logger.info('Refseq sequence fasta (%s) exists already'%fa_fn)
  return fa_fn 
#-----------------------------------------------------------

#-----------------------------------------------------------  
def seq_trim(qual_str, min_qual) :
  counter = 0
  while ord(qual_str[counter])-33 < min_qual :
    counter += 1
    if counter == len(qual_str) : break
  return counter  
#----------------------------------------------------------- 

#----------------------------------------------------------- 
def get_seq_readname(read) :
  end = '1'
  if read.is_read2 : end = '2'
  return read.qname + "/" + end
#----------------------------------------------------------- 

#----------------------------------------------------------- 
def trim_coords(qual_str, min_qual):
  q = []
  coords = [0,len(qual_str)]
  start = seq_trim(qual_str, min_qual)
  if start == len(qual_str) :
    return (0,0,0)
  else :
    end = len(qual_str) - seq_trim(qual_str[::-1], min_qual)
    lngth = end - start
    return (start,end,lngth)
#----------------------------------------------------------- 

#----------------------------------------------------------- 
def trim_qual(read, min_qual, min_len):
#  print read
#  print read.seq
#  print read.qual
  qual_str = read.qual
  q = []
  coords = [0,len(qual_str)]
  start = seq_trim(qual_str, min_qual)
  if start == len(qual_str) :
    return None
  else :
    end = len(qual_str) - seq_trim(qual_str[::-1], min_qual)
    lngth = end - start
    if lngth < min_len : return None
    nseq = read.seq[start:end]
    nqual = qual_str[start:end]
    read.seq = nseq
    read.qual = nqual
    return read
#----------------------------------------------------------- 

#----------------------------------------------------------- 
def fq_line(read, indel_only, min_len, trim=True) :
  add_val = '0'
  if indel_only : add_val = '1'
  lineout = None
  if trim : read = trim_qual(read, 5, min_len)
  if read : 
    lineout = "@"+get_seq_readname(read) + "_" + add_val + "\n" + read.seq + "\n+\n" + read.qual + "\n"
  return lineout
#----------------------------------------------------------- 

#----------------------------------------------------------- 
def get_overlap_index_nomm(a,b) :
  i = 0
  while a[i:] != b[:len(a[i:])]:
    i += 1
  return i
#----------------------------------------------------------- 

#----------------------------------------------------------- 
def seq_complexity(seq, N) :
  nmers = {}
  total_possible = len(seq) - 2
  for i in range(len(seq) - (N - 1)) :
    nmers[str(seq[i:i+N]).upper()] = True
  complexity = round((float(len(nmers))/float(total_possible))*100,4)
  return complexity
#-----------------------------------------------------------

#-----------------------------------------------------------
def get_overlap_index_mm(a, b) :
  i = 0
  nmismatch = [0,0]
  match = False
  while i < len(a) and not match :
    nmismatch = [0,0]
    c = 0
    match_len = min(len(a[i:]), len(b[:len(a[i:])]))
    for aa,bb in zip(a[i:],b[:len(a[i:])]) :
      if aa != bb : 
        nmismatch[0] += 1
        nmismatch[1] += 1
      else :
        nmismatch[0] = 0
      if nmismatch[0] > 1 or nmismatch[1] > 3 : 
        break  
      c += 1
#    print c, match_len, i, c== match_len, a[i:], b[:len(a[i:])]
    if c == match_len : 
      match = True
    i += 1
  return i-1
#-----------------------------------------------------------

#----------------------------------------------------------- 
def get_read_kmers(seq,l,skmers) :
  kmers = []
  i = 0
  while (i+l) <= len(seq) :
    k = seq[i:i+l]
    #if k in skmers :
    kmers.append(k)
    i += 1
  return list(set(kmers)&set(skmers))
#----------------------------------------------------------- 

#----------------------------------------------------------- 
def get_overlap_index(a,b) :
  i = 0
  nmismatch = 10
  while nmismatch > 1 :
    nmismatch = 0
    for aa,bb in zip(a[i:],b[:len(a[i:])]) :
      if aa != bb : nmismatch += 1
    i += 1
#  while a[i:] != b[:len(a[i:])]:
#    i += 1
  return i-1
#----------------------------------------------------------- 

#----------------------------------------------------------- 
def server_ready(f) :
  logger = logging.getLogger('root')
  while not os.path.exists(f) :
    logger.info("Waiting for log file %s"%f)
    time.sleep(10)
  ready = False
  f = open(f,'r')
  flines = f.readlines()
  for line in flines : 
    if line.find('Server ready for queries') > -1 : ready = True
  return ready
#----------------------------------------------------------- 
# End of utility functions
################################################################

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Class params
# Tracks the analysis-level parameters.
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
class params :
  def __init__(self,config_d) :
    self.opts = config_d
    self.gene_annotations = None
    self.targets = {}
    self.paths = {}
    self.logger = logging.getLogger('root')
    self.repeat_mask = None
    self.set_params()

  def set_targets(self,gene_list) :
    region_list = None
    if gene_list :
      region_list = []
      hl = open(gene_list,'r')
      for line in hl.readlines() :
        line = line.strip()
        region_list.append(line.upper())

    self.logger.info('Parsing target list')
    # Gene file
    # TODO: Check to make sure there aren't duplicate genes.
    targets_f = open(self.opts['targets_bed_file'],'rU')
    cur_region = ['',[]]
    for target in targets_f.readlines() :
      # Each target is formatted like a bed, chr bp1 bp2 name
      target = target.strip()
      targetsplit = target.split()
      chrm, bp1, bp2, name = targetsplit[0:4]
      if region_list :
        if name.upper() not in region_list : continue

      feature = None
      if len(targetsplit) > 4 : feature = targetsplit[4]

      if name.upper() not in self.targets : self.targets[name.upper()] = []
      self.targets[name.upper()].append((chrm,int(bp1),int(bp2),name,feature))
    self.logger.info('%d targets'%len(self.targets))

  def set_params(self) :
    self.logger.info('Setting up parameters')

    if 'indel_size' not in self.opts : self.opts['indel_size'] = 0
    var_filter = self.opts['var_filter']
    if var_filter == 'all' :
      self.opts['var_filter'] = ['indel','rearrangement','trl']
    else :
      self.opts['var_filter'] = var_filter.split(",")
      if 'indel' in  self.opts['var_filter'] or 'rearrangement' in self.opts['var_filter'] or 'trl' in self.opts['var_filter'] :
        self.logger.info('Variant filters %s set, only reporting these variants'%','.join(self.opts['var_filter']))
      else :
        self.logger.debug('Variant filter options %s are not valid, using all'%','.join(self.opts['var_filter']))
        self.opts['var_filter'] = ['indel','rearrangement','trl']

    # Log all parameters passed in, warn for poor paths
    for param in self.opts :
      value = self.opts[param]
      self.logger.info('%s = %s'%(param,value))

    self.set_targets(self.opts['gene_list'])
    # Set gene annotations
    self.gene_annotations = anno()
    self.gene_annotations.add_genes(self.opts['gene_annotation_file'])

    if 'other_regions_file' in self.opts :
      self.gene_annotations.add_regions(self.opts['other_regions_file'])

    # Setup analysis directories
    self.paths['analysis'] = os.path.abspath(os.path.normpath(self.opts['analysis_dir']))
    self.paths['output'] = os.path.join(self.paths['analysis'],'output')
    if 'targets_dir' in self.opts :
      self.paths['targets'] = os.path.abspath(os.path.normpath(self.opts['targets_dir']))
    else :
      self.paths['targets'] = os.path.join(self.paths['analysis'],'targets')
    self.paths['ref_data'] = os.path.abspath(os.path.normpath(self.opts['reference_data_dir']))

    for path in self.paths :
      self.logger.info('Creating %s directory (%s)'%(path,self.paths[path]))
      if not os.path.exists(self.paths[path]) : os.makedirs(self.paths[path])

    # Set repeats if specified
    if not self.opts['keep_repeat_regions'] and 'repeat_mask_file' in self.opts:
      self.logger.info('Storing all repeats by chrom from file %s'%self.opts['repeat_mask_file'])
      self.repeat_mask = setup_rmask_all(self.opts['repeat_mask_file'])

  def start_blat_server(self) :
    blat_bin_dir = os.path.split(self.opts['blat'])[0]
    self.opts['gfserver'] = os.path.join(blat_bin_dir,'gfServer')
    self.opts['gfclient'] = os.path.join(blat_bin_dir,'gfClient')
    fatotwobit = os.path.join(blat_bin_dir,'faToTwoBit')
    self.opts['reference_fasta_dir'] = os.path.split(self.opts['reference_fasta'])[0]
    ref_fasta_name = os.path.basename(self.opts['reference_fasta']).split(".fa")[0]
    # Check if 2bit is there.
    self.opts['blat_2bit'] = os.path.join(self.opts['reference_fasta_dir'],ref_fasta_name+".2bit")
    if not os.path.exists(self.opts['blat_2bit']) :
      self.logger.info('Creating 2bit from %s reference fasta'%ref_fasta_name+".fa")
      # Create 2bit requires faToTwoBit
      curdir = os.getcwd()
      os.chdir(self.opts['reference_fasta_dir'])
      cmd = '%s %s %s'%(fatotwobit,ref_fasta_name+".fa",ref_fasta_name+".2bit")
      p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
      output, errors = p.communicate()
      os.chdir(curdir)

    curdir = os.getcwd()
    os.chdir(self.opts['reference_fasta_dir'])
    # Start gfServer, change dir to 2bit file, gfServer start localhost 8000 .2bit 
    self.opts['blat_port'] = random.randint(8000,9500)
    self.opts['gfserver_log'] = os.path.join(self.paths['output'],'gfserver_%d.log'%self.opts['blat_port'])
    cmd = '%s -canStop -log=%s -stepSize=5 start localhost %d %s &'%(self.opts['gfserver'], self.opts['gfserver_log'], self.opts['blat_port'],ref_fasta_name+".2bit")
    self.logger.info("Starting gfServer %s"%cmd)
    p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
    start_time = time.time()
    while not server_ready(self.opts['gfserver_log']) :
      new_time = time.time()
      wait_time = new_time - start_time
      if wait_time > 1000 :
        self.logger.info('gfServer wait time exceeded ~15 minutes, exiting')
        sys.exit()
      self.logger.info('Waiting for blat gfServer to load reference seq')
      time.sleep(60)
    self.logger.info('Server ready!')
    os.chdir(curdir)

  def get_kmer_size(self) : 
    return int(self.opts['kmer_size'])

  def get_min_segment_length(self, type) :
    return int(self.opts[type + '_minseg_len'])

  def get_sr_thresh(self, type) :
    if type == 'min' :
      return min(self.get_sr_thresh('trl'), self.get_sr_thresh('rearrangement'), self.get_sr_thresh('indel'))
    else : 
      if type == 'trl' : 
        return int(self.opts['trl_sr_thresh'])
      elif type == 'rearrangement' :
        return int(self.opts['rearr_sr_thresh'])
      elif type == 'indel' :
        return int(self.opts['indel_sr_thresh'])
# End of params class
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Class fastq_read
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
class fq_read :
  def __init__(self, header, seq, qual, indel_only) :
    self.id = header
    self.seq = str(seq)
    self.qual = str(qual)
    self.used = False
    self.dup = False
    self.indel_only = indel_only

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Class FastqFile
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
class FastqFile(object) :
  def __init__(self,f) :
    if isinstance(f,str) :
      f = open(f)
      self._f = f
  def __iter__(self) :
    return self

  def next(self) :
    header, seq, qual_header, qual = [self._f.next() for _ in range(4)]
    header = header.strip()
    inst,lane,tile,x,y_end = header.split(':')
    seq = seq.strip()
    qual = qual.strip()
    bc = None
    y = y_end
    if y.find('/') > -1 :
      y, end = y.split('/')
    if y.find('#') > -1 :
      y, bc = y.split('#')
    header_dict = {'inst':inst,
                   'lane':int(lane),
                   'tile':int(tile),
                   'x':int(x),
                   'y':int(y),
                   'end':end,
                   'bc': bc}
    return (header,seq,qual)
# End FastqFile class
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@a
# Class anno
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
class anno :
  def __init__(self) :
    self.genes = {}
    self.logger = logging.getLogger('root')

  def add_genes(self, gene_fn) :
    self.logger.info('Adding gene annotations from %s'%gene_fn)
    gene_f = open(gene_fn,'r')
    gene_flines = gene_f.readlines()
    for line in gene_flines[1:] :
      line = line.strip()
      linesplit = line.split()
      chrom = linesplit[2]
      start = int(linesplit[4])
      end = int(linesplit[5])
      geneid = linesplit[12]
      if geneid in self.genes :
        if start <= self.genes[geneid][1] and end >= self.genes[geneid][2] : self.genes[geneid] = [chrom,start,end]
      else : self.genes[geneid] = [chrom,start,end]
    gene_f.close()

  def add_regions(self, regions_bed_fn) :
    region_f = open(regions_bed_fn, 'rU')
    region_lines = region_f.readlines()
    for line in region_lines :
      line = line.strip()
      chrom, start, end, name = line.split()
      if name not in self.genes : self.genes[name] = [chrom,int(start),int(end)]
    self.logger.info('Adding in %d other target regions'%len(region_lines))
  
  def set_gene(self, chrom, pos) :
    ann_genes = []
    if chrom.find('chr') == -1 : chrom = 'chr'+str(chrom) 
    for g in self.genes :
      gs = self.genes[g][1]
      ge = self.genes[g][2]
      if chrom == self.genes[g][0] :
        if len(pos) == 1 :
          if int(pos[0]) >= gs and int(pos[0]) <= ge :
            ann_genes.append(g)
            break
        else :
          # Find genes between pos1 and pos2
          if (int(pos[0]) >= gs and int(pos[0]) <= ge) or (int(pos[1]) >= gs and int(pos[1]) <= ge) :
            ann_genes.append(g)
    if len(ann_genes) == 0 : ann_genes = ['intergenic']
    return ",".join(ann_genes)

# End of anno class
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

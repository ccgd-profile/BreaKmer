#! /usr/bin/local/python

import sys
import re
import logging
import os
import olc
from collections import OrderedDict
from utils import *
#----------------------------------------------------------- 
def setup_contigs(mer, fq_recs, kmer_len, akmers, buff) :
  logger = logging.getLogger('root')
  ct = None
  mer_reads = find_reads(mer, fq_recs.items(), set())
  buff.add_used_mer(mer)
  for read_vals in mer_reads :
    read, mer_pos, bool, rlen, nreads = read_vals
    buff.add_used_read(read.id)
    if not ct : 
      ct = contig(mer, read, mer_pos, nreads, kmer_len)
      buff.add_contig(read, ct)
    else :
      ct.check_read(mer, akmers.get_count(mer), read, mer_pos, nreads, akmers.smers_set, 'setup')
  if ct : 
#    print 'Setup done'
    ct.finalize(fq_recs, akmers, buff, 'setup')
#----------------------------------------------------------- 

#-----------------------------------------------------------
def init_assembly(mers, fq_recs, kmer_len, rc_thresh, read_len) :
  logger = logging.getLogger('root')
  kmer_clusters = []
  if len(mers) == 0 : 
    return kmer_clusters

  ks = kmers()
  for mer in mers :
#    if mers[mer] > 1 :  
    ks.add_kmer(mer,mers[mer])

  buff = buffer()
  akmers = ks.get_all_kmer_values()
  while akmers.has_mers() :
    akmers.update_smer_set()
    mer, count = akmers.mers.items()[0]
    if count < 2 : 
      continue
    logger.info('Initiating kmer %s, found in %d reads'%(mer,count))
    setup_contigs(mer, fq_recs, kmer_len, akmers, buff)
    while len(buff.contigs) > 0 :
      ct = buff.get_contig()
      ct.grow(fq_recs, ks, akmers, kmer_len, buff)
      if ct.get_total_read_support() < int(rc_thresh) or ct.get_contig_len() <= read_len :
#        print 'Failed contig'
        logger.info('Contig did not meet the read count threshold %d, with %d or contig length (%d) < read_len (%d)'%(rc_thresh, len(ct.reads), len(ct.aseq.seq), read_len))
      else : 
#        print 'Contig done', ct.aseq.seq, len(ct.reads), len(ct.aseq.seq)
        logger.info('Adding contig to buffer')
        kmer_clusters.append(ct)
    buff.remove_kmers(akmers)
    buff.remove_reads(fq_recs)
#    print len(akmers.mers), 'sample kmers left'
  return kmer_clusters
#-----------------------------------------------------------

#----------------------------------------------------------- 
def same_reads(seq1, seq2) :
  same = False
  aln = olc.nw(seq1, seq2)
  if aln[3] == 0 and aln[5] == 0 and aln[6] > 0.95*(len(seq1)) :
    same = True
  return same
#----------------------------------------------------------- 

#------------------------------------------------------- 
# Check if seq2 is a subseq of seq1 
def subseq(seq1, seq2) : 
  aln = olc.nw(seq2, seq1)
  seq2_sub = (False,None)
  if aln[2] == len(seq2) and aln[3] == 0 and aln[6] >= (0.85*(len(seq2))) : 
    if len(seq2) < len(seq1) :
      seq2_sub = (True,None)
    else :
      seq2_sub = (True,aln[6]) 
  else :
    seq2_sub = (False,aln[6])
#  print seq2_sub, aln
  return seq2_sub
#------------------------------------------------------- 
 
#------------------------------------------------------- 
def sim_seqs(seq1, b_read) : 
  sim = False
  if not b_read.redundant : 
    seq2 = b_read.read.seq
    if same_reads(seq1, seq2) or subseq(seq2, seq1) : 
      sim = True 
  return sim 
#------------------------------------------------------- 
 
#------------------------------------------------------- 
def read_search(mer, read_values) : 
  r = (None, None, None) 
  seq, reads = read_values
  x = re.search(mer, seq) 
  if x : r = (reads[0], x.start(), True, len(reads[0].seq), len(reads)) 
  return r 
#------------------------------------------------------- 
 
#------------------------------------------------------- 
def find_reads(mer, read_items, used_reads, order='for' ): 
  mer_reads = [] 
  mr = filter(lambda x: x[2], map(read_search, [mer]*len(read_items), read_items))
  ids = map(lambda x: x[0].id, mr)
  filt_ids = set(ids) - set(used_reads)
  matched_reads = filter(lambda x: (x[0].id in filt_ids), mr)
#  matched_reads = filter(lambda y: y[2] and (not y[0].id in used_reads), map(read_search, [mer]*len(reads), reads)) 
  if order == 'rev' : 
    mer_reads = sorted(matched_reads, key=lambda z: (-z[1],-z[3])) 
  else : 
    mer_reads = sorted(matched_reads, key=lambda z: (z[1],-z[3])) 
  return mer_reads 
#------------------------------------------------------- 
 
#------------------------------------------------------- 
def get_read_kmers_ordered(seq, l, skmers, order='for') : 
  kmers = [] 
#  i = 0 
  m = len(seq)/2 
  kmers = map(lambda x: (seq[x:x+l],x,int(x<m), abs(x-m), order), range(0,(len(seq)-l)))
  ks = set(map(lambda x: x[0], kmers))
#  while (i+l) <= len(seq) :  
#    k = seq[i:i+l] 
#    kmers.append((k, i, int(i<m), abs(i-m), order)) 
#    i += 1 
#  kmers = filter(lambda x: x[0] in set(skmers), kmers) 
  ss = ks & skmers
  kmers = filter(lambda x: x[0] in ss, kmers)
  if order == 'rev' :  
    kmers.reverse() 
  elif order == 'mid' : 
    kmers = sorted(kmers, key=lambda x: (x[2], x[3])) 
  return kmers 
#------------------------------------------------------- 
 
#------------------------------------------------------- 
def get_read_kmers(seq, l, skmers) : 
#  kmers = [] 
#  i = 0 
#  while (i+l) <= len(seq) :  
#    kmers.append(seq[i:i+l]) 
#    i += 1 
  kmers = set(map(lambda x: seq[x:x+l], range(0,(len(seq)-l)))) 
  l = kmers & skmers 
  return l 
#-------------------------------------------------------
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Class assembly_counts
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
class assembly_counts :
  def __init__(self, read, nreads) :
    self.indel_only = [0]*len(read.seq)
#    self.all = [0]*len(read.seq)
    self.others = [0]*len(read.seq)
    self.set_counts(0, len(read.seq), nreads, read.indel_only)

  def get_counts(self, p1, p2, sv_type) :
    counts = []
    if sv_type == 'indel' or sv_type == 'rearr' :
      if p1 == p2 : counts = self.indel_only[p1] + self.others[p1]
      else : counts = map(lambda (x,y): x+y, zip(self.indel_only[p1:p2], self.others[p1:p2]))
    else :
      if p1 == p2 : counts = self.others[p1]
      else : 
        counts = self.others[p1:p2]
    return counts

  def get_total_reads(self) :
    return max(self.indel_only) + max(self.others) 

  def set_superseq(self, read, nreads, start, end) :
    # The counts vectors have to grow by the new seq length with nreads, then add back
    # the previous counts into the right positions indicated by start, end
    tmp_indel_only = [0]*len(read.seq)
    tmp_others = [0]*len(read.seq)
    if read.indel_only : 
      tmp_indel_only = [nreads]*len(read.seq)
    else :
      tmp_others = [nreads]*len(read.seq)
    tmp_indel_only[start:end] = map(lambda (x,y): x+y, zip(tmp_indel_only[start:end], self.indel_only))
    tmp_others[start:end] = map(lambda (x,y): x+y, zip(tmp_others[start:end], self.others))
    self.indel_only = tmp_indel_only
    self.others = tmp_others

  def set_counts(self, start, end, nreads, indel_only) :
    if indel_only :
      self.indel_only[start:end] = map(lambda x: x+nreads, self.indel_only[start:end])
    else :
      self.others[start:end] = map(lambda x: x+nreads, self.others[start:end])

  def extend_counts(self, l, nreads, indel_only, direction) :
    fill_counts = [0]*l
    ecounts = [nreads]*l
    if indel_only :
      if direction == 'post' : 
        self.indel_only.extend(ecounts)
        self.others.extend(fill_counts)
      else :
        ecounts.extend(self.indel_only)
        self.indel_only = ecounts
        fill_counts.extend(self.others)
        self.others = fill_counts
    else :
      if direction == 'post' :
        self.indel_only.extend(fill_counts)
        self.others.extend(ecounts)
      else :
        ecounts.extend(self.others)
        self.others = ecounts
        fill_counts.extend(self.indel_only)
        self.indel_only = fill_counts

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Class assembly_seq
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
class assembly_seq :
  def __init__(self, read, nreads) :
    self.seq = read.seq
#    print 'Started assembly seq', read.id, read.seq
    self.counts = assembly_counts(read, nreads)
  #*******************************************
  def set_superseq(self, read, nreads, start, end) :
#    print 'Set superseq', read.id, read.seq
    self.seq = read.seq
#    print 'New assembly seq', self.seq
    self.counts.set_superseq(read, nreads, start, end)
  #*******************************************
  def add_subseq(self, start, end, nreads, indel_only) :
#    print 'Add subseq'
    self.counts.set_counts(start, end, nreads, indel_only)
#    self.counts[start:end] = map(lambda x: x+nreads, self.counts[start:end])
  #*******************************************
  def add_postseq(self, post_seq, start, end, nreads, indel_only) :
#    print 'Add postseq'
#    print 'Before', self.seq
#    print 'Postseq', post_seq
    self.seq += post_seq
#    print 'After', self.seq
    self.counts.set_counts(start, end, nreads, indel_only)
    self.counts.extend_counts(len(post_seq), nreads, indel_only, 'post')
#    self.counts[start:end] = map(lambda x: x+nreads, self.counts[start:end])
#    postcounts = [nreads]*len(post_seq)
#    self.counts.extend(postcounts)
  #*******************************************
  def add_preseq(self, pre_seq, start, end, nreads, indel_only) :
#    print 'Adding preseq'
#    print 'Before', self.seq
#    print 'Preseq', pre_seq
    self.seq = pre_seq + self.seq
#    print 'After', self.seq
    self.counts.set_counts(start, end, nreads, indel_only)
    self.counts.extend_counts(len(pre_seq), nreads, indel_only, 'pre')
#    self.counts[start:end] = map(lambda x: x+nreads, self.counts[start:end]) 
#    precounts = [nreads]*len(pre_seq)
#    precounts.extend(self.counts)
#    self.counts = precounts
  #*******************************************

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Class kmers
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
class kmers : 
  def __init__(self) :
    self.s = []
  #*******************************************
  def add_kmer(self,mer,count) :
    if len(set(mer)) > 1 :
      self.s.append((int(count), mer, kmer(mer,count)))
  #*******************************************
  def get_all_kmer_values(self) :
    mers_sorted = sorted(self.s, key=lambda x: (int(x[0]), x[1]), reverse=True)
    am = akmers()
    for mer in mers_sorted :
      am.add_mer(mer[1], mer[0])
    return am
  #*******************************************

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Class kmer
# Tracks kmer variables
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
class kmer :
  def __init__(self,val,count) :
    self.value = val
    self.count = count
    self.done = False

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Class akmer 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
class akmers :
  def __init__(self) :
    self.mers = OrderedDict()
    self.smers_set = set()
  #*******************************************
  def get_mers(self) :
    return set(self.mers.keys())
  #*******************************************
  def update_smer_set(self) :
    self.smers_set = set(self.mers.keys())
  #*******************************************
  def add_mer(self, mer, count) :
    self.mers[mer] = count
   #*******************************************
  def remove_mer(self, mer) :
    del self.mers[mer]
  #*******************************************
  def has_mers(self) :
    if len(self.mers) > 0 and max(self.mers.values()) > 1 :
      return True
    else :
      return False
  #*******************************************
  def get_count(self, mer) :
    return self.mers[mer]
  #*******************************************

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Class buffer
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
class buffer :
  def __init__(self) :
    self.used_mers = set()
    self.used_reads = set()
    self.contigs = OrderedDict()
  #*******************************************
  def add_contig(self, read, contig) :
    if read.id not in self.contigs and not read.used :
      self.contigs[read.id] = contig
      read.used = True
  #*******************************************
  def remove_contig(self, read_id) :
    if read_id in self.contigs :
      del self.contigs[read_id]
  #*******************************************
  def get_contig(self) :
    read_id = self.contigs.keys()[0]
    ct = self.contigs[read_id]
    del self.contigs[read_id]
    return ct
  #*******************************************
  def add_used_read(self, read_id) :
    self.used_reads.add(read_id)
  #*******************************************
  def add_used_mer(self, mer) :
    self.used_mers.add(mer)
  #*******************************************
  def remove_kmers(self, skmers) :
    map(skmers.remove_mer, list(self.used_mers))
    self.used_mers = set()
  #*******************************************
  def remove_reads(self, fq_reads) :
    del_used = filter(lambda x: x in fq_reads, list(self.used_reads))
    map(fq_reads.__delitem__, del_used)
    self.used_reads = set()
  #*******************************************

# Class b_read
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
class b_read :
  def __init__(self, read, redundant, checked, aligned) :
    self.read = read
    self.redundant = redundant
    self.align_checked = checked
    self.aligned = aligned

# Class read_batch
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
class read_batch :
  def __init__(self, read, mer_pos) :
    self.delete = set()
    self.alt = []
    self.batch_reads = [b_read(read, False, True, True)]
    self.mer_pos_d = {mer_pos:[0]}
  #*******************************************
  def set_last_read_aligned(self) :
    self.batch_reads[-1].aligned = True
  #*******************************************
  def clean(self, fq_reads, buff, last_keep_read) :
    map(fq_reads.__delitem__, map(lambda x: x.seq, list(self.delete)))
    for read_id in filter(lambda x: x in buff.contigs, list(self.delete)) :
      if not buff.contigs[read_.id].setup :
        del buff.contigs[read.id]
    self.delete = set()
    self.alt = []
    self.batch_reads = [last_keep_read]
    self.mer_pos_d = {}
  #*******************************************
  def check_mer_read(self, pos, read) :
    check = True
    redund_read = False
    add_read = True
    add_to_pos_d = False 

    if add_read :
      if add_to_pos_d :
        if pos not in self.mer_pos_d : 
          self.mer_pos_d[pos] = []
        self.mer_pos_d[pos].append(len(self.batch_reads))
      self.batch_reads.append(b_read(read,False,check,False))
    return check 
  #*******************************************

# Class cluster_contig
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
class contig :
  def __init__(self, kmer_val, read, mer_pos, nreads, kmer_len) :
    self.reads = set()
    self.aseq = assembly_seq(read, nreads)
    self.kmer_locs = []
    self.kmers = [] 
    self.checked_kmers = [kmer_val]
    self.kmer_len = kmer_len
    self.buffer = set([read.id])
    self.rb = read_batch(read, mer_pos)
    self.setup = False
  #*******************************************
  def get_total_read_support(self) :
    return self.aseq.counts.get_total_reads()
  #*******************************************
  def get_contig_len(self) :
    return len(self.aseq.seq)
  #*******************************************
  def set_kmer_locs(self) :
    self.kmer_locs = [0]*len(self.aseq.seq)
    for kmer in self.kmers :
      kmer_pos = self.aseq.seq.find(kmer[0])
      self.kmer_locs[kmer_pos:(kmer_pos+self.kmer_len)] = map(lambda x: x+1, self.kmer_locs[kmer_pos:(kmer_pos+self.kmer_len)])  
  #*******************************************
  def get_kmer_locs(self) :
    return self.kmer_locs
  #*******************************************
  def get_contig_seq(self) :
    return self.aseq.seq
  #*******************************************
  def get_contig_counts(self) :
    return self.aseq.counts
  #*******************************************
  def check_align(self, kread, mer, nreads, skmers, type='setup') :
    match = False
    v1 = olc.nw(self.aseq.seq, kread.seq)
    v2 = olc.nw(kread.seq, self.aseq.seq)
#    print 'check_align()'
#    print 'Consensus seq', self.aseq.seq
#    print 'Read seq', kread.seq
#    print 'Resolve counts' , resolve_counts
#    print v1
#    print v2
    min_score = float(min(len(self.aseq.seq), len(kread.seq) )) / 4.0
    ident1 = round(float(v1[6]) / float(v1[2] - v1[3]), 2)
    ident2 = round(float(v2[6]) / float(v2[2] - v2[3]), 2)
#    print 'Min score', min_score, ', identity of overlap segments', ident1, ident2

    if (v1[6] < min_score or ident1 < 0.90) and (v2[6] < min_score or ident2 < 0.90) :
      return False
    if v1[6] == v2[6] and v1[3] == 0 and v1[5] == 0 and len(self.aseq.seq) == len(kread.seq) :
#      print 'Consensus and read sequence are the same'
      return True
    if v1[6] == v2[6] :
      match = True
      if len(self.aseq.seq) < len(kread.seq) or (v1[2] == len(self.aseq.seq) and v1[3] == 0) :
        # consensus sequence is a subseq of kread
#        max_counts = max(self.aseq.counts)
#        self.aseq = assembly_seq(kread, nreads)
#        print 'Consensus is subseq of kread'
        self.aseq.set_superseq(kread, nreads, v1[5], v1[4])
#        self.aseq.add_subseq(v1[5], v1[4], max_counts, kread.indel_only) 
        if type == 'grow' : 
          self.set_kmers(skmers)
      elif len(kread.seq) < len(self.aseq.seq) or (v2[2] == len(kread.seq) and v2[3] == 0):
#        print 'Consensus contains read seq'
        self.aseq.add_subseq(v2[5], v2[4], nreads, kread.indel_only)
      else :
        match = False
        indx11 = v1[0].replace('-','').find(mer)
        indx12 = v1[1].replace('-','').find(mer)
        indx21 = v2[0].replace('-','').find(mer)
        indx22 = v2[1].replace('-','').find(mer)
        if indx11 > -1 and indx12 > -1 :
          if (indx21 == -1 and indx22 == -1) or (abs(indx21 - indx22) > abs(indx11 - indx12)) :
            match = True
            self.contig_overlap_read(v1, kread, nreads, skmers, type)
        elif indx21 > -1 and indx22 > -1 :
          if (indx11 == -1 and indx12 == -1) or (abs(indx21 - indx22) < abs(indx11 - indx12)) :
            match = True
            self.read_overlap_contig(v2, kread, nreads, skmers, type)  
    elif v1[6] > v2[6] :
      match = True
      self.contig_overlap_read(v1, kread, nreads, skmers, type)
    else :
      # Read hangs left
      match = True
      self.read_overlap_contig(v2, kread, nreads, skmers, type)
    return match
  
  def contig_overlap_read(self, aln, kread, nreads, skmers, type) :
    # Consensus hangs left (consensus end overlaps read beginning) 
    if aln[2] == len(self.aseq.seq) and aln[3] == 0 : 
      # consensus sequence is a subseq of kread
#      max_counts = max(self.aseq.counts)
#      self.aseq = assembly_seq(kread, nreads)
#      print 'Consensus is subseq of kread'
      self.aseq.set_superseq(kread, nreads, aln[5], aln[4])
#      self.aseq.add_subseq(aln[5], aln[4], max_counts, kread.indel_only)
      if type == 'grow' : self.set_kmers(skmers)
    else :
#      print 'Pre', self.aseq.seq
#      print 'Overlap', v1[0], v1[1]
#      print 'Post', kread.seq[v1[4]:]
      post_seq = kread.seq[aln[4]:]
      nseq = self.aseq.seq[(len(self.aseq.seq)-(self.kmer_len-1)):] + post_seq
      self.aseq.add_postseq(post_seq, aln[3], aln[2], nreads, kread.indel_only)
      if type == 'grow' : 
#        print 'Added postseq', nseq
        nkmers = get_read_kmers_ordered(nseq, self.kmer_len, skmers, 'for')
#        print 'New kmers', nkmers
        self.kmers.extend(nkmers)
#      print 'New consensus, postseq', self.aseq.seq
  
  def read_overlap_contig(self, aln, kread, nreads, skmers, type) :
    if aln[2] == len(kread.seq) and aln[3] == 0 :
#      print 'Consensus contains read seq'
      self.aseq.add_subseq(aln[5], aln[4], nreads, kread.indel_only)
    else :
#      print 'Pre', kread.seq[0:v2[3]]
#      print 'Overlap', v2[0], v2[1]
#      print 'Post', self.aseq.seq[v2[4]:]
      pre_seq = kread.seq[0:aln[3]]
      nseq = pre_seq + self.aseq.seq[0:(self.kmer_len-1)]
      self.aseq.add_preseq(kread.seq[0:aln[3]], aln[5], aln[4], nreads, kread.indel_only)
      if type == 'grow' :
#        print 'Added preseq', nseq
        nkmers = get_read_kmers_ordered(nseq, self.kmer_len, skmers, 'rev')
#        print 'New kmers', nkmers
        self.kmers.extend(nkmers)
#      print 'New consensus, preseq', self.aseq.seq

  def set_kmers(self, skmers) :
    self.setup = True
    self.kmers = get_read_kmers_ordered(str(self.aseq.seq), self.kmer_len, skmers, 'mid')
  
  def check_read(self, mer, mer_count, read, mer_pos, nreads, skmers, type='setup') :
    hit = ''
    self.buffer.add(read.id)
    check_read_align = self.rb.check_mer_read(mer_pos, read)
    if check_read_align :
      match = self.check_align(read, mer, nreads, skmers, type)
      if match :
        hit = 'remove'
        read.used = True
        self.rb.set_last_read_aligned()
      elif mer_count > 2 and not read.used :
        self.rb.alt.append((read, nreads))
      else :
        self.rb.delete.add(read)
    return hit

  def check_alt_reads(self, akmers, buff) :
    new_contigs = []
    mer_set = set()
    for read in self.rb.alt :
      alt_kmers = get_read_kmers(read[0].seq, self.kmer_len, akmers.smers_set)
      x = set(alt_kmers) - set(self.kmers) - buff.used_mers - mer_set
      if len(x) > 0 :
        for mer in list(x) :
          read_count = akmers.get_count(mer)
          if read_count > 1 :
            mer_pos = read[0].seq.find(mer)
            new_contigs.append((read[0], contig(mer, read[0], mer_pos, read[1], self.kmer_len)))
            mer_set = mer_set | x
            break
    return new_contigs

  def finalize(self, fq_reads, akmers, buff, source='setup') :
    # Set kmers
    if source == 'setup' :
      self.set_kmers(akmers.smers_set)

    # Get alternate read kmers and see if any are different from contig kmers.
    new_contigs = self.check_alt_reads(akmers, buff)
    for nc in new_contigs :
      buff.add_contig(nc[0], nc[1])

    rm_reads = map(lambda y: y.read, filter(lambda x: x.redundant, self.rb.batch_reads))
    keep_reads = filter(lambda x: x.aligned and not x.redundant, self.rb.batch_reads)
    add_reads = map(lambda y: y.read, keep_reads)
    self.reads = self.reads | set(add_reads)
    self.reads = self.reads - set(rm_reads)
    self.rb.clean(fq_reads, buff, keep_reads[-1])

  def refresh_kmers(self) :
    return filter(lambda x: x[0] not in set(self.checked_kmers), self.kmers)

  def get_mer_reads(self, kmer_values, read_items) :
    mer, pos, less_than_half, dist_half, order = kmer_values
    read_order = 'for'
    if order == 'mid' :
      if less_than_half == 0 :
        read_order = 'rev'
    elif order == 'for' :
      read_order = 'rev'

    reads = find_reads(mer, read_items, self.buffer, read_order)
    return reads

  def grow(self, fq_recs, ks, akmers, kmer_len, buff) :
    logger = logging.getLogger('root')
    if not self.setup :
      self.set_kmers(akmers.smers_set)

#    print 'Growing consensus seq', self.aseq.seq
#    print 'Contig kmers', len(self.kmers)
#    print 'Contig checked kmers', len(self.checked_kmers)
#    print 'Contig reads', len(self.reads) #[x.id for x in list(self.reads)]

    nkmers = self.refresh_kmers()
    while len(nkmers) > 0 :
      kiter = 0
      for kmer_v in nkmers :
        mer, pos, less_than_half, dist_half, order = kmer_v
        reads = self.get_mer_reads(kmer_v, fq_recs.items())
        buff.add_used_mer(mer)
#        print 'Found reads from kmer', mer, [x[0].id + " " + str(x[0].used) for x in reads], len(reads)
        for read_v in reads :
          read, mer_pos, bool, rlen, nreads = read_v
          buff.add_used_read(read.id)
          hit = self.check_read(mer, akmers.get_count(mer), read, mer_pos, nreads, akmers.smers_set, 'grow')
          if hit == 'remove' :
            buff.remove_contig(read.id)
        self.finalize(fq_recs, akmers, buff, 'grow')
        self.checked_kmers.append(mer)
        kiter += 1
      nkmers = self.refresh_kmers()
      logger.debug("%d kmers left to check"%len(nkmers))
    self.set_kmer_locs() 
#    print 'Done', self.aseq.seq, len(self.reads)
#    print 'Kmers', self.kmer_locs
    logger.info('Contig done with contig seq %s. Supported by %d read(s).'%(self.aseq.seq, len(self.reads)))
    logger.info('Read IDs: %s'%(",".join([x.id for x in list(self.reads)])))

#      print 'Contig kmers', len(self.kmers), self.kmers
#      print 'Checked kmers', len(self.checked_kmers)
#      print len(akmers.get_mers()), 'Skmers remain'
#      print len(nkmers), 'kmers to check'
#    print '!'*20
#    print 'Contig done', self.aseq.seq 
#    print 'Indel only counts', self.aseq.counts.indel_only
#    print 'Other counts', self.aseq.counts.others
#    print 'Reads', [x.id for x in self.reads]
#    print '!'*20
#    print '\n'*2

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if __name__ == '__main__' :
  args = sys.argv[1:]
  in_fn, fq_in, rc_thresh, kmer_len = args
  rc_thresh = int(rc_thresh)
  kmer_len=int(kmer_len)
  fa_f = open(in_fn,"rU")
  fq_recs = {}
  read_len = 0
  for header,seq,qual in FastqFile(fq_in) : 
    fr = fq_read(header, seq, qual)
    read_len = max(read_len, len(fr.seq)) 
    seq = fr.seq
    if seq not in fq_recs :
      fq_recs[seq] = []
    fq_recs[seq].append(fr)

  iter = 0
  mers = {}
  for line in fa_f.readlines() :
    line = line.strip()
    linesplit = line.split()
    mer = linesplit[0]
    read_count = int(linesplit[1])
    mers[mer] = read_count
    iter += 1
    if (iter % 200000) == 0 : print iter

  contigs = init_assembly(mers, fq_recs, kmer_len, rc_thresh, read_len)

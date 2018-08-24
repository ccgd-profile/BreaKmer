#! /usr/bin/local/python

import sys
import re
import logging
import os
import olc
from collections import OrderedDict
from utils import *
from storage import open_demuxer

#----------------------------------------------------------- 
def setup_contigs(mer, fq_recs, kmer_len, akmers, buff) :
  logger = logging.getLogger('root')
  ct = None
  mer_reads = find_reads(mer, fq_recs.values(), set())
  buff.add_used_mer(mer)
  for read_vals in mer_reads :
    read, mer_pos, bool, rlen = read_vals
    buff.add_used_read(read.id)
    if not ct : 
      ct = contig(mer, read, mer_pos, kmer_len)
      buff.add_contig(read, ct)
    else :
      ct.check_read(mer, akmers.get_count(mer), read, mer_pos, akmers.smers_set, 'setup')
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
    if mers[mer] > 1 :  
      ks.add_kmer(mer,mers[mer])

  buff = buffer()
  akmers = ks.get_all_kmer_values()
  while akmers.has_mers() :
    akmers.update_smer_set()
    mer, count = akmers.mers.items()[0]
    logger.info('Initiating kmer %s, found in %d reads'%(mer,count))
    setup_contigs(mer, fq_recs, kmer_len, akmers, buff)
    while len(buff.contigs) > 0 :
      ct = buff.get_contig()
      ct.grow(fq_recs, ks, akmers, kmer_len, buff)
      if len(ct.reads) < int(rc_thresh) or len(ct.aseq.seq) < read_len :
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
  if aln[2] == len(seq2) and aln[3] == 0 and aln[6] >= (0.90*(len(seq2))) : 
    if len(seq2) < len(seq1) :
      seq2_sub = (True,None)
    else :
      seq2_sub = (True,aln[6]) 
  else :
    seq2_sub = (False,aln[6])
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
def read_search(mer, read) : 
  r = (None, None, None) 
  x = re.search(mer, read.seq) 
  if x : r = (read, x.start(), True, len(read.seq)) 
  return r 
#------------------------------------------------------- 
 
#------------------------------------------------------- 
def find_reads(mer, reads, used_reads, order='for' ): 
  mer_reads = [] 
  mr = filter(lambda x: x[2], map(read_search, [mer]*len(reads), reads))
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
# Class assembly_seq
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
class assembly_seq :
  def __init__(self, seq) :
    self.seq = seq
    self.counts = [1]*len(seq)
  #*******************************************
  def add_subseq(self, start, end) :
    self.counts[start:end] = map(lambda x: x+1, self.counts[start:end])
  #*******************************************
  def add_postseq(self, post_seq, start, end) :
    self.counts[start:end] = map(lambda x: x+1, self.counts[start:end])
    postcounts = [1]*len(post_seq)
    self.counts.extend(postcounts)
    self.seq += post_seq
  #*******************************************
  def add_preseq(self, pre_seq, start, end) :
    self.counts[start:end] = map(lambda x: x+1, self.counts[start:end]) 
    precounts = [1]*len(pre_seq)
    precounts.extend(self.counts)
    self.counts = precounts
    self.seq = pre_seq + self.seq
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
    mers_sorted = sorted(self.s, key=lambda x: int(x[0]), reverse=True)
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
  def clean(self, fq_reads, buff) :
    map(fq_reads.__delitem__, list(self.delete))
    for read_id in filter(lambda x: x in buff.contigs, list(self.delete)) :
      if not buff.contigs[read_id].setup :
        del buff.contigs[read_id]
    self.delete = set()
    self.alt = []
    self.batch_reads = []
    self.mer_pos_d = {}
  #*******************************************
  def check_mer_read(self, pos, read) :
    check = True
    redund_read = False
    add_read = True
    add_to_pos_d = False 

    if not pos in self.mer_pos_d : 
      add_to_pos_d = True
    else :
      x = filter(lambda y: y, map(sim_seqs, [read.seq]*len(self.mer_pos_d[pos]), [self.batch_reads[x] for x in self.mer_pos_d[pos]]))
      if len(x) > 0 :
        self.delete.add(read.id)
        add_read = False
        check = False
      else : 
        add_to_pos_d = True

    if add_read :
      if len(self.batch_reads) > 0 :
        ss1 = subseq(self.batch_reads[-1].read.seq, read.seq)
        if ss1[0] and not ss1[1] : 
          self.delete.add(read.id)
          add_read = False
          check = False
        else :
          ss2 = subseq(read.seq, self.batch_reads[-1].read.seq)  
          if ss2[0] and not ss2[1] : 
            self.delete.add(self.batch_reads[-1].read.id)
            self.batch_reads[-1].redundant = True
          elif (ss1[0] and ss1[1]) or (ss2[0] and ss2[1]) :
            if ss1[0] and ((ss1[1] > ss2[1]) or (ss1[1]==ss2[1])) : 
              # Read is subseq
              self.delete.add(read.id)
              add_read = False
              check = False            
            elif ss2[0] and ((ss2[1] > ss1[1]) or (ss1[1]==ss2[1])) :
              # Previous read is subseq
              self.delete.add(self.batch_reads[-1].read.id)
              self.batch_reads[-1].redundant = True

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
  def __init__(self, kmer_val, read, mer_pos, kmer_len) :
    self.reads = set()
    self.aseq = assembly_seq(read.seq)
    self.kmers = [] 
    self.checked_kmers = [kmer_val]
    self.kmer_len = kmer_len
    self.buffer = set([read.id])
    self.rb = read_batch(read, mer_pos)
    self.setup = False
  #*******************************************
  def get_contig_seq(self) :
    return self.aseq.seq
  #*******************************************
  def get_contig_counts(self) :
    return self.aseq.counts
  #*******************************************
  def check_align(self, kread, mer, skmers, type='setup') :
    match = False
    v1 = olc.nw(self.aseq.seq, kread.seq)
    v2 = olc.nw(kread.seq, self.aseq.seq)
#    print 'check_align()'
#    print 'Consensus seq', self.aseq.seq
#    print 'Read seq', kread.seq
#    print v1
#    print v2
    min_score = float(min(len(self.aseq.seq), len(kread.seq) )) / 4.0
    ident1 = float(v1[6]) / float(v1[2] - v1[3])
    ident2 = float(v2[6]) / float(v2[2] - v2[3])
#    print 'Min score', min_score, ', identity of overlap segments', ident1, ident2

    if (v1[6] < min_score or ident1 < 0.9) and (v2[6] < min_score or ident2 < 0.9) :
      return False
    if v1[6] == v2[6] and v1[3] == 0 and v1[5] == 0 :
#      print 'Consensus and read sequence are the same'
      return True
    if v1[6] == v2[6] :
      match = True
      if len(self.aseq.seq) < len(kread.seq) or (v1[2] == len(self.aseq.seq) and v1[3] == 0) :
        # consensus sequence is a subseq of kread
        self.aseq = assembly_seq(kread.seq)
#        print 'Consensus is subseq of kread'
        self.aseq.add_subseq(v1[5], v1[4]) 
        if type == 'grow' : 
          self.set_kmers(skmers)
      elif len(kread.seq) < len(self.aseq.seq) or (v2[2] == len(kread.seq) and v2[3] == 0):
#        print 'Consensus contains read seq'
        self.aseq.add_subseq(v2[5], v2[4])
      else :
        match = False
        indx11 = v1[0].replace('-','').find(mer)
        indx12 = v1[1].replace('-','').find(mer)
        indx21 = v2[0].replace('-','').find(mer)
        indx22 = v2[1].replace('-','').find(mer)
        if indx11 > -1 and indx12 > -1 :
          if (indx21 == -1 and indx22 == -1) or (abs(indx21 - indx22) > abs(indx11 - indx12)) :
            match = True
            self.contig_overlap_read(v1, kread, skmers, type)
        elif indx21 > -1 and indx22 > -1 :
          if (indx11 == -1 and indx12 == -1) or (abs(indx21 - indx22) < abs(indx11 - indx12)) :
            match = True
            self.read_overlap_contig(v2, kread, skmers, type)  
    elif v1[6] > v2[6] :
      match = True
      self.contig_overlap_read(v1, kread, skmers, type)
    else :
      # Read hangs left
      match = True
      self.read_overlap_contig(v2, kread, skmers, type)
    return match
  
  def contig_overlap_read(self, aln, kread, skmers, type) :
    # Consensus hangs left (consensus end overlaps read beginning) 
    if aln[2] == len(self.aseq.seq) and aln[3] == 0 : 
      # consensus sequence is a subseq of kread
      self.aseq = assembly_seq(kread.seq)
#      print 'Consensus is subseq of kread'
      self.aseq.add_subseq(aln[5], aln[4])
      if type == 'grow' : self.set_kmers(skmers)
    else :
#      print 'Pre', self.aseq.seq
#      print 'Overlap', v1[0], v1[1]
#      print 'Post', kread.seq[v1[4]:]
      post_seq = kread.seq[aln[4]:]
      nseq = self.aseq.seq[(len(self.aseq.seq)-(self.kmer_len-1)):] + post_seq
      self.aseq.add_postseq(post_seq, aln[3], aln[2])
      if type == 'grow' : 
#        print 'Added postseq', nseq
        nkmers = get_read_kmers_ordered(nseq, self.kmer_len, skmers, 'for')
#        print 'New kmers', nkmers
        self.kmers.extend(nkmers)
#      print 'New consensus, postseq', self.aseq.seq
  
  def read_overlap_contig(self, aln, kread, skmers, type) :
    if aln[2] == len(kread.seq) and aln[3] == 0 :
#      print 'Consensus contains read seq'
      self.aseq.add_subseq(aln[5], aln[4])
    else :
#      print 'Pre', kread.seq[0:v2[3]]
#      print 'Overlap', v2[0], v2[1]
#      print 'Post', self.aseq.seq[v2[4]:]
      pre_seq = kread.seq[0:aln[3]]
      nseq = pre_seq + self.aseq.seq[0:(self.kmer_len-1)]
      self.aseq.add_preseq(kread.seq[0:aln[3]], aln[5], aln[4])
      if type == 'grow' :
#        print 'Added preseq', nseq
        nkmers = get_read_kmers_ordered(nseq, self.kmer_len, skmers, 'rev')
#        print 'New kmers', nkmers
        self.kmers.extend(nkmers)
#      print 'New consensus, preseq', self.aseq.seq

  def set_kmers(self, skmers) :
    self.setup = True
    self.kmers = get_read_kmers_ordered(str(self.aseq.seq), self.kmer_len, skmers, 'mid')
  
  def check_read(self, mer, mer_count, read, mer_pos, skmers, type='setup') :
    hit = ''
    self.buffer.add(read.id)
    if self.rb.check_mer_read(mer_pos, read) :
      match = self.check_align(read, mer, skmers, type)
      if match :
        hit = 'remove'
        read.used = True
        self.rb.set_last_read_aligned()
      elif mer_count > 2 and not read.used :
        self.rb.alt.append(read)
      else :
        self.rb.delete.add(read.id)
    return hit

  def check_alt_reads(self, akmers, buff) :
    new_contigs = []
    mer_set = set()
    for read in self.rb.alt :
      alt_kmers = get_read_kmers(read.seq, self.kmer_len, akmers.smers_set)
      x = set(alt_kmers) - set(self.kmers) - buff.used_mers - mer_set
      if len(x) > 0 :
        for mer in list(x) :
          read_count = akmers.get_count(mer)
          if read_count > 1 :
            mer_pos = read.seq.find(mer)
            new_contigs.append((read, contig(mer, read, mer_pos, self.kmer_len)))
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

    add_reads = map(lambda y: y.read, filter(lambda x: x.aligned and not x.redundant, self.rb.batch_reads))
    self.reads = self.reads | set(add_reads)
    self.rb.clean(fq_reads, buff)

  def refresh_kmers(self) :
    return filter(lambda x: x[0] not in set(self.checked_kmers), self.kmers)

  def get_mer_reads(self, kmer_values, fq_recs) :
    mer, pos, less_than_half, dist_half, order = kmer_values
    read_order = 'for'
    if order == 'mid' :
      if less_than_half == 0 :
        read_order = 'rev'
    elif order == 'for' :
      read_order = 'rev'

    reads = find_reads(mer, fq_recs.values(), self.buffer, read_order)
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
        reads = self.get_mer_reads(kmer_v, fq_recs)
        buff.add_used_mer(mer)
#        print 'Found reads from kmer', mer, [x[0].id + " " + str(x[0].used) for x in reads], len(reads)
        for read_v in reads :
          read, mer_pos, bool, rlen = read_v
          buff.add_used_read(read.id)
          hit = self.check_read(mer, akmers.get_count(mer), read, mer_pos, akmers.smers_set, 'grow')
          if hit == 'remove' :
            buff.remove_contig(read.id)
        self.finalize(fq_recs, akmers, buff, 'grow')
        self.checked_kmers.append(mer)
        kiter += 1
      nkmers = self.refresh_kmers()
      logger.debug("%d kmers left to check"%len(nkmers))
    logger.info('Contig done with contig seq %s. Supported by %d read(s).'%(self.aseq.seq, len(self.reads)))
    logger.info('Read IDs: %s'%(",".join([x.id for x in list(self.reads)])))

#      print 'Contig kmers', len(self.kmers), self.kmers
#      print 'Checked kmers', len(self.checked_kmers)
#      print len(akmers.get_mers()), 'Skmers remain'
#      print len(nkmers), 'kmers to check'
#    print '!'*20
#    print 'Contig done', self.aseq.seq, len(self.reads), 'reads' # [x.id for x in self.reads], len(self.reads), 'reads'
#    print '!'*20
#    print '\n'*2

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if __name__ == '__main__' :
  args = sys.argv[1:]
  in_fn, fq_in, rc_thresh, kmer_len = args
  rc_thresh = int(rc_thresh)
  kmer_len=int(kmer_len)
  fa_f = open_demuxer(in_fn,"rU")
  fq_recs = {}
  read_len = 0
  for header,seq,qual in FastqFile(fq_in) : 
    fr = fq_read(header, seq, qual)
    read_len = max(read_len, len(fr.seq)) 
    fq_recs[fr.id] = fr

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

#! /usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
import glob
import logging
import time
from Bio import SeqIO
import subprocess
from pysam import *

__author__ = "Ryan Abo"
__copyright__ = "Copyright 2015, Ryan Abo"
__email__ = "ryanabo@gmail.com"
__license__ = "MIT"


def run_cutadapt(cutadapt, cutadapt_config_f, input_fn, output_fn, logging_src):
    """
    """
    cutadapt_parameters = stringify(cutadapt_config_f)
    cmd = '%s %s %s %s > %s' % (sys.executable, cutadapt, cutadapt_parameters, input_fn, output_fn)
    log(logging_src, 'debug', 'Cutadapt system command %s' % cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    output, errors = p.communicate()
    log(logging_src, 'debug', 'Cutadapt output %s' % output)
    log(logging_src, 'debug', 'Cutadapt errors %s' % errors)
    return output, errors


def log(name, level, msg):
        """Write log message to the appropriate level.
        """
        logger = logging.getLogger(name)
        if level == 'info':
                logger.info(msg)
        elif level == 'debug':
                logger.debug(msg)
        elif level == 'error':
                logger.error(msg)


def stringify(fn):
        """Turn file contents into a space delimited string"""

        str = []
        for line in open(fn, 'rU').readlines():
                line = line.strip()
                str.append(line)
        return ' '.join(str)


def create_ref_test_fa(target_fa_in, test_fa_out):
        """
        """
        if not os.path.isfile(get_marker_fn(test_fa_out)):
                fa_in = open(target_fa_in, "rU")
                fa_out = open(test_fa_out, "w")

                record = SeqIO.read(fa_in, "fasta")
                ref_target_seq = str(record.seq)
                end = min(len(ref_target_seq), 1500)
                start = max(0, len(ref_target_seq)-1500)
                fa_out.write(">"+record.id + "_start\n" + ref_target_seq[0:end] + "\n>" + record.id + "_end\n" + ref_target_seq[start:len(ref_target_seq)] + "\n")
                fa_out.close()

                cmd = 'touch %s' % get_marker_fn(test_fa_out)
                p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                output, errors = p.communicate()
                return True
        else:
                return False


def get_altref_genecoords(blat_path, altref_fa, query_fa_fn, chr, out_fn):
        altref_twobit = os.path.splitext(altref_fa)[0] + ".2bit"
        blat_db = altref_twobit + ":" + str(chr)
        cmd = "%s -noHead %s %s %s" % (blat_path, blat_db, query_fa_fn, out_fn)
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        output, errors = p.communicate()

        coords = [[0, 0], [0, 0], False]
        blat_res = open(out_fn, 'rU')
        hits = [False, False]
        for line in blat_res.readlines():
                line = line.strip()
                linesplit = line.split()
                res_id = linesplit[9]
                if res_id.find("start") > -1:
                        if coords[0][0] < int(linesplit[0]):
                                coords[0][0] = int(linesplit[0])
                                coords[0][1] = int(linesplit[15])
                                hits[0] = True
                elif res_id.find("end") > -1:
                        if coords[1][0] < int(linesplit[0]):
                                coords[1][0] = int(linesplit[0])
                                coords[1][1] = int(linesplit[16])
                                hits[1] = True
        coords[2] = hits[0] and hits[1]
        blat_res.close()
        return coords


def test_cutadapt(fq_fn, cutadapt_bin, cutadapt_config):
        fq_clean = os.path.basename(fq_fn).split('.')[0] + "_cleaned.fq"
        fq_clean_fn = os.path.join(os.path.dirname(fq_fn), fq_clean)
        cutadapt_params = stringify(cutadapt_config)
        cmd = '%s %s %s %s > %s' % (sys.executable, cutadapt_bin, cutadapt_params, fq_fn, fq_clean_fn)
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        output, errors = p.communicate()
        rc = p.returncode
        if rc != 0:
                return (None, rc)
        else:
                return (fq_clean_fn, rc)


def test_jellyfish(jfish_bin, fa_fn, analysis_dir):
    cmd = '%s --version' % jfish_bin
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    output, errors = p.communicate()
    jfish_version = int(output.split()[1].split('.')[0])

    kmer_size = 15
    count_fn = os.path.join(analysis_dir, "test_jellyfish_counts")
    cmd = '%s count -m %d -s %d -t %d -o %s %s' % (jfish_bin, kmer_size, 100000000, 8, count_fn, fa_fn)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    output, errors = p.communicate()
    if p.returncode != 0:
        return ("Jellyfish counts", p.returncode)

    if jfish_version < 2:
            count_fn += '_0'
    dump_fn = os.path.join(analysis_dir, "test_jellyfish_dump")
    cmd = '%s dump -c -o %s %s' % (jfish_bin, dump_fn, count_fn)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    output, errors = p.communicate()
    if p.returncode != 0:
            return ("Jellyfish dump", p.returncode)
    return ("Jellyfish", 0)


def write_test_fq(fq_fn):
    fq_f = open(fq_fn, 'w')
    fq_f.write("@H91H9ADXX140327:1:2102:19465:23489/2\nCACCCCCACTGAAAAAGATGAGTATGCCTGCCGTGTGAACCATGTGACTTTACAATCTGCATATTGGGATTGTCAGGGAATGTTCTTAAAGATC\n+\n69EEEFBAFBFABCCFFBEFFFDDEEHHDGH@FEFEFCAGGCDEEEBGEEBCGBCCGDFGCBBECFFEBDCDCEDEEEAABCCAEC@>>BB?@C\n@H91H9ADXX140327:2:2212:12198:89759/2\nTCTTGTACTACACTGAATTCACCCCCACTGAAAAAGATGAGTATGCCTGCCGTGTGAACCATGTGACTTTACAATCTGCATATTGGGATTGTCAGGGA\n+\nA@C>C;?AB@BBACDBCAABBDDCDDCDEFCDDDDEBBFCEABCGDBDEEF>@GBGCEDGEDGCGFECAACFEGDFFGFECB@DFGCBABFAECEB?=")
    fq_f.close()


def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
        if is_exe(exe_file):
            return exe_file
    return None


def calc_contig_complexity(seq, N=3, w=6):
    cmers = []
    for i in range(len(seq)):
        s = max(0, i-w)
        e = min(len(seq), i+w)
        cmer = count_nmers(seq[s:e], N)
        n = len(cmer)
        nmod = float(n) / float(e-s)
        cmers.append(round(nmod, 2))
    cmers_mean = sum(cmers) / len(cmers)
    return cmers_mean, cmers


def count_nmers(seq, N):
    nmers = {}
    for i in range(len(seq) - (N - 1)):
        mer = str(seq[i:i+N]).upper()
        if mer not in nmers:
            nmers[mer] = 0
        nmers[mer] += 1
    return nmers


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


def filter_by_feature(brkpts, query_region, keep_intron_vars):
    """
    """
    in_filter = False
    span_filter = False
    if not keep_intron_vars:
        in_vals, span_vals = check_intervals(brkpts, query_region)
        if in_vals[0]:
            if 'exon' not in in_vals[1]:
                in_filter = True
        else:
            in_filter = True
        if span_vals[0]:
            if 'exon' not in span_vals[1]:
                span_filter = True
        else:
            span_filter = True
    return in_filter, span_filter


def check_intervals(breakpts, query_region):
    """
    """
    in_values = [False, [], []]
    span_values = [False, [], []]
    for bp in breakpts:
        for interval in query_region[4]:
            if (int(bp) >= (interval[1] - 20)) and (int(bp) <= (interval[2] + 20)):
                in_values[1].append(interval[4])
                in_values[0] = True
                in_values[2].append(interval)
            if (interval[2] <= max(breakpts)) and (interval[1] >= min(breakpts)):
                span_values[0] = True
                span_values[1].append(interval[4])
                span_values[2].append(interval)
    return in_values, span_values


def setup_logger(logFnPath, name):
    """Creates the logger object and associated text file to use throughout
    the analysis.
    It first creates a log.txt file in the specified analyis directory as the
    FileHandler. The console handler is then formatted to report the time,
    name of the source, level of the message and the message.
    Args:
        log_fn_path: Absolute path for the directory that will contain the
        log file.
        name: The name of the package to initially setup the logger object.
    Returns:
        Nothing is returned.
    """
    outputPath = os.path.abspath(logFnPath)
    if not os.path.exists(outputPath):
        os.makedirs(outputPath)

    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    # FileHandler
    fileHandle = logging.FileHandler(os.path.join(outputPath, 'log.txt'), mode='w')
    fileHandle.setLevel(logging.DEBUG)
    # ConsoleHandler
    consoleHandle = logging.StreamHandler()
    consoleHandle.setLevel(logging.ERROR)

    formatter = logging.Formatter(fmt='%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    fileHandle.setFormatter(formatter)
    consoleHandle.setFormatter(formatter)

    logger.addHandler(fileHandle)
    logger.addHandler(consoleHandle)


def check_repeat_regions(coords, repeat_locs):
    """
    """
    start, end = coords
    seg_len = float(end-start)
    in_repeat = False
    rep_overlap = 0.0
    rep_coords = []
    filter_reps_edges = [False, False]
    for rloc in repeat_locs:
        rchr, rbp1, rbp2, rname = rloc
        if (rbp1 >= start and rbp1 <= end) or (rbp2 >= start and rbp2 <= end) or (rbp1 <= start and rbp2 >= end):
            in_repeat = True
            rep_overlap += float(min(rbp2, end)-max(rbp1, start))
            rep_coords.append((rbp1, rbp2))
            # Simple or low complexity seq repeat for filtering
            if rname.find(")n") > -1 or rname.find("_rich") > -1:
                if (rbp1 <= start and rbp2 >= start):
                    filter_reps_edges[0] = True
                elif (rbp1 <= end and rbp2 >= end):
                    filter_reps_edges[1] = True
    roverlap = round((float(min(rep_overlap, seg_len)) / float(seg_len)) * 100, 2)
    return in_repeat, roverlap, rep_coords, filter_reps_edges


def get_marker_fn(fn):
    return os.path.join(os.path.split(fn)[0], "."+os.path.basename(fn))


def run_jellyfish(fa_fn, jellyfish, kmer_size):
    logger = logging.getLogger('root')
    file_path = os.path.split(fa_fn)[0]
    file_base = os.path.basename(fa_fn)
    dump_fn = os.path.join(file_path, file_base + "_" + str(kmer_size) + "mers_dump")
    dump_marker_fn = get_marker_fn(dump_fn)
    if not os.path.isfile(dump_marker_fn):
        if not os.path.exists(fa_fn):
            logger.info('%s does not exist.'%fa_fn)
            dump_fn = None
            return dump_fn

        cmd = '%s --version'%jellyfish
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        output, errors = p.communicate()
        jfish_version = int(output.split()[1].split('.')[0])
        logger.info('Using jellyfish version %d'%jfish_version)

        count_fn = os.path.join(file_path, file_base + "_" + str(kmer_size) + "mers_counts")
        logger.info('Running %s on file %s to determine kmers'%(jellyfish,fa_fn)) 
        cmd = '%s count -m %d -s %d -t %d -o %s %s'%(jellyfish,kmer_size,100000000,8,count_fn,fa_fn)
        logger.info('Jellyfish counts system command %s'%cmd)
        p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
        output, errors = p.communicate()
        logger.info('Jellyfish count output %s'%output)
        logger.info('Jellyfish count errors %s'%errors)

        if jfish_version < 2: count_fn += '_0'
        cmd = '%s dump -c -o %s %s'%(jellyfish,dump_fn,count_fn)
        logger.info('Jellyfish dump system command %s'%cmd)
        p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
        output, errors = p.communicate()
        logger.info('Jellyfish dump output %s'%output)
        logger.info('Jellyfish dump errors %s'%errors)
        cmd = 'touch %s'%dump_marker_fn
        p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
        output, errors = p.communicate()  
        logger.info('Completed jellyfish dump %s, touching marker file %s'%(dump_fn,dump_marker_fn))
        count_fns = glob.glob(os.path.join(file_path, "*mers_counts*"))
        for cf in count_fns:
            os.remove(cf)
    else:
        logger.info('Jellfish already run and kmers already generated for target.')
    return dump_fn


def setup_ref_data(setup_params):
     genes = setup_params[0]
     rep_mask, ref_fa, altref_fa_fns, ref_path, jfish_path, blat_path, kmer_size = setup_params[1]

     logger = logging.getLogger('root')

     for gene in genes:
         chr, bp1, bp2, name, intvs = gene
         gene_ref_path = os.path.join(ref_path,name)
         if rep_mask: 
             logger.info('Extracting repeat mask regions for target gene %s.'%name)
             setup_rmask(gene, gene_ref_path, rep_mask)
         
         logger.info('Extracting refseq sequence for %s, %s:%d-%d'%(name, chr, bp1, bp2))
         directions = ['forward', 'reverse']
         for dir in directions:
             target_fa_fn = os.path.join(gene_ref_path, name + '_' + dir + '_refseq.fa')
             ref_fn = extract_refseq_fa(gene, gene_ref_path, ref_fa, dir, target_fa_fn)
             run_jellyfish(ref_fn, jfish_path, kmer_size)

         if altref_fa_fns: 
             if not create_ref_test_fa(os.path.join(gene_ref_path, name + '_forward_refseq.fa'), os.path.join(gene_ref_path, name + '_start_end_refseq.fa')):
                 return

             altref_fns = []
             alt_iter = 1
             altref_fas = altref_fa_fns.split(',')
             for altref in altref_fas:
                 for dir in directions:
                        fn = os.path.join(gene_ref_path, name + '_' + dir + '_altrefseq_' + str(alt_iter) + '.fa')
                        marker_fn = get_marker_fn(fn) 
                        if not os.path.isfile(marker_fn):
                            altref_fns.append((altref, fn, alt_iter))
                 alt_iter += 1
            
             if len(altref_fns) > 0:
                 altref_fas = altref_fa_fns.split(',')
                 alt_iter = 1
                 for i in range(len(altref_fns)):
                     alt_gene_coords = get_altref_genecoords(blat_path, altref_fns[i][0], os.path.join(gene_ref_path, name + '_start_end_refseq.fa'), chr, os.path.join(gene_ref_path, name + '_altref_blat_' + str(altref_fns[i][2]) + '.psl'))
                     if not alt_gene_coords[2]:
                         logger.info("No sequence for target gene %s in %s, no reference kmers extracted."%(name, altref_fns[i][0]))
                         alt_iter += 1
                         continue
                     gene = (chr, alt_gene_coords[0][1], alt_gene_coords[1][1], name, intvs)
                     target_fa_fn = altref_fns[i][1] #os.path.join(gene_ref_path, name + '_' + dir + '_altrefseq_' + str(alt_iter) + '.fa')
                     ref_fn = extract_refseq_fa(gene, gene_ref_path, altref_fns[i][0], dir, target_fa_fn)
                     run_jellyfish(ref_fn, jfish_path, kmer_size)
                 os.remove(os.path.join(gene_ref_path, name + '_start_end_refseq.fa'))


def get_fastq_reads(fn, sv_reads):
    read_len = 0
    filtered_fq_fn = fn.split(".fastq")[0] + "_filtered.fastq"
    filt_fq = open(filtered_fq_fn, 'w')
    fq_recs = {}
#  f = open(fn,'r')
#  fq_recs = list(SeqIO.parse(f,'fastq'))
    for header,seq,qual in FastqFile(fn): 
        qname_split = header.lstrip("@").split("_")
        indel_only = qname_split[-1]
        qname = "_".join(qname_split[0:len(qname_split)-1])
        if qname in sv_reads: 
            oseq, sc_seqs, clip_coords, indel_meta = sv_reads[qname]
            cleaned_seq = seq 
            old_seq = oseq.seq
            add = True
            if str(cleaned_seq) != str(old_seq) and sc_seqs:
                sc_clips = sc_seqs['clipped']
                idx = old_seq.find(cleaned_seq)
                trimmed_seq = ''
                if idx == 0: 
                    trimmed_seq = old_seq[len(cleaned_seq):len(old_seq)]
                else: trimmed_seq = old_seq[0:idx]    
                sc_lens = 0
                for sc_seq in sc_clips: 
                    sc_lens += len(sc_seq)
                    if trimmed_seq.find(sc_seq) > -1: 
                        add = False
                if len(cleaned_seq) == (len(old_seq) - sc_lens):
                    for sc_seq in sc_clips: 
                        if cleaned_seq.find(sc_seq) == -1:
                            # Don't add, just trimmed clipped portion.
                            add = False
#    else: print qname, 'not in sv reads'
        if add:
            filt_fq.write(header + "\n" + seq + "\n+\n" + qual + "\n") 
            fr = fq_read(header, seq, qual, indel_meta)
            read_len = max(read_len, len(fr.seq)) 
            seq = fr.seq
            if seq not in fq_recs: 
                fq_recs[seq] = []
            fq_recs[seq].append(fr)
    filt_fq.close()
    return filtered_fq_fn, fq_recs, read_len


def get_fastq_reads_old(fn, sv_reads):
    read_len = 0
    fq_recs = {}
    f = open(fn,'r')
#  fq_recs = list(SeqIO.parse(f,'fastq'))
    for header,seq,qual in FastqFile(fn): 
        qname = header.lstrip("@")
        if qname in sv_reads: 
            oseq, sc_seqs, clip_coords = sv_reads[qname]
            cleaned_seq = seq 
            old_seq = oseq.seq
            add = True
            if str(cleaned_seq) != str(old_seq) and sc_seqs:
                idx = old_seq.find(cleaned_seq)
                trimmed_seq = ''
                if idx == 0: 
                    trimmed_seq = old_seq[len(cleaned_seq):len(old_seq)]
                else: trimmed_seq = old_seq[0:idx]    
                sc_lens = 0
                for sc_seq in sc_seqs: 
                    sc_lens += len(sc_seq)
                    if trimmed_seq.find(sc_seq) > -1: 
                        add = False
                if len(cleaned_seq) == (len(old_seq) - sc_lens):
                    for sc_seq in sc_seqs: 
                        if cleaned_seq.find(sc_seq) == -1:
                            # Don't add, just trimmed clipped portion.
                            add = False
#    else: print qname, 'not in sv reads'
        if add:
            fr = fq_read(header, seq, qual)
            read_len = max(read_len, len(fr.seq))   
            fq_recs[fr.id] = fr
    return fq_recs, read_len


def load_kmers(fns, kmers):
        if not fns:
                return kmers

        fns = fns.split(",")
        for fn in fns:
                f = open(fn, 'rU')
                for line in f.readlines():
                        line = line.strip()
                        mer, count = line.split()
                        if mer not in kmers:
                                kmers[mer] = 0
                        kmers[mer] += int(count)
        return kmers


def setup_rmask_all(rmask_fn):
    logger = logging.getLogger('root')

    rmask = {}
    f = open(rmask_fn,'rU')
    flines = f.readlines()
    for line in flines:
        line = line.strip()
        rchr,rbp1,rbp2,rname = line.split("\t")[0:4]
        rchr = rchr.replace('chr','')
        if rchr not in rmask: rmask[rchr] = []
        rmask[rchr].append((rchr,int(rbp1),int(rbp2),rname))
    return rmask


def setup_rmask(gene_coords, ref_path, rmask_fn):
    logger = logging.getLogger('root')
    chrom, s, e, name, intervals = gene_coords
    mask_out_fn = os.path.join(ref_path,name+'_rep_mask.bed')
    marker_fn = get_marker_fn(mask_out_fn)

    rmask = []
    if not os.path.isfile(marker_fn):
        fout = open(mask_out_fn,'w')
        f = open(rmask_fn,'rU')
        flines = f.readlines()
        for line in flines:
            line = line.strip()
            rchr,rbp1,rbp2,rname = line.split("\t")[0:4]
            rchr = rchr.replace('chr','')
            if rchr == chrom:
                if int(rbp1) >= int(s) and int(rbp2) <= int(e): 
                    fout.write("\t".join([str(x) for x in [rchr,int(rbp1),int(rbp2),rname]])+"\n")
                    rmask.append((rchr,int(rbp1),int(rbp2),rname))
        f.close()
        fout.close()

        cmd = 'touch %s'%marker_fn
        p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
        output, errors = p.communicate()  
        logger.info('Completed writing repeat mask file %s, touching marker file %s'%(mask_out_fn,marker_fn))
    else:
        rep_f = open(mask_out_fn,'rU')
        rep_flines = rep_f.readlines()
        for line in rep_flines:
            line = line.strip()
            rchr,rbp1,rbp2,rname = line.split()
            rmask.append((rchr,int(rbp1),int(rbp2),rname))
    return rmask


def extract_refseq_fa(gene_coords, ref_path, ref_fa, direction, target_fa_fn):
    logger = logging.getLogger('root')

    chrom, s, e, name, intervals = gene_coords
#  fa_fn = os.path.join(ref_path, name+'_'+direction+'_refseq.fa')
    marker_fn = get_marker_fn(target_fa_fn) 

    if not os.path.isfile(marker_fn):
        ref_d = SeqIO.to_dict(SeqIO.parse(ref_fa, 'fasta'))
        seq_str = ''
        seq = ref_d[chrom].seq[(s-200):(e+200)]
        if direction == "reverse":
            seq_str = str(seq.reverse_complement())
        else:
            seq_str = str(seq)
        fa = open(target_fa_fn,'w')
        fa.write(">"+name+"\n"+seq_str+"\n")
        fa.close()
        cmd = 'touch %s'%marker_fn
        p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
        output, errors = p.communicate()
        logger.info('Completed writing refseq fasta file %s, touching marker file %s'%(target_fa_fn, marker_fn))
    else:
        logger.info('Refseq sequence fasta (%s) exists already'%target_fa_fn)
    return target_fa_fn 


def seq_trim(qual_str, min_qual):
    counter = 0
    while ord(qual_str[counter])-33 < min_qual:
        counter += 1
        if counter == len(qual_str): break
    return counter  


def get_seq_readname(read):
    end = '1'
    if read.is_read2: end = '2'
    return read.qname + "/" + end


def trim_coords(qual_str, min_qual):
    q = []
    coords = [0,len(qual_str)]
    start = seq_trim(qual_str, min_qual)
    if start == len(qual_str):
        return (0,0,0)
    else:
        end = len(qual_str) - seq_trim(qual_str[::-1], min_qual)
        lngth = end - start
        return (start,end,lngth)


def trim_qual(read, min_qual, min_len):
    qual_str = read.qual
    q = []
    coords = [0,len(qual_str)]
    start = seq_trim(qual_str, min_qual)
    if start == len(qual_str):
        return None
    else:
        end = len(qual_str) - seq_trim(qual_str[::-1], min_qual)
        lngth = end - start
        if lngth < min_len: return None
        nseq = read.seq[start:end]
        nqual = qual_str[start:end]
        read.seq = nseq
        read.qual = nqual
        return read


def fq_line(read, indel_only, min_len, trim=True):
    add_val = '0'
    if indel_only: add_val = '1'
    lineout = None
    if trim: read = trim_qual(read, 5, min_len)
    if read: 
        lineout = "@"+get_seq_readname(read) + "_" + add_val + "\n" + read.seq + "\n+\n" + read.qual + "\n"
    return lineout


def get_overlap_index_nomm(a,b):
    i = 0
    while a[i:] != b[:len(a[i:])]:
        i += 1
    return i


def seq_complexity(seq, N):
    nmers = {}
    total_possible = len(seq) - 2
    for i in range(len(seq) - (N - 1)):
        nmers[str(seq[i:i+N]).upper()] = True
    complexity = round((float(len(nmers))/float(total_possible))*100,4)
    return complexity


def get_overlap_index_mm(a, b):
    i = 0
    nmismatch = [0,0]
    match = False
    while i < len(a) and not match:
        nmismatch = [0,0]
        c = 0
        match_len = min(len(a[i:]), len(b[:len(a[i:])]))
        for aa,bb in zip(a[i:],b[:len(a[i:])]):
            if aa != bb: 
                nmismatch[0] += 1
                nmismatch[1] += 1
            else:
                nmismatch[0] = 0
            if nmismatch[0] > 1 or nmismatch[1] > 3: 
                break  
            c += 1
#    print c, match_len, i, c== match_len, a[i:], b[:len(a[i:])]
        if c == match_len: 
            match = True
        i += 1
    return i-1


def get_read_kmers(seq,l,skmers):
    kmers = []
    i = 0
    while (i+l) <= len(seq):
        k = seq[i:i+l]
        #if k in skmers:
        kmers.append(k)
        i += 1
    return list(set(kmers)&set(skmers))


def get_overlap_index(a,b):
    i = 0
    nmismatch = 10
    while nmismatch > 1:
        nmismatch = 0
        for aa,bb in zip(a[i:],b[:len(a[i:])]):
            if aa != bb: nmismatch += 1
        i += 1
#  while a[i:] != b[:len(a[i:])]:
#    i += 1
    return i-1


def server_ready(f):
    logger = logging.getLogger('root')
    while not os.path.exists(f):
        logger.info("Waiting for log file %s"%f)
        time.sleep(10)
    ready = False
    f = open(f,'r')
    flines = f.readlines()
    for line in flines: 
        if line.find('Server ready for queries') > -1: ready = True
    return ready


class fq_read:
        def __init__(self, header, seq, qual, indel_only):
                self.id = header
                self.seq = str(seq)
                self.qual = str(qual)
                self.used = False
                self.dup = False
                self.indel_only = indel_only


class FastqFile(object):
    def __init__(self,f):
        if isinstance(f,str):
            f = open(f)
            self._f = f
    def __iter__(self):
        return self

    def next(self):
        header, seq, qual_header, qual = [self._f.next() for _ in range(4)]
        header = header.strip()
        inst,lane,tile,x,y_end = header.split(':')
        seq = seq.strip()
        qual = qual.strip()
        bc = None
        y = y_end
        if y.find('/') > -1:
            y, end = y.split('/')
        if y.find('#') > -1:
            y, bc = y.split('#')
        header_dict = {'inst':inst,
                                     'lane':int(lane),
                                     'tile':int(tile),
                                     'x':int(x),
                                     'y':int(y),
                                     'end':end,
                                     'bc': bc}
        return (header,seq,qual)


class Annotation:
    def __init__(self):
        self.genes = {}
        self.logging_name = 'breakmer.utils.Annotation'

    def add_genes(self, gene_fn):
        log(self.logging_name, 'info', 'Adding gene annotations from %s'%gene_fn)
        gene_f = open(gene_fn,'r')
        gene_flines = gene_f.readlines()
        for line in gene_flines[1:]:
            line = line.strip()
            linesplit = line.split()
            chrom = linesplit[2]
            start = int(linesplit[4])
            end = int(linesplit[5])
            geneid = linesplit[12]
            if geneid in self.genes:
                if start <= self.genes[geneid][1] and end >= self.genes[geneid][2]: self.genes[geneid] = [chrom,start,end]
            else: self.genes[geneid] = [chrom,start,end]
        gene_f.close()

    def add_regions(self, regions_bed_fn):
        region_f = open(regions_bed_fn, 'rU')
        region_lines = region_f.readlines()
        for line in region_lines:
            line = line.strip()
            chrom, start, end, name = line.split()
            if name not in self.genes: self.genes[name] = [chrom,int(start),int(end)]
        log(self.logging_name, 'info', 'Adding in %d other target regions'%len(region_lines))
    
    def set_gene(self, chrom, pos):
        ann_genes = []
        if chrom.find('chr') == -1: chrom = 'chr'+str(chrom) 
        for g in self.genes:
            gs = self.genes[g][1]
            ge = self.genes[g][2]
            if chrom == self.genes[g][0]:
                if len(pos) == 1:
                    if int(pos[0]) >= gs and int(pos[0]) <= ge:
                        ann_genes.append(g)
                        break
                else:
                    # Find genes between pos1 and pos2
                    if (int(pos[0]) >= gs and int(pos[0]) <= ge) or (int(pos[1]) >= gs and int(pos[1]) <= ge):
                        ann_genes.append(g)
        if len(ann_genes) == 0: ann_genes = ['intergenic']
        return ",".join(ann_genes)
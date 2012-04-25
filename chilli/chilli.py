#!/usr/bin/env python
# -*- coding:utf-8 -*-
from __future__ import division
'''Usefull function collected by Wubin Qu


by Wubin Qu <quwubin@gmail.com>,
Copyright @ 2010, All Rights Reserved.
'''

Author  = 'Wubin Qu  <quwubin@gmail.com>'
Date    = 'Jun-27-2011 18:55:37'
Version = '1.0'

import sys
import os
import argparse

D2n_dic = dict(A=0, T=3, C=2, G=1, a=0, t=3, c=2, g=1)
n2D_dic = {0:'A', 3:'T', 2:'C', 1:'G', 0:'a', 3:'t', 2:'c', 1:'g'}

def get_opt():
    '''Get options'''
    global args
    parser = argparse.ArgumentParser(description='Descriptions about the script.')
    parser.add_argument('-i', '--infile', nargs='?', type=argparse.FileType('r'),
		    default=sys.stdin, help='Input file to be processed', required=True)
    parser.add_argument('-o', '--outfile', nargs='?', type=argparse.FileType('w'),
		    default=sys.stdout, help='Output file name for storing the results')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
    args = parser.parse_args()

def DNA2int(seq):
    '''convert a sub-sequence/seq to a non-negative integer'''
    result = 0
    seq_len = len(seq)
    for i, letter in enumerate(seq):
        result += D2n_dic[letter] * 4 ** (seq_len - i - 1)
    return result

def DNA2int_2_strand(seq):
    '''convert a sub-sequence/seq to a non-negative integer'''
    plus_mer = 0
    minus_mer = 0
    length = len(seq) - 1 
    for i, letter in enumerate(seq):
        plus_mer += D2n_dic[letter] * 4 ** (length - i)
        minus_mer += (3 - D2n_dic[letter]) * 4 ** i

    return plus_mer, minus_mer

def baseN(num, b):
    '''convert non-negative decimal integer n to
    equivalent in another base b (2-36)'''
    return ((num == 0) and  '0' ) or ( baseN(num // b, b).lstrip('0') + "0123456789abcdefghijklmnopqrstuvwxyz"[num % b])

def int2DNA(num, k):
    ''''''
    seq = baseN(num ,4)
    return 'A'* (k-len(seq))+(''.join([n2D_dic[int(base)]  for base in seq]))

def session(parent_dir=None):
    '''Create session directory'''
    import sys
    import uuid
    if not parent_dir:
	if sys.platform == 'win32':
	    print 'Not support now'
	    exit()
	else:
	    parent_dir = '/tmp'

    session_id = str(uuid.uuid4())
    session_dir = os.path.join(parent_dir, session_id)
    os.mkdir(session_dir)
    os.chmod(session_dir, 0777)
    return session_dir

def print_seq(seq, width=80, linesep=''):
    '''Print biological sequence with constant width, such as 80'''
    import os

    line = ''
    start = 0
    line_list = []
    width = int(width)
    if not linesep:
        linesep = os.linesep
    while start <= len(seq):
        seq_line = seq[start : (start + width)]
        start = start + width
        line_list.append(seq_line)

    line = linesep.join(line_list)
    return line 

def format_show_time(unformated_time=False):
    '''Format show time'''
    import time
    ISOTIMEFORMAT = '%Y-%m-%d %X'
    if not unformated_time:
        return time.strftime(ISOTIMEFORMAT)
    else:
        return time.strftime(ISOTIMEFORMAT, unformated_time)

def seconds2min_sec(seconds):
    '''Format the senconds as Minutes:Seconds'''

    if seconds < 1:
        return '%.2f' % seconds

    m, s = divmod(int(seconds), 60)

    return m, s

def cal_dimer(seq1, seq2, mv=50, dv=1.5, d=50, n=0.25, align_mode='ANY', thermo_path='lib/primer3_config/'):
    '''Run ntthal'''
    import subprocess
    cmd = 'ntthal -mv %s -dv %s -d %s -n %s -s1 %s -s2 %s -a %s -path %s' \
            % (mv, dv, d, n, seq1, seq2, align_mode, thermo_path)

    out, align = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    try:
        fields = out.splitlines()[0].split()
        delta_g = float(fields[-4]) / 1000
        tm = float(fields[-1])
        align = '\n'.join(out.splitlines()[1:])
    except:
        return False, False, False

    return delta_g, tm, align

def cal_hairpin(seq, mv=50, dv=1.5, d=50, n=0.25, align_mode='HAIRPIN', thermo_path='lib/primer3_config/'):
    '''Run ntthal'''
    import subprocess
    cmd = 'ntthal -mv %s -dv %s -d %s -n %s -s1 %s -a %s -path %s' \
            % (mv, dv, d, n, seq, align_mode, thermo_path)

    out, align = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    try:
        fields = out.splitlines()[0].split()
        delta_g = float(fields[-4]) / 1000
        tm = float(fields[-1])
        align = '\n'.join(out.splitlines()[1:])
    except:
        return False, False, False

    return delta_g, tm, align

def cal_tm(seq, mv=50, dv=1.5, d=50, n=0.25, tp=1, sc=1):
    '''Run oligotm'''
    import subprocess
    cmd = 'oligotm -mv %s -dv %s -d %s -n %s -tp %s -sc %s %s' \
            % (mv, dv, d, n, tp, sc, seq)
    out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

    if err:
        print err
        exit()

    return float(out.strip())

def run_blastn(query, subject):
    ''''''
    import subprocess
    cmd = 'blastn -task blastn -query %s -subject %s -outfmt=7' % (query, subject)
    out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

    if err:
        print err
        exit()

    return out

def run_blastdbcmd(gene_id, db_list, output, start=None, stop=None):
    '''Run blastdbcmd'''
    import subprocess
    uniq = uuid.uuid4()
    if output == 'file':
        if not (start and stop):
            out_file = os.path.join('session_dir', '%s.txt' % (gene_id))
        else:
            out_file = os.path.join('session_dir', '%s_%s_%s.txt' % (gene_id, start, stop))

        if os.path.isfile(out_file):
            return out_file, ''

    else:
        out_file = ''



    for db in db_list:
        db = '/opt/ncbi_db/%s' % db
        if not (start and stop):
            cmd = 'bin/blastdbcmd -db %s -entry %s' % (db, gene_id)
        else:
            cmd = 'bin/blastdbcmd -db %s -entry %s -range %s-%s' % (db, gene_id, start, stop)

        if out_file:
            cmd = '%s -out %s' % (cmd, out_file)

        seq, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        if err:
            continue
        else:
            break

    if out_file:
        return out_file, err
    else:
        return seq, err

def detecting_consecutive_integers(integer_list):
    '''Return consecutive integer list'''
    import itertools
    import operator

    results = []

    for k, g in itertools.groupby(enumerate(integer_list), lambda (i,x):i-x):
        c_i = map(operator.itemgetter(1), g)
        if len(c_i) > 1:
            results.append('%s-%s' % (c_i[0], c_i[-1]))
        else:
            results.append(str(c_i[0]))

    return results

def parse_blast7_best(fh):
    '''Parse BLASTN/BLASTP outfmt "7" output file'''
    if not isinstance(fh, file):
        fh = fh.splitlines()

    import re
    records = {}
    for line in fh:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        fields = line.split()
        query_id = fields[0]
        query_id = re.sub('lcl\|', '', query_id)
        sbjct_id = fields[1]

        if query_id == sbjct_id:
            continue

        # Here is an example:
        # Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
        if query_id not in records:
	    records[query_id] = {
		    'sbjct_id' : sbjct_id,
		    'identity' : float(fields[2]),
		    'aln_len' : int(fields[3]),
		    'mismatches' : int(fields[4]),
		    'gap_opens' : int(fields[5]),
		    'q_start' : int(fields[6]),
		    'q_end' : int(fields[7]),
		    's_start' : int(fields[8]),
		    's_end' : int(fields[9]),
		    'evalue' : float(fields[10]),
		    'bit' : float(fields[11]),
		    }
	else:
	    pass
        

    return records

def parse_blastn7(fh):
    '''Parse BLASTN/BLASTP, BLAT(BLAST8) outfmt "7" output file'''
    if not isinstance(fh, file):
        fh = fh.splitlines()

    import re
    records = {}
    for line in fh:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        fields = line.split()
        query_id = fields[0]
        query_id = re.sub('lcl\|', '', query_id)
        sbjct_id = fields[1]

        if query_id == sbjct_id:
            continue

        if query_id not in records:
            records[query_id] = {}
        
        if sbjct_id not in records[query_id]:
            records[query_id][sbjct_id] = []

        # Here is an example:
        # Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
        records[query_id][sbjct_id].append(
	    {
		'identity' : float(fields[2]),
		'aln_len' : int(fields[3]),
		'mismatches' : int(fields[4]),
		'gap_opens' : int(fields[5]),
		'q_start' : int(fields[6]),
		'q_end' : int(fields[7]),
		's_start' : int(fields[8]),
		's_end' : int(fields[9]),
		'evalue' : float(fields[10]),
		'bit' : float(fields[11]),
	    }
	)

    return records

def parse_blastn7_to_array(fh):
    '''Parse BLASTN/BLASTP, BLAT(BLAST8) outfmt "7" output file'''
    if not isinstance(fh, file):
        fh = fh.splitlines()

    import re
    records = {}
    for line in fh:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        fields = line.split()
        query_id = fields[0]
        query_id = re.sub('lcl\|', '', query_id)
        sbjct_id = fields[1]

        if query_id == sbjct_id:
            continue

        if query_id not in records:
            records[query_id] = {}
        
        # Here is an example:
        # Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
        records[query_id].append(
        {
        'sbjct_id' : sbjct_id,
        'identity' : float(fields[2]),
        'aln_len' : int(fields[3]),
        'mismatches' : int(fields[4]),
        'gap_opens' : int(fields[5]),
        'q_start' : int(fields[6]),
        'q_end' : int(fields[7]),
        's_start' : int(fields[8]),
        's_end' : int(fields[9]),
        'evalue' : float(fields[10]),
        'bit' : float(fields[11]),
        }
    )

    return records

def cal_coverage(range_list):
    '''Cal the coverage region for A seq blastn against B seq'''
    from operator import itemgetter

    positions = []
    coverage = []
    for start, end in range_list:
        positions.append(start)
        positions.append(end)

    positions = list(set(positions))
    positions.sort()
    #print(positions)

    for i in xrange(len(positions) - 1):
        for start, end in range_list:
            if positions[i] >= start and positions[i+1] <= end:
                coverage.append([positions[i], positions[i+1]])
                break

    coverage = list(set([str(item) for item in coverage]))
    coverage = [eval(item) for item in coverage]
    coverage.sort(key = itemgetter(0))
    #print coverage

    if len(coverage) <= 1:
        merge_coverage = coverage
    else:
        merge_coverage = []
        merge_coverage.append(coverage[0])
        for i in xrange(1, len(coverage)):
            if merge_coverage[-1][1] == coverage[i][0]:
                merge_coverage[-1][1] = coverage[i][1]
            else:
                merge_coverage.append(coverage[i])

    return merge_coverage

def print2stderr(msg):
    '''Print to sys.stderr'''
    import sys

    print >> sys.stderr, msg

def perm_unique(elements):
    '''permutations with unique values
    http://stackoverflow.com/questions/6284396/permutations-with-unique-values 
    '''
    eset=set(elements)
    listunique = [[i,elements.count(i)] for i in eset]
    u=len(elements)
    return perm_unique_helper(listunique,[0]*u,u-1)

def perm_unique_helper(cl,lst,d):
    '''Sub-function of perm_unique'''
    if d < 0:
	yield list(lst)
    else:
	for (ind,(j,i)) in enumerate(cl):
	    if i>0:
		lst[d]=j
		cl[ind][1]-=1
		for g in  perm_unique_helper(cl,lst,d-1):
		    yield g

                cl[ind][1]+=1

def cal_GC_content(seq, length=None):
    '''Calculate the GC content of a sequence'''

    seq = seq.lower()
    if not length:
        length = len(seq)

    c_num = seq.count('c')
    g_num = seq.count('g')
    n_num = seq.count('n')

    GC_content = (c_num + g_num + n_num * 0.5) / length * 100
    return GC_content

def sn_creator_with_drop_back(sn_cutoff_list):
    '''sn_list creator, return sn_list_iter
    2011-12-12 0:04
    '''
    import random
    max_sn_list = 1
    for sn_cutoff in sn_cutoff_list:
        max_sn_list *= sn_cutoff

    sn_list_set = set()

    run_num = 0
    seed_cutoff = 1
    while True:
        sn_list = []
        seed_grow_threshold = 1
        for index, cutoff in enumerate(sn_cutoff_list):
            sn_list.append(random.randint(1, min(cutoff, seed_cutoff)))
            seed_grow_threshold *= min(cutoff, seed_cutoff)

        if str(sn_list) in sn_list_set:
            continue

        sn_list_set.add(str(sn_list))
        run_num += 1
        if run_num >= seed_grow_threshold:
            seed_cutoff += 1

        if run_num >= max_sn_list:
            return

        yield [x-1 for x in sn_list]

def sn_creator_no_drop_back(sn_cutoff_list):
    '''sn_list creator, return sn_list_iter
    2011-12-12 0:04
    '''
    sn_list = []
    while True:
        if len(sn_list) < 1:
            sn_list = [0] * len(sn_cutoff_list)
            yield sn_list
            continue

        for index, cutoff in enumerate(sn_cutoff_list):
            sn_list[index] += 1
            if sn_list[index] > cutoff - 1:
                return

        yield sn_list

def set_cache(data, filename):
    '''Filename should be full path'''
    import os
    import shelve

    d = shelve.open(filename)
    d[os.path.basename(filename)] = data
    d.close()

def get_cache(filename):
    '''Filename should be full path'''
    import os
    import shelve

    d = shelve.open(filename)
    data = d[os.path.basename(filename)]
    d.close()

    return data

def run_primer3(p3_file, primer3_core=None):
    '''Run primer3_core (V2.3.0'''
    import subprocess
    if not primer3_core:
        print2stderr("primer3_core need set first")

    cmd = '%s < %s' % (primer3_core, p3_file)
    
    try:
        out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    except:
        print2stderr('Error: run primer3_core')

    if err:
        print2stderr(err)

    return out

def get_align(a, b):
    '''Align ||||'''
    if len(a) == 0 or len(b) == 0:
        return ''

    if len(a) != len(b):
        print2stderr("%s and %s don't have the same length" % (a, b))

    align = []
    for i in range(len(a)):
        if a[i] == b[i]:
            align.append('|')
        else:
            align.append('x')

    return ''.join(align)

def random_string(length):
    import string
    import random
    rstr = ''.join(random.choice(string.letters) for i in xrange(length))
    return rstr


def main ():
    '''Main'''
    #get_opt()
    a = detecting_consecutive_integers([1, 2, 4, 5, 6, 7, 11, 13, 14, 15])
    print a

if __name__ == '__main__':
    main()


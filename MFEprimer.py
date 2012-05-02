#!/usr/bin/env python
# -*- coding: utf-8 -*- #   
from __future__ import division
'''MFEprimer: Multiple factor evaluation of the specificity of PCR primers'''

Program = 'MFEprimer'
Author = 'Wubin Qu, Yang Zhou & Chenggang Zhang, BIRM, China'
AuthorEmail = 'Wubin Qu <quwubin@gmail.com>'
Date = 'Nov-8-2011'
License = 'Please contact Wubin Qu <quwubin@gmail.com>'
Version = 'V2.0'

import sys, os
import time
import math
import textwrap
import subprocess
import platform
import argparse
from operator import itemgetter
import sqlite3
from pprint import pprint

from chilli import Seq
from chilli import SeqCheck
from chilli import TmDeltaG
from chilli import chilli
from chilli import FastaFormatParser
from chilli import DegenerateSeqConvetor

src_path = os.path.split(os.path.realpath(sys.argv[0]))[0]
bin_path = os.path.join(src_path, 'bin', platform.architecture()[0])

global degenerate
degenerate = 'no'

nn_mm_data = {
    'GA' : ['GC'],
    'GG' : ['GT', 'GC'],
    'GT' : ['GC'],

    'CG' : ['CC', 'CA'],
    'CT' : ['CC', 'CA'],

    'AG' : ['AC'],
    }

def get_opt():
    '''Check and parsing the opts'''
    parser = argparse.ArgumentParser(prog='MFEprimer', description='MFEprimer: A fast and thermodynamics-based PCR primer specificity checking program.', usage='%(prog)s.py [options] -i primers.fasta -d Human.genomic.fasta')
    parser.add_argument('-i', '--infile', nargs='?', type=argparse.FileType('r'),
	    default=sys.stdin, help='[Required] Primer sequences for specificity checking. [File]', required=True)
    parser.add_argument('-d', '--database', nargs='+', type=str, 
	    help='[Required] Database name for specificity checking.', required=True)
    parser.add_argument('-o', '--outfile', nargs='?', type=argparse.FileType('w'),
	    default=sys.stdout, help='[Optional] Output file name for storing the results, default is screen. [File]', required=False)

    parser.add_argument('-k', '--k_value', nargs='?', type=int,
	    default=9, help='[Optional] K value, must be identical to the k value when indexing the database. [Integer]', required=False)

    parser.add_argument('--amplicon', action='store_true',
	    help='[Optional] Produce the amplicons sequence in Fasta format, only works for normal output format (not tabular).')

    parser.add_argument('--tab', action='store_true',
	    help='[Optional] Output in tabular format.')

    parser.add_argument('-v', '--version', action='version', version='%(prog)s 2.0')

    filter_group = parser.add_argument_group('Results filter settings:' ,'set these arguments for filtering the results.')

    filter_group.add_argument('--ppc', nargs='?', type=float, default=30, help='[Optional] Lower limit of the PPC. Legal value: [0~99]. Default = 30. [Float]', required=False)

    filter_group.add_argument('--size_start', nargs='?', help='[Optional] Lower limit of the expected amplicon size range in bp, default = 50. [Integer]', default=50, type=int)
    filter_group.add_argument('--size_stop', nargs='?', help='[Optional] Upper limit of the expected amplicon size range in bp, default = 5000. [Integer]', default=2000, type=int)

    filter_group.add_argument('--tm_start', nargs='?', help='[Optional] Lower limit of the Melting temperature (Tm) [Celsius degree], default = 0. [Float]', default=0, type=float)
    filter_group.add_argument('--tm_stop', nargs='?', help='[Optional] Upper limit of the Melting temperature (Tm) [Celsius degree], default = unlimited. [Float]', default=sys.maxint, type=float)

    filter_group.add_argument('--dg_start', nargs='?', help='[Optional] Lower limit of the Gibbs free energy for the last 5 resides at 3\' end of the primer [kcal/mol], default = unlimited. [Float]', default=-sys.maxint, type=float)
    filter_group.add_argument('--dg_stop', nargs='?', help='[Optional] Upper limit of the Gibbs free energy for the last 5 resides at 3\' end of the primer [kcal/mol], default = 0. [Float]', default=0, type=float)

    thermo_group = parser.add_argument_group('Experimental settings:' ,'set these arguments for your real experimental values.')
    thermo_group.add_argument('--mono_conc', nargs='?', type=float,
	    default=50, help='[Optional] Concentration of monovalent cations [mM], default = 50.0 [Float]', required=False)
    thermo_group.add_argument('--diva_conc', nargs='?', type=float,
	    default=1.5, help='[Optional] Concentration of divalent cations [mM], default = 1.5 [Float]', required=False)
    thermo_group.add_argument('--oligo_conc', nargs='?', type=float,
	    default=50, help='[Optional] Concentration of annealing oligos [nM], default = 50.0 [Float]', required=False)
    thermo_group.add_argument('--dntp_conc', nargs='?', type=float,
	    default=0.25, help='[Optional] Concentration of dNTPs [nM], default = 0.25 [Float]', required=False)

    options = parser.parse_args()
    if options.ppc < 0 or options.ppc > 100:
	print 'Error: Illegal value for ppc'
	exit()

    if options.size_start > options.size_stop or options.size_start < 0:
	print 'Illegal value for size_start or size_stop'
	exit()

    if options.tm_start > options.tm_stop or options.tm_start < 0:
	print 'Illegal value for tm_start or tm_stop'
	exit()

    if options.dg_start > options.dg_stop or options.dg_start > 0:
	print 'Illegal value for dg_start or DeltagG_stop'
	exit()

    return options

def print_head(out, options):
    '''Format head information for the output'''
    linesep = os.linesep

    out.append('%s %s [%s]%s' % (Program, Version, Date, linesep))
    #out.append(''*4 + Author + linesep * 2)

    reference = 'Wubin Qu, Zhiyong Shen, Dongsheng Zhao, Yi Yang and Chenggang Zhang. (2009) MFEprimer: multiple factor evaluation of the specificity of PCR primers, Bioinformatics, 25(2), 276-278'

    reference = textwrap.fill(reference, 80)
    out.append('%s.%s' % (reference, linesep*2))

    return out

def primer_analysis(product, options, oligos, session_dir, fcdict, db):
    '''Analysis the candidate forward and reverse primer and check whether they can amplify an amplicon'''
    mid_seq_id_list = []
    tmp_list = []
    amp_list = []
    filter_product = []

    # Forward primer
    for i in xrange(len(product)):
        # Reverse primer
        #print i
        amp = product[i]
        hid = amp['hid']
        pid = amp['pid']
        mid = amp['mid']
        f_len = amp['plen']
        r_len = amp['mlen']
        pseq = amp['pseq']
        mseq = amp['mseq']
        size = amp['size']
        f_3_pos = amp['f3_pos'] 
        r_3_pos = amp['r3_pos']

        p_qseq = amp['p_qseq']
        p_aseq = amp['p_aseq']
        p_sseq = amp['p_sseq']
        p_tail = amp['p_tail']

        m_qseq = amp['m_qseq']
        m_aseq = amp['m_aseq']
        m_sseq = amp['m_sseq']
        m_tail = amp['m_tail']

        p_Tm = amp['p_Tm']
        p_DeltaG = amp['p_DeltaG']
        m_Tm = amp['m_Tm']
        m_DeltaG = amp['m_DeltaG']
        #print 2*i, 2*i + 1
        p_3_DeltaG = TmDeltaG.calDeltaG(p_qseq[-5:], Seq.complement(p_sseq[-5:]), mono_conc=options.mono_conc, diva_conc=options.diva_conc, dntp_conc=options.dntp_conc)
        m_3_DeltaG = TmDeltaG.calDeltaG(m_qseq[:5], Seq.complement(m_sseq[:5]), mono_conc=options.mono_conc, diva_conc=options.diva_conc, dntp_conc=options.dntp_conc)

        # Filter DeltaG
        if p_3_DeltaG < float(options.dg_start) or p_3_DeltaG > float(options.dg_stop):
            continue
        if m_3_DeltaG < float(options.dg_start) or m_3_DeltaG > float(options.dg_stop):
            continue

        ppc = cal_PPC(len(p_qseq), f_len, len(m_qseq), r_len)
        # Filter by PPC
        if ppc < options.ppc:
            continue
        
        mid_seq_id = '%s:%s-%s' % (hid, f_3_pos, r_3_pos)
        mid_seq_id_list.append(mid_seq_id)


        ave_Tm = (p_Tm + m_Tm) / 2 # For sort
        to_be_added = (ave_Tm, ppc, p_3_DeltaG, m_3_DeltaG)
        tmp_list.append(to_be_added)
        filter_product.append(amp)

    mid_seq_list = get_mid_seq(mid_seq_id_list, options, session_dir, db)
    for i in xrange(len(mid_seq_list)):
        mid_seq = mid_seq_list[i]
        (ave_Tm, ppc, p_3_DeltaG, m_3_DeltaG) = tmp_list[i]
        amp = filter_product[i]
        pid = amp['pid']
        mid = amp['mid']

        hid = int(amp['hid'])
        real_hid = fcdict[str(hid)]['id']
        hdesc = fcdict[str(hid)]['desc']
        amp_graphic = draw_graphical_alignment_primer(amp, oligos, options, mid_seq)
        size = amp['size']
        amp['p_3_DeltaG'] = p_3_DeltaG
        amp['m_3_DeltaG'] = m_3_DeltaG
        amp['real_hid'] = real_hid
        amp['hdesc'] = hdesc
        amp['mid_seq'] = mid_seq
        amp['amp_graphic'] = amp_graphic
        amp_list.append([ave_Tm, ppc, size, amp])

    return amp_list

def cal_PPC(f_match, p_len, r_match, m_len):
    '''Cal PPC parameter'''
    ave = (f_match + r_match) / 2
    stdev = math.sqrt((f_match - ave)**2 + (r_match - ave)**2)
    cv = 1 - stdev/ave
    ppc = f_match / p_len * r_match / m_len * cv * 100

    return ppc

def format_output_primer(amp_list, oligos, options, start_time, session_dir):
    '''Format output in primer task'''
    linesep = os.linesep

    out = []

    out = print_head(out, options)

    ID_list = []
    for i in xrange(len(oligos)):
        ID_list.append(oligos[i]['id'])

    query_line = textwrap.fill('Query = %s' % ('; '.join(ID_list)), 80)
    out.append(query_line)
    out.append('        %s primer sequences' % (len(oligos)))

    out.append(linesep)

    out.append('Database = %s' % textwrap.fill(', '.join([os.path.basename(db) for db in options.database]), 80))
    #out.append('        %s sequences' % (len(fcdict)))
    out.append(linesep)

    out.append('Reports Begining'.ljust(80, '.'))
    out.append(linesep * 2)

    amp_num = len(amp_list)
    if amp_num > 1:
        out.append('Distribution of %s potential PCR amplicons predicted by MFEprimer-2.0 on the query primers' % amp_num)
    else:
        out.append('Distribution of %s potential PCR amplicon predicted by MFEprimer-2.0 on the query primers' % amp_num)

    out.append(linesep)
    out.append('[Sorted by average Tm in descending order]')
    out.append('FP '.rjust(69) + 'RP '.rjust(8) + 'FP '.rjust(8) + 'RP '.rjust(8) + 'FP '.rjust(8) + 'RP '.rjust(8))
    # Δ takes two characters position
    out.append('Size'.rjust(53) + 'PPC '.rjust(8) + 'Tm '.rjust(8) + 'Tm '.rjust(8) + ('%sG' % u'\u0394').rjust(7) + ('%sG' % u'\u0394').rjust(7) + ('3\'%sG' % u'\u0394').rjust(8) + ('3\'%sG' % u'\u0394').rjust(7))
    out.append('Primers producing potential PCR products:'.ljust(42) + '(bp)'.rjust(11) + '(%) '.rjust(8) + u'\u2103'.rjust(6) + u'\u2103'.rjust(7) + '(kcal/mol)'.center(22) + '(kcal/mol)'.center(14))
    out.append(linesep)

    detail_line = []
    fa_file = []
    sn = 0
    #amp_list.append([amp_len, ave_Tm, p, m, ppc, amp_graphic, mid_seq, real_hid, hdesc])
    amp_list.sort(key=itemgetter(1, 2), reverse=True)
    for ave_Tm, ppc, amp_len, amp in amp_list:
        sn = sn + 1
        hid = amp['real_hid']
        desc = '%s: %s' % (sn, hid)

        amp_len = amp['size']

        p_qid = amp['pid']
        f_len = amp['plen']
        pseq = amp['pseq']
        f_3_pos = amp['f3_pos'] 
        p_3_DeltaG = amp['p_3_DeltaG']
        p_qseq = amp['p_qseq']
        p_aseq = amp['p_aseq']
        p_sseq = amp['p_sseq']
        p_tail = amp['p_tail']
        p_Tm = amp['p_Tm']
        p_DeltaG = amp['p_DeltaG']
        p_sb = f_3_pos - len(p_aseq) + 1

        m_qid = amp['mid']
        r_len = amp['mlen']
        mseq = amp['mseq']
        r_3_pos = amp['r3_pos']
        m_3_DeltaG = amp['m_3_DeltaG']
        m_qseq = amp['m_qseq']
        m_aseq = amp['m_aseq']
        m_sseq = amp['m_sseq']
        m_tail = amp['m_tail']
        m_Tm = amp['m_Tm']
        m_DeltaG = amp['m_DeltaG']
        m_se = r_3_pos + len(m_aseq)

        amp_graphic = amp['amp_graphic']
        mid_seq = amp['mid_seq']
        real_hid = amp['real_hid']
        hdesc = amp['hdesc']

        amp_seq = p_tail + p_qseq + mid_seq + m_qseq + m_tail
        amp_GC = chilli.cal_GC_content(amp_seq, + amp_len)

        if len(desc) > 42:
            desc = desc[:42] + '...' 

        if p_qid == m_qid:
            ppc = '-%.1f' % ppc
        else:
            ppc = '%.1f' % ppc

        out.append(desc.ljust(42) + (str(amp_len)).rjust(11) + ppc.rjust(8) + ('%.1f' % p_Tm).rjust(8) + ('%.1f' % m_Tm).rjust(8) + ('%.1f' % p_DeltaG).rjust(8) + ('%.1f' % m_DeltaG).rjust(8) + ('%.1f' % p_3_DeltaG).rjust(8) + ('%.1f' % m_3_DeltaG).rjust(8))

        if not hdesc:
            detail_line.append('%s: %s + %s ==> %s%s' % (sn, p_qid, m_qid, hid, linesep))
            fa_desc = '>%s %s + %s %s' % (sn, p_qid, m_qid, hid)
        else:
            detail_line.append('%s: %s + %s ==> %s %s%s' % (sn, p_qid, m_qid, hid, hdesc, linesep))
            fa_desc = '>%s %s %s %s %s' % (sn, p_qid, m_qid, hid, hdesc)

        detail_line.append('  ' + 'PPC = %s%%, Size = %s bp, GC content = %.1f%%' % (ppc, amp_len, amp_GC))
        detail_line.append('  ' + 'FP: Tm = %.1f (%s), %sG = %.1f (kcal/mol), 3\'%sG = %.1f (kcal/mol)' % (p_Tm, u'\u2103', u'\u0394', p_DeltaG, u'\u0394', p_3_DeltaG))
        detail_line.append('  ' + 'RP: Tm = %.1f (%s), %sG = %.1f (kcal/mol), 3\'%sG = %.1f (kcal/mol)' % (m_Tm, u'\u2103', u'\u0394', m_DeltaG, u'\u0394', m_3_DeltaG))
	detail_line.append('  ' + 'Binding sites: %s(%s/%s) ... %s(%s/%s)' % (p_sb, len(p_aseq), f_len, m_se, len(m_aseq), r_len))
        detail_line.append(linesep)
        detail_line.append(amp_graphic + linesep)
        fa_seq = chilli.print_seq(amp_seq, 80)
        fa_file.append(fa_desc)
        fa_file.append(fa_seq)
        detail_line.append(fa_desc + linesep + fa_seq + linesep)

    #out = []
    out.append(linesep)
    out.append('Details for the primers binding to the DNA template') 
    out.append('[Sorted by average Tm in descending order]' + linesep)
    for i in xrange(len(detail_line)):
        line = detail_line[i]
        out.append(line)

    out.append(linesep*2)
    out = print_foot(out, options, start_time)

    out = os.linesep.join(out)
    options.outfile.write(out.encode('utf-8'))

    if options.amplicon:
        try:
            out_file = options.outfile.name + '.fa'
            fh = open(out_file, 'w')
        except:
            msg = 'Error: can not open %s for write' % out_file
            print2stderr(msg)

        fh.write(os.linesep.join(fa_file))
        fh.close()

def print_foot(out, options, start_time):
    '''Output the foot information in the output file'''
    linesep = os.linesep
    out.append('#' * 80 + linesep)
    elapsed_time = time.time() - start_time
    
    out.append('InFile: '.rjust(42) + os.path.basename(options.infile.name))
    out.append('OutFile: '.rjust(42) + os.path.basename(options.outfile.name))


    if options.amplicon:
        out.append('OutAmpliconFile: '.rjust(42) + '%s.fa' % (options.outfile.name))

    out.append('Database: '.rjust(42) + textwrap.fill(', '.join([os.path.basename(db) for db in options.database]), 80))
    out.append('Degenerate: '.rjust(42) + degenerate.capitalize())
    out.append('k value: '.rjust(42) + str(options.k_value) + linesep)

    out.append('Concentration of monovalent cations [mM]: '.rjust(42) + str(options.mono_conc))
    out.append('Concentration of divalent cations [mM]: '.rjust(42) + str(options.diva_conc))
    out.append('Concentration of annealing oligos [nM]: '.rjust(42) + str(options.oligo_conc))
    out.append('Concentration of dNTPs [mM]: '.rjust(42) + str(options.dntp_conc) + linesep)
    out.append('PPC cutoff [%]: '.rjust(42) + str(options.ppc))
    out.append('Size start [bp]: '.rjust(42) + str(options.size_start))
    out.append('Size stop [bp]: '.rjust(42) + str(options.size_stop))

    out.append(('Tm start [%s]: ' % (u'\u2103')).rjust(41) + str(options.tm_start))

    if options.tm_stop == sys.maxint:
        tm_stop = 'Unlimited'
    else:
        tm_stop = options.tm_stop

    if options.dg_start == -sys.maxint:
        dg_start = 'Unlimited'
    else:
        dg_start = options.dg_start

    out.append(('Tm stop [%s]: ' % (u'\u2103')).rjust(41) + str(tm_stop))
    # It seems that ΔG only occupied 1.5 space
    out.append(('%sG start [kcal/mol]: ' % (u'\u0394')).rjust(41) + str(dg_start))
    out.append(('%sG stop [kcal/mol]: ' % (u'\u0394')).rjust(41) + str(options.dg_stop))

    #out.append(linesep)
    #out.append('Command and options: ')
    #full_cmd = ' '.join(sys.argv)
    #out.append('  ' + textwrap.fill(full_cmd, 80))

    out.append(linesep)
    if elapsed_time < 1:
        out.append('Time used: ' + '%.2f s' % elapsed_time)
    else:
        out.append('Time used: ' + '%s m %s s' % chilli.seconds2min_sec(int(elapsed_time)))
    out.append('Date: %s' % chilli.format_show_time())
    out.append(linesep)

    return out

def draw_graphical_alignment_primer(amp, oligos, options, mid_seq):
    '''Draw the graphical alignment for each of the amplicons'''
    pseq = amp['pseq']
    p_ID = amp['pid']
    p_qseq = amp['p_qseq']
    p_aseq = amp['p_aseq']
    p_sseq = amp['p_sseq']
    p_tail = amp['p_tail']
    f_len = amp['plen']
    f_3_pos = amp['f3_pos'] 
    p_qb = len(p_tail) + 1
    p_qe = f_len
    p_sb = f_3_pos - len(p_aseq) + 1

    mseq = amp['mseq']
    m_ID = amp['mid']
    m_qseq = amp['m_qseq']
    m_aseq = amp['m_aseq']
    m_sseq = amp['m_sseq']
    m_tail = amp['m_tail']
    r_len = amp['mlen']
    r_3_pos = amp['r3_pos']
    m_qe = len(m_tail) + 1
    m_qb = r_len
    m_se = r_3_pos + len(m_aseq)

    mid_seq_len = len(mid_seq)
    if mid_seq_len >= 40:
        link_seq = '%s%s...%s%s' % (p_sseq, mid_seq[:20], mid_seq[-20:], m_sseq)
    else:
        half_len = int(mid_seq_len / 2)
        link_seq = '%s%s...%s%s' % (p_sseq, mid_seq[:half_len], mid_seq[-half_len:], m_sseq)
    line0 = ' '*len(p_tail) + '5\' %s 3\'' % link_seq
    line1 = ' '*(3+len(p_tail)) + str(p_sb).ljust(len(link_seq) - len(m_sseq)) + m_aseq
    line2 = '3\' '.rjust(len(line1) - len(m_aseq)) + m_qseq + m_tail + ' 5\''
    if m_qe <= 2:
        line3 = str(m_qb).rjust(len(line1) - len(m_aseq)+1) + str(m_qe).rjust(len(m_qseq)-1) 
    else:
        line3 = str(m_qb).rjust(len(line1) - len(m_aseq)+1) + str(m_qe).rjust(len(m_qseq)-1) + '1'.rjust(len(m_tail))
    line4 = (m_ID + '<<<').rjust(len(line2))
    line_2 = '5\' ' + p_tail + p_qseq + ' 3\''
    line_1 = ' '*(3+len(p_tail)) + p_aseq + ' '*(len(link_seq) - len(p_aseq) - 1) + str(m_se)
    if p_qb <= 2:
        line_3 = ' '*(len(p_tail)+3) + str(p_qb) + ' '*(len(p_qseq) - len(str(p_qb)) - 1) + str(p_qe)
    else:
        line_3 = ' '*3 + '1' + ' '*(len(p_tail) - 1) + str(p_qb) + ' '*(len(p_qseq) - len(str(p_qb)) - 1) + str(p_qe)
    line_4 = '>>>' + p_ID

    lines = [line_4, line_3, line_2, line_1, line0, line1, line2, line3, line4]

    return os.linesep.join(lines)

def get_mid_seq(mid_seq_id_list, options, session_dir, db):
    '''Get the amp sequence using twoBitToFa from Blat suite'''
    mid_seq_list = []
    if len(mid_seq_id_list) == 0:
        return mid_seq_list

    tmp_name = '%s.mid_seq_id_list.txt.tmp' % (os.path.basename(db))
    mid_seq_id_list_file = os.path.join(session_dir, tmp_name)
    fh = open(mid_seq_id_list_file, 'w')
    fh.write(os.linesep.join(mid_seq_id_list))
    fh.close()

    tmp_name = '%s.twoBitToFa_output.txt.tmp' % (os.path.basename(db))
    outfile = os.path.join(session_dir, tmp_name)
    twoBitDB = db + '.2bit'

    cmd = '%s%stwoBitToFa -seqList=%s %s %s' % (bin_path, os.sep, mid_seq_id_list_file, twoBitDB, outfile)

    try:
        out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    except:
        msg = 'Error: twoBitToFa running error.'
	print msg
        print2stderr(msg)
    
    if err:
	print err
        print2stderr(err)

    try:
        fh = open(outfile)
    except:
        msg = 'Error: twoBitToFa running error: no output file produced.'
	print msg
        print2stderr(msg)

    records = FastaFormatParser.parse(fh)
    mid_seq_list = []
    for record in records:
        mid_seq_list.append(record['seq'].lower())

    fh.close()
    return mid_seq_list

def print2stderr(msg):
    '''Print the error message to STDERR'''
    print >> sys.stderr, msg
    exit()

def get_pos_data(data):
    ''''''
    pos_dict = {}
    for record in data.split(';'):
        fields = record.split(':')
        pos_dict[int(fields[0])] = [int(pos) for pos in fields[1].split(',')]

    return pos_dict

def get_position(options, mer_id, db):
    '''Get position of the mer from SQLite3 database'''
    plus_pos = {}
    minus_pos = {}
    dbname = db + '.sqlite3.db'
    conn = sqlite3.connect(dbname)
    cur = conn.cursor()
    query = "select plus, minus from pos where mer_id=%s" %  mer_id
    cur.execute(query)
    try:
	(plus, minus) = cur.fetchone()
    except:
	print "Error found when retrieving position values from indexed database"
	print "Is the k-value right?"
	exit()

    if plus:
	plus_pos = get_pos_data(plus)

    if minus:
	minus_pos = get_pos_data(minus)

    cur.close()
    conn.close()
    return plus_pos, minus_pos

def check_infile(options):
    '''Check and return Oligos'''
    err_or_degenerate = SeqCheck.fasta_format_check(options.infile)
    if err_or_degenerate in ['yes', 'no']:
	global degenerate
	degenerate = err_or_degenerate
    else:
	print2stderr(err_or_degenerate)

    options.infile.seek(0)

    if degenerate == 'no':
        oligos = FastaFormatParser.parse(options.infile)
    else:
        oligos = DegenerateSeqConvetor.convert(FastaFormatParser.parse(options.infile))

    options.infile.close()

    return oligos

def primer_process(options, session_dir, fcdict, db, oligos):
    '''Primer Process'''
    #options.processor = int(options.processor)
    oligo_pos = []
    oligo_id_list = []
    #t1 = time.time()
    for oligo in oligos:
        primer_seq = oligo['seq']
        oligo_id_list.append(oligo['id'])

        #mer = primer_seq[-options.k:]
        mer = primer_seq[-options.k_value:]
        mer_id = chilli.DNA2int(mer)

        # p for plus strand, m for minus strand
        p_pos_list, m_pos_list = get_position(options, mer_id, db)

        oligo_pos.append({
            'p_list' : p_pos_list,
            'm_list' : m_pos_list,
        })

    #t2 = time.time()
    #cost = t2 - t1
    #print cost

    product = []
    binding_range = []
    binding_primer = []
    for i in xrange(len(oligos)):
        p_list = oligo_pos[i]['p_list']
        p_oligo_length = oligos[i]['size']
        for k in xrange(len(oligos)):
            m_list = oligo_pos[k]['m_list']
            m_oligo_length = oligos[k]['size']

            for j in p_list.iterkeys():
                hid = str(j) # Because the database has been re-formated
                try:
                    p_pos = p_list[j]
                    m_pos = m_list[j]
                except:
                    continue

                for p in p_pos:
                    left = get_pos_range(p, m_pos)
                    for pos_index in xrange(left, len(m_pos)):
			f3_pos = p+1
			r3_pos = m_pos[pos_index]

                        # The amplicon size <= p.len + m.len
			if f3_pos >= r3_pos:
			    continue

                        product_size = p_oligo_length + m_pos[pos_index] - p + m_oligo_length - 1

                        if product_size < options.size_start:
                            continue
                        if product_size > options.size_stop:
                            break

                        p_start = p - p_oligo_length + 1
                        if p_start < 0:
                            p_start = 0

                        m_stop = r3_pos + m_oligo_length
                        if m_stop > fcdict[hid]['size']:
                            m_stop = fcdict[hid]['size'] 

                        binding_range.append('%s:%s-%s' % (hid, p_start, p + 1))
                        # Reverse: Correction for reverse
                        binding_range.append('%s:%s-%s' % (hid, r3_pos, m_stop))

                        amp = {
                            'hid' : hid,
                            'pid' : oligos[i]['id'],
                            'mid' : oligos[k]['id'],
                            'plen' : p_oligo_length,
                            'mlen' : m_oligo_length,
                            'pseq' : oligos[i]['seq'],
                            'mseq' : Seq.rev_com(oligos[k]['seq']),
                            'size' : product_size,
                            'f3_pos' : f3_pos,
                            'r3_pos' : r3_pos,
                        }

                        product.append(amp)

    return product, binding_range

def get_pos_range(value, list_com):
    '''Based on Binary Search'''
    left = 0
    right = len(list_com) - 1
    #print list_com
    if right == 1:
        return left
    while (left < right):
        mid = (left + right) // 2
        #print mid, left, right
        if value == list_com[mid]:
            left = right = mid
        elif value > list_com[mid]:
            left = mid + 1
        else: 
            right = mid

    return left

def extract_by_twoBitToFa(id_pos_range_list, options, session_dir, db):
    '''Get the amp sequence using twoBitToFa from Blat suite'''
    seq_list = []
    if len(id_pos_range_list) == 0:
        return seq_list

    tmp_name = '%s.id_pos_range_list.txt.tmp' % (os.path.basename(db))
    id_pos_range_list_file = os.path.join(session_dir, tmp_name)
    fh = open(id_pos_range_list_file, 'w')
    fh.write(os.linesep.join(id_pos_range_list))
    fh.close()

    tmp_name = '%s.twoBitToFa_output.txt.tmp' % (os.path.basename(db))
    outfile = os.path.join(session_dir, tmp_name)
    twoBitDB = db + '.2bit'

    cmd = '%s%stwoBitToFa -seqList=%s %s %s' % (bin_path, os.sep, id_pos_range_list_file, twoBitDB, outfile)

    try:
        out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    except:
        msg = 'Error: twoBitToFa running error.'
        print2stderr(msg)
    
    if err:
        print2stderr(err)

    try:
        fh = open(outfile)
    except:
        msg = 'Error: twoBitToFa running error: no output file produced.'
        print2stderr(msg)

    records = FastaFormatParser.parse(fh)
    seq_list = []
    for record in records:
        seq_list.append(record['seq'].lower())

    fh.close()
    return seq_list

def Thermodynamics_alignment(fp, ts, primer_type):
    '''Alignment'''
    fp = fp.upper()
    ts = ts.upper()
    aseq = [':',] # Alignment seq
    if primer_type == 'forward':
	fp = fp[::-1]
	ts = ts[::-1]

    i = 1
    while i < len(ts):
        fp_di = fp[(i-1) : (i+1)]
        ts_di = ts[(i-1) : (i+1)]
	if fp_di == ts_di:
            aseq.append(':')
	else:
	    if fp_di in nn_mm_data:
	        if ts_di in nn_mm_data[fp_di]:
                    try:
                        if fp[i+1] == ts[i+1]:
                            aseq.append('.:')
			    i += 1
		        else:
		            break
                    except:
                        break
		else:
		    break
	    else:
	        break

	i += 1

    aseq = ''.join(aseq)
    if primer_type == 'forward':
        return aseq[::-1]
    else:
        return aseq

def get_align_seq(seq_list, options, product):
    '''Alignment seq'''
    filter_product = []
    for i in xrange(len(product)):
        amp = product[i]
        pseq = amp['pseq']
        mseq = amp['mseq']
        pts = seq_list[2*i] # Forward primer target sequence
        mts = seq_list[2*i+1] # Reverse primer target sequence
        # ts for target sequence
        #p_aseq = Watson_Click_alignment(pseq, pts, 'forward')
        p_aseq = Thermodynamics_alignment(pseq, pts, 'forward')
        p_aseq_len = len(p_aseq)
        p_qseq = pseq[-p_aseq_len:].upper()
        p_sseq = pts[-p_aseq_len:].upper()
        p_tail = pseq[:-p_aseq_len]

        p_thermo =  TmDeltaG.Cal(p_qseq, Seq.complement(p_sseq), mono_conc=options.mono_conc, diva_conc=options.diva_conc, dntp_conc=options.dntp_conc)
        p_DeltaG = p_thermo.DeltaG
        p_Tm = p_thermo.Tm
        if p_Tm < float(options.tm_start) or p_Tm > float(options.tm_stop):
            continue

        #m_aseq = Watson_Click_alignment(mseq, mts, 'reverse')
        m_aseq = Thermodynamics_alignment(mseq, mts, 'reverse')
        m_aseq_len = len(m_aseq)

        m_qseq = mseq[:m_aseq_len].upper()
        m_sseq = mts[:m_aseq_len].upper()
        r_len = amp['mlen']
        if m_aseq_len == r_len:
            m_tail = ''
        else:
            m_tail = mseq[m_aseq_len:]

        m_thermo =  TmDeltaG.Cal(m_qseq, Seq.complement(m_sseq), mono_conc=options.mono_conc, diva_conc=options.diva_conc, dntp_conc=options.dntp_conc)
        m_DeltaG = m_thermo.DeltaG
        m_Tm = m_thermo.Tm
        if m_Tm < float(options.tm_start) or m_Tm > float(options.tm_stop):
            continue

        amp['p_qseq'] = p_qseq
        amp['p_aseq'] = p_aseq
        amp['p_sseq'] = p_sseq
        amp['p_tail'] = p_tail

        amp['m_qseq'] = m_qseq
        amp['m_aseq'] = m_aseq
        amp['m_sseq'] = m_sseq
        amp['m_tail'] = m_tail

        amp['p_Tm'] = p_Tm
        amp['p_DeltaG'] = p_DeltaG
        amp['m_Tm'] = m_Tm
        amp['m_DeltaG'] = m_DeltaG
        filter_product.append(amp)

    return filter_product

def tab_out(amp_list, oligos, options, start_time, session_dir):
    '''Format output in primer task'''
    # amp_id, fp_id, rp_id, ppc, size, gc, fp_tm, fp_dg, rp_tm, rp_dg, seq, hit_id
    options.outfile.write("AmpID\tFpID\tRpID\tHitID\tPPC\tSize\tAmpGC\tFpTm\tRpTm\tFpDg\tRpDg\tBindingStart\tBindingStop\tAmpSeq\n")
    sn = 0
    amp_list.sort(key=itemgetter(1, 2), reverse=True)
    for ave_Tm, ppc, amp_len, amp in amp_list:
        sn = sn + 1

        amp_len = amp['size']

	# p for plus, m for minus primer # History reason

        p_qid = amp['pid']
        f_3_pos = amp['f3_pos'] 
        p_qseq = amp['p_qseq']
        p_aseq = amp['p_aseq']
        p_tail = amp['p_tail']
        p_Tm = amp['p_Tm']
        p_DeltaG = amp['p_DeltaG']
        p_sb = f_3_pos - len(p_aseq) + 1

        m_qid = amp['mid']
        r_3_pos = amp['r3_pos']
        m_qseq = amp['m_qseq']
        m_aseq = amp['m_aseq']
        m_tail = amp['m_tail']
        m_Tm = amp['m_Tm']
        m_DeltaG = amp['m_DeltaG']
        m_se = r_3_pos + len(m_aseq)

        mid_seq = amp['mid_seq']
        real_hid = amp['real_hid']

        amp_seq = p_tail + p_qseq + mid_seq + m_qseq + m_tail
        amp_GC = chilli.cal_GC_content(amp_seq, + amp_len)

        if p_qid == m_qid:
            ppc = '-%.1f' % ppc
        else:
            ppc = '%.1f' % ppc

	options.outfile.write("%d\t%s\t%s\t%s\t%s\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\t%d\t%s\n" % (sn, p_qid, m_qid, real_hid, ppc, amp_len, amp_GC, p_Tm, m_Tm, p_DeltaG, m_DeltaG, p_sb, m_se, amp_seq))


def process_primer(options, session_dir):
    '''Primer task'''
    if not isinstance(options.database, list):
        options.database = [options.database]

    if not isinstance(options.infile, file):
        options.infile = open(options.infile)

    oligos = check_infile(options)
    amp = []
    for db in options.database:
        fcdict_cache = db + '.uni'
        fcdict = chilli.get_cache(fcdict_cache)
        product, binding_range = primer_process(options, session_dir, fcdict, db, oligos)
        seq_list = extract_by_twoBitToFa(binding_range, options, session_dir, db)
        filter_product = get_align_seq(seq_list, options, product)
        amp_list = primer_analysis(filter_product, options, oligos, session_dir, fcdict, db)
        amp.extend(amp_list)

    return amp, oligos

def main():
    '''Main control function'''
    options = get_opt()

    session_dir = chilli.session()

    start_time = time.time()
    amp_list, oligos = process_primer(options, session_dir)
    if options.tab:
	tab_out(amp_list, oligos, options, start_time, session_dir)
    else:
	format_output_primer(amp_list, oligos, options, start_time, session_dir)


if __name__ == '__main__':
    main()

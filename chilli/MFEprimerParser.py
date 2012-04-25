#!/usr/bin/env python
# -*- coding:utf-8 -*-
from __future__ import division
'''A Parser for MFEprimer output file


by Wubin Qu <quwubin@gmail.com>,
Copyright @ 2010, All Rights Reserved.
'''

Author  = 'Wubin Qu  <quwubin@gmail.com>'
Date    = 'Jul-07-2011 12:42:35'

import sys
import os
import re
import argparse
from pprint import pprint

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

def parse(fh):
    '''Parse main module'''

    amp = []

    # Filter the first description lines
    while True:
        try:
            line = fh.next()
        except:
            break

        if line.startswith('Details'):
            break

    while True:
        try:
            line = fh.next()
        except:
            break

        line = line.strip()
        if not line:
            continue

        # The end
        if line.startswith('#'):
            break

        # The begin
        if re.match('\d+:', line.strip()):
            amp_id = line.split(':')[0].strip()
            hit_id = line.split('==>')[1].strip()
            continue

        if re.match('PPC', line):
            fields = line.split(',')
            ppc = fields[0].split('=')[1].strip('% ')
            size = fields[1].split('=')[1].strip('bp ')
            gc = fields[2].split('=')[1].strip('% ')
            continue

        if re.match('FP', line):
            fields = line.split(',')
            fp_tm = fields[0].split('=')[1].split()[0].strip()
            fp_dg = fields[1].split('=')[1].split()[0].strip()
            continue

        if re.match('RP', line):
            fields = line.split(',')
            rp_tm = fields[0].split('=')[1].split()[0].strip()
            rp_dg = fields[1].split('=')[1].split()[0].strip()
            continue

        if line.startswith('>>>'):
            fp_id = line[3:]
            continue

        if line.strip().endswith('<<<'):
            rp_id = line[:-3]
            continue
        # Handle amplicon sequence in fasta format
        if re.match('>\d+', line):
            seq = []
            while True:
                try:
                    line = fh.next().strip()
                except:
                    break

                if not line:
                    amp.append(AmpCreator(amp_id, fp_id, rp_id, float(ppc), 
                                          int(size), float(gc), float(fp_tm), 
                                          float(fp_dg), float(rp_tm), float(rp_dg), 
                                          ''.join(seq), hit_id))

                    break
                else:
                    seq.append(line)

    return amp

class AmpCreator:
    def __init__(self, amp_id, fp_id, rp_id, ppc, size, gc, fp_tm, fp_dg, rp_tm, rp_dg, seq, hit_id):
        self.id = amp_id
        self.seq = seq
        self.gc = gc
        self.ppc = ppc
        self.size = size
        self.fp_id = fp_id
        self.rp_id = rp_id
        self.fp_tm = fp_tm
        self.rp_tm = rp_tm
        self.fp_dg = fp_dg
        self.rp_dg = rp_dg
        self.hit_id = hit_id

def main ():
    '''Main'''
    get_opt()
    amp_list = parse(args.infile)
    for amp in amp_list:
        print amp.id, amp.fp_id, amp.rp_id
        print amp.size, amp.gc, amp.fp_tm, amp.rp_tm, amp.fp_dg, amp.rp_dg
        print amp.seq

if __name__ == '__main__':
    main()


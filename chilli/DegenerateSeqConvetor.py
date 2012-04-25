#!/usr/bin/env python
# -*- coding:utf-8 -*-
'''Convert degenerate sequence to multiple normal sequences.

Usually used in degenerate primer specificity checking.

by Wubin Qu <quwubin@gmail.com>,
Copyright @ 2010, All Rights Reserved.
'''

Author  = 'Wubin Qu  <quwubin@gmail.com> CZlab, BIRM, China'
Date    = 'Apr-13-2010 09:21:25'
License = 'GPL v3'
Version = '1.0'

import sys
import os
from optparse import OptionParser
import FastaFormatParser as FFP
import chilli

iupac_dict = {
    'a' : 'a',
    'A' : 'A',
    't' : 't',
    'T' : 'T',
    'c' : 'c',
    'C' : 'C',
    'g' : 'g',
    'G' : 'G',
    'R' : ('G', 'A'),
    'r' : ('g', 'a'),
    'Y' : ('T', 'C'),
    'y' : ('t', 'c'),
    'S' : ('G', 'C'),
    's' : ('g', 'c'),
    'W' : ('T', 'A'),
    'w' : ('t', 'a'),
    'K' : ('G', 'T'),
    'k' : ('g', 't'),
    'M' : ('A', 'C'),
    'm' : ('a', 'c'),
    'D' : ('G', 'T', 'A'),
    'd' : ('g', 't', 'a'),
    'H' : ('T', 'A', 'C'),
    'h' : ('t', 'a', 'c'),
    'B' : ('G', 'T', 'C'),
    'b' : ('g', 't', 'c'),
    'V' : ('G', 'A', 'C'),
    'v' : ('g', 'a', 'c'),
    'N' : ('G', 'A', 'T', 'C'),
    'n' : ('g', 'a', 't', 'c'),
    'I' : ('G', 'A', 'T', 'C'),
    'i' : ('g', 'a', 't', 'c'),
}

def get_opt():
    '''Handle options'''
    usage = 'Usage: %prog [options]'

    version = '%prog Version: ' + '%s [%s]' % (Version, Date)
    parser = OptionParser(usage=usage, version=version)
    parser.add_option('-i', '--infile', dest='infile', help='Input file anme. [String]')
    parser.add_option('-o', '--outfile', dest='outfile', default=sys.stdout, help='Output file anme. [String]')
    [options, args] = parser.parse_args()

    if len(args) > 1:
        parser.error('Incorrect argument, add" "-h" for help.')

    if not options.infile:
        parser.error('Input degenerate sequence file needed, add" "-h" for help.')

    return options

def iupac2normal(seq, prefixes=['']):
    '''Returns a list containing degenerate nucleic acid sequences based on the IUPAC table
    Reference: http://lists.open-bio.org/pipermail/biopython/2006-September/003190.html 
    '''
    if len(seq) == 0:
        return prefixes
    else:
        first = seq[0]
        last_seq = seq[1:]
        new_prefixes = []
        for prefix in prefixes:
            if first in iupac_dict:
                for base in iupac_dict[first]:
                    new_prefixes.append(prefix + base)
            else:
                chilli.print2stderr('Error: not recognized base: %s.' % (first))

        return iupac2normal(last_seq, prefixes = new_prefixes)

def convert(records):
    '''Convert'''
    result_fasta_array = []
    for record in records:
        id = record['id']
        desc = record['desc']
        seq = record['seq']
        normal_seq_list = iupac2normal(seq)
        sn = 0
        for normal_seq in normal_seq_list:
            sn += 1
            result_fasta_array.append(
                {
                    'id' : '%s_%s' % (id, sn),
                    'desc' : desc,
                    'seq' : normal_seq,
		    'size' : len(seq),
                }
            )

    return result_fasta_array

def format2fasta(result_fasta_array):
    '''Return the result as the fasta format'''
    out = ''
    for fasta in result_fasta_array:
        id = fasta['id']
        desc = fasta['desc']
        seq = fasta['seq']
        out += '>%s %s%s%s%s' % (id, desc, os.linesep, seq, os.linesep)

    return out

def main ():
    '''Main'''
    options = get_opt()
    try:
        fh = open(options.infile)
    except:
        chilli.print2stderr('Error: can not open file: %s' % (options.infile))

    records = FFP.parse(fh)
    result_fasta_array = convert(records)
    out = format2fasta(result_fasta_array)
    if options.outfile == sys.stdout:
        print out
    else:
        fo = open(options.outfile, 'w')
        fo.write(out)
        fo.close()


if __name__ == '__main__':
    main()


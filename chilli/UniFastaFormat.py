#!/usr/bin/env python
# -*- coding:utf-8 -*-
'''Convert the fasta format database to an standard
format with "num_id original desc line" desc line


by Wubin Qu <quwubin@gmail.com>,
Copyright @ 2010, All Rights Reserved.
'''

Author  = 'Wubin Qu  <quwubin@gmail.com> CZlab, BIRM, China'
Date    = 'Dec-22-2010 13:57:46'
License = 'Please contact Wubin Qu <quwubin.com>'
Version = '1.0'

import sys
import os
from optparse import OptionParser
import chilli
import FastaSimpleFormatParser


def get_opt():
    '''Handle options'''
    usage = 'Usage: %prog [options]'

    version = '%prog Version: ' + '%s [%s]' % (Version, Date)
    parser = OptionParser(usage=usage, version=version)
    parser.add_option('-i', '--infile', dest='infile', help='Input file anme. [String]')
    [options, args] = parser.parse_args()

    if len(args) > 1:
        parser.error('Incorrect argument, add" "-h" for help.')

    if not options.infile:
        parser.error('Input file needed, add" "-h" for help.')

    return options

def print2stderr(msg):
    '''Print msg to sys.stderr and exit the program'''
    print >> sys.stderr, msg
    exit()

def convert(infile, outfile, cache_name):
    '''Convert'''
    fh = open(infile)
    fo = open(outfile, 'w')

    pre_fcdict = {}
    sn = 0
    for line in fh:
        line = line.strip()
        if line.startswith('>'):
            fasta_id, sep, desc = line.partition(' ')
            pre_fcdict[str(sn)] = {
                'id' : fasta_id[1:80],
                'desc' : desc[:240],
            }
            line = '>%s %s' % (sn, line[1:])
            sn += 1

	fo.write(line + os.linesep)

    fh.close()
    fo.close()

    if len(pre_fcdict) < 1:
        print2stderr('No Fasta format sequences in the database')

    records = FastaSimpleFormatParser.parse(open(outfile))
    fcdict = {}
    for record in records:
        id = record['id']
        desc = pre_fcdict[id]['desc']
        size = record['size']
        fcdict[id] = {
            'id' : pre_fcdict[id]['id'],
            'desc' : desc,
            'size' : size,
        }

    chilli.set_cache(fcdict, cache_name)

def main ():
    '''Main'''
    options = get_opt()
    convert_db = options.infile + '.unifasta'
    cache_name = options.infile + '.uni'
    convert(options.infile, convert_db, cache_name)

if __name__ == '__main__':
    main()


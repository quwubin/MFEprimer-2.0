#!/usr/bin/env python
# -*- coding: utf-8 -*- #   
'''A parser for fasta-format sequence file for bioinformatics

by Wubin Qu <quwubin@gmail.com>
Copyright @ 2010, All Rights Reserved.

USAGE:

    import FastaFormatParser
    options, args = get_opt()
    infile = args[0]
    records = FastaFormatParser.parse(open(infile))
    for record in records:
        print record['id']
        print record['desc']
        print record['seq']

ChangeLog

2011-9-19
    
    * Add cStringIO

2010-2-23

    * v1.0 released.


2010-4-2

    * Bug fixed for blank lines in the fasta format file.

'''

Author = 'Wubin Qu <quwubin@gmail.com>, BIRM, China'
Date = '2010-2-23'
License = 'GPL v3'
Version = '1.2'

from optparse import OptionParser

def get_opt():
    '''Options and args'''
    usage = 'Usage: %prog fasta_format_file'
    version = '%prog' + ' %s [%s]' % (Version, Date)
    parser = OptionParser(usage=usage, version=version)
    [options, args] = parser.parse_args()
    if len(args) != 1:
        parser.error('Incorrect argument, add "-h" for help.')

    return options, args

def parse(fh):
    '''The Parser'''
    records = []

    if not isinstance(fh, file):
        import StringIO
        if isinstance(fh, list):
            import os
            fh = os.linesep.join(fh)

        fh = StringIO.StringIO(fh)

    # Remove the comment and blank lines before the first record
    line = ''
    while True:
        try:
            line = fh.next().strip()
        except StopIteration:
            break

        if not line: continue # Blank line

        if line.startswith('>'):
            break

    # The records
    eof = False
    while True:
        # Reach to the end of the file
        if eof:
            break

        if line.startswith('>'):
            id, sep, desc = line[1:].partition(' ')

        seq_list = []
        while True:
            try:
                line = fh.next().strip()
            except StopIteration:
                eof = True
                break

            if line.startswith('>'):
                break

            seq_list.append(line)

        seq = ''.join(seq_list)
        record = {
            'id' : id,
            'desc' : desc,
            'seq' : seq.lower(),
            'size' : len(seq),
        }
        records.append(record)

    return records

def main ():
    '''Main'''
    options, args = get_opt()
    infile = args[0]
    records = parse(open(infile))
    for record in records:
        print record['id']
        print record['desc']
        print record['seq']
        print record['size']

if __name__ == '__main__':
    main()


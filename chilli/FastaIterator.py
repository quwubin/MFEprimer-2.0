#!/usr/bin/env python
# -*- coding: utf-8 -*- #   
'''A Fasta format Parser return Iterator

by Wubin Qu <quwubin@gmail.com>
Copyright @ 2010-2011, All Rights Reserved.

USAGE:

    import FastaIterator
    infile = args[1]
    records = FastaIterator.parse(open(infile))
    for record in records:
        print record.id
        print record.desc
        print record.seq

ChangeLog

2011-11-16
    * Use iterator

2010-2-23

    * v1.0 released.


2010-4-2

    * Bug fixed for blank lines in the fasta format file.

'''

Author = 'Wubin Qu <quwubin@gmail.com>'
Date = '2011-11-16'
License = 'Please contact Wubin Qu <quwubin@gmail.com>'
Version = '1.0'

import sys

def parse(fh):
    '''
    A Fasta-format Parser return Iterator
    '''
    if not isinstance(fh, file):
        import StringIO
        if isinstance(fh, list):
            import os
            fh = os.linesep.join(fh)

        fh = StringIO.StringIO(fh)

    # Remove the comment and blank lines before the first record
    while True:
        line = fh.readline()
        if not line: return # Blank line
	
	line = line.strip()

        if line.startswith('>'):
            break

    while True:
        if not line.startswith('>'):
            raise ValueError("Records in Fasta files should start with '>' character")

        id, sep, desc = line[1:].partition(' ')

        seq_lines = []
        line = fh.readline()
        while True:
            if not line: break

	    line = line.strip()

            if line.startswith('>'):
                break

	    if not line: 
		line = fh.readline()
		continue

            seq_lines.append(line.replace(' ', '').replace("\r", ''))
            line = fh.readline()

        yield Iterator(id, desc, ''.join(seq_lines))

        if not line: return

    assert False, 'Should not reach this line'

class Iterator:
    '''Create the class'''
    def __init__(self, id, desc, seq):
        self.id = id
        self.desc = desc
        self.seq = seq
        self.size = len(seq)

def main ():
    '''Main'''
    infile = sys.argv[1]
    records = parse(open(infile))
    i = 0
    for a in records:
	i += 1
        print a.id, a.desc, a.size, i
        print a.seq


if __name__ == '__main__':
    main()


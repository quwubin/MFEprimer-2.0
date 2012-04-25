#!/usr/bin/env python
# -*- coding: utf-8 -*- #   
'''A module for validating the format of biological sequences, such 
as Fasta format, etc.

Currently, only Fasta-format supported.

by Wubin Qu <quwubin@gmail.com>
Copyright @ 2010, All Rights Reserved.
'''

Author = 'Wubin Qu <quwubin@gmail.com>, BIRM, China'
Date = '2010-3-11'
License = 'Please contact Wubin Qu <quwubin.com>'
Version = '1.0'

import sys, re, os
from DegenerateSeqConvetor import iupac_dict

def print2stderr(msg):
    '''Print to STDERR'''
    print >> sys.stderr, msg
    exit()

def fasta_format_check(fh, err_path=sys.stderr):
    '''The Parser'''
    # Remove the comment and blank lines before the first record
    degenerate = 'no'

    if isinstance(fh, list):
        import StringIO
        fh = StringIO.StringIO(os.linesep.join(fh))

    line = ''
    while True:
        try:
            line = fh.next().strip()
        except StopIteration:
            break

        if not line: continue # Blank line

        if line.startswith('>'):
            break
	else:
	    msg = 'Illegal line: %s%sThere is no description line or the description line in Fasta-format should begin with symbol ">"' 
	    if err_path == sys.stderr:
		print2stderr(msg % (line, '\n'))
	    else:
		return msg % (line, '<br />')

    # The records
    eof = False
    while True:
        # Reach to the end of the file
        if eof:
            break

        if line.startswith('>'):
	    if len(line) <= 1:
		msg = 'Illegal line: %s%sThere should be a unique sequence name after the great symble.'
                if err_path == sys.stderr:
		    print2stderr(msg % (line, '\n'))
                else:
		    return msg % (line, '<br />')


        seq_list = []
        while True:
            try:
                line = fh.next().strip()
            except StopIteration:
                eof = True
                break

            if line.startswith('>'):
                break

	    if re.search('[^iatcgnryswkmdhbvIATCGNRYSWKMDHBV]+', line):
		msg = 'Illegal line: %s%sThe sequence lines should not contain the ambiguous bases.'
		if err_path == sys.stderr:
		    print2stderr(msg % (line, '\n'))
		else:
		    return msg % (line, '<br />')

	    if re.search('[^atcgATCG]+', line):
		degenerate = 'yes'

            seq_list.append(line)

        seq = ''.join(seq_list)
	if len(seq) <= 0:
	    msg = 'Illegal line: %s%sThe sequence lines should after the annotation line.'
	    if err_path == sys.stderr:
		print2stderr(msg % (line, '\n'))
	    else:
		return msg % (line, '<br />')

    # If not good, it would never reach here
    return degenerate


def main ():
    '''Main'''
    infile = sys.argv[1]
    #err = fasta_format_check(open(infile))
    line_list = open(infile).readlines()
    err = fasta_format_check(line_list)
    if not err:
	print 'Good: File %s is FASTA-format.' % infile

if __name__ == '__main__':
    main()


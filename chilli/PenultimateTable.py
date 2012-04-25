#!/usr/bin/env python
# -*- coding:utf-8 -*-
'''ARMS primer additional mismatch selection table


by Wubin Qu <quwubin@gmail.com>,
Copyright @ 2010, All Rights Reserved.
'''

Author  = 'Wubin Qu  <quwubin@gmail.com>, China'
Date    = 'Nov-13-2010 23:25:42'
Version = '1.0'

import sys
import os
from optparse import OptionParser

'''A penultimate mismatch table from Little, S. (1997) ARMS analysis of point mutations.'''
pt = {
    'AA' : {'A' : 'A', 'G' : 'G', 'C' : 'A', 'T' : 'G'},
    'AG' : {'A' : 'C', 'G' : 'T', 'C' : 'A', 'T' : 'G'},
    'AC' : {'A' : 'G', 'G' : 'A', 'C' : 'C', 'T' : 'T'},
    'TT' : {'A' : 'C', 'G' : 'T', 'C' : 'A', 'T' : 'G'},
    'TG' : {'A' : 'G', 'G' : 'A', 'C' : 'T', 'T' : 'C'},
    'TC' : {'A' : 'C', 'G' : 'T', 'C' : 'A', 'T' : 'G'},
    'CC' : {'A' : 'C', 'G' : 'T', 'C' : 'A', 'T' : 'G'},
    'GG' : {'A' : 'A', 'G' : 'G', 'C' : 'A', 'T' : 'G'},

    'GA' : {'A' : 'C', 'G' : 'T', 'C' : 'A', 'T' : 'G'},
    'CA' : {'A' : 'G', 'G' : 'A', 'C' : 'C', 'T' : 'T'},
    'GT' : {'A' : 'G', 'G' : 'A', 'C' : 'T', 'T' : 'C'},
    'CT' : {'A' : 'C', 'G' : 'T', 'C' : 'A', 'T' : 'G'},
}

def get_opt():
    '''Handle options'''
    usage = 'Usage: %prog [options]'

    version = '%prog Version: ' + '%s [%s]' % (Version, Date)
    parser = OptionParser(usage=usage, version=version)
    parser.add_option('-i', '--infile', dest='infile', help='Input file name. [String]')
    parser.add_option('-o', '--outfile', dest='outfile', help='Output file name. [String]')
    [options, args] = parser.parse_args()

    if len(args) > 1:
        parser.error('Incorrect argument, add" "-h" for help.')

    if not options.infile:
        pass

    return options

def main ():
    '''Main'''
    #options = get_opt()
    print pt

if __name__ == '__main__':
    main()


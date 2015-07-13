#!/usr/bin/env python
from __future__ import division

import os
import sys
import datetime
from time import time
from optparse import OptionParser
import sqlite3
import FastaIterator

D2n_dic = dict(A=0, T=3, C=2, G=1, a=0, t=3, c=2, g=1)
n2D_dic = {0:'A', 3:'T', 2:'C', 1:'G', 0:'a', 3:'t', 2:'c', 1:'g'}

def optget():
    '''parse options'''
    parser = OptionParser()
    parser.add_option("-f", "--file", dest = "filename", help = "DNA file in fasta to be indexed")
    parser.add_option("-k", "--k", dest = "k", type='int', help = "K mer , default is 9", default = 9)

    (options, args) = parser.parse_args()

    if not options.filename:
        print_usage()
        exit()	

    return options

def print_usage():
    print '''
%s: Index DB for MFEprimer-2.0

Usage:

    %s -f human.genomic -k 9

Author: Wubin Qu <quwubin@gmail.com>
Last updated: 2012-4-24
    ''' % (os.path.basename(sys.argv[0]), os.path.basename(sys.argv[0]))

def get_memory_percent():
    '''Print Memory information'''
    import os
    try:
	import psutil
    except:
	print '''psutil module needed.

	You can download and install it from here: http://code.google.com/p/psutil/
	'''
	exit()

    p = psutil.Process(os.getpid())
    return p.memory_percent()

def insert_db(conn, mer_count, plus, minus):
    for mer_id in xrange(mer_count):
        conn.execute("insert into pos (mer_id, plus, minus) values (?, ?, ?)", \
                [mer_id, plus[mer_id], minus[mer_id]])

    conn.commit()

def update_db(conn, mer_count, plus, minus):
    for mer_id in xrange(mer_count):
        (plus_data, minus_data) = conn.execute("select plus, minus from pos where mer_id=?", [mer_id]).fetchone()
        if plus_data:
            if plus[mer_id]:
                plus_data += ';%s' % plus[mer_id]
            else:
                pass
        else:
            plus_data = plus[mer_id]

        if minus_data:
            if minus[mer_id]:
                minus_data += ';%s' % minus[mer_id]
            else:
                pass
        else:
            minus_data = minus[mer_id]

        conn.execute("update pos set plus=?, minus=? where mer_id=?", \
                [plus_data, minus_data, mer_id])

    conn.commit()

def baseN(num, b):
    '''convert non-negative decimal integer n to
    equivalent in another base b (2-36)'''
    return ((num == 0) and  '0' ) or ( baseN(num // b, b).lstrip('0') + "0123456789abcdefghijklmnopqrstuvwxyz"[num % b])

def int2DNA(num, k):
    seq = baseN(num, 4)
    return 'A' * (k-len(seq)) + (''.join([n2D_dic[int(base)] for base in seq]))

def DNA2int_2(seq):
    '''convert a sub-sequence/seq to a non-negative integer'''
    plus_mer = 0
    minus_mer = 0
    length = len(seq) - 1 
    for i, letter in enumerate(seq):
        plus_mer += D2n_dic[letter] * 4 ** (length - i)
        minus_mer += (3 - D2n_dic[letter]) * 4 ** i

    return plus_mer, minus_mer

def DNA2int(seq):
    '''convert a sub-sequence/seq to a non-negative integer'''
    plus_mer = 0
    length = len(seq) - 1 
    for i, letter in enumerate(seq):
        plus_mer += D2n_dic[letter] * 4 ** (length - i)

    return plus_mer

def index(filename, k):
    ''''''
    start = time()

    mer_count = 4**k

    dbname = '.'.join(filename.split('.')[:-1]) + '.sqlite3.db'

    conn = sqlite3.connect(dbname)
    cur = conn.cursor()
    cur.executescript('''
    drop table if exists pos;
    create table pos(
    mer_id integer primary key, 
    plus text,
    minus text
    );''')

    plus = ['']*mer_count
    minus = ['']*mer_count

    is_empty = False
    is_db_new = True

    for record in FastaIterator.parse(open(filename)):
        is_empty = False
        print record.id

        fasta_seq = record.seq
	#print 'Time used: ', time() - start

        plus_mer_list = [''] * mer_count
        minus_mer_list = [''] * mer_count

        i_max = len(fasta_seq) - k
        i = 0
        kmer = fasta_seq[:k]
        while i < i_max:
            #print i, len(fasta_seq), i_max
            #print kmer
            try:
                plus_mer_id, minus_mer_id = DNA2int_2(kmer)
            except:
                #print 'Unrecognized base: %s' % fasta_seq[i+k]
                # Skip the unrecognized base, such as 'N'
                i += 1
                kmer = kmer[1:] + fasta_seq[i+k-1]
                continue

            if plus_mer_list[plus_mer_id]:
                plus_mer_list[plus_mer_id] += ',%i' % (i+k-1)
            else:
                plus_mer_list[plus_mer_id] = str(i+k-1)

            if minus_mer_list[minus_mer_id]:
                minus_mer_list[minus_mer_id] += ',%i' % (i)
            else:
                minus_mer_list[minus_mer_id] = str(i)

            i += 1
            kmer = kmer[1:] + fasta_seq[i+k-1]
            if not i % 100000:
                print "%s: %.2f%%, %s" % (record.id, i/i_max*100, str(datetime.timedelta(seconds=(time() - start))))
        else:
            pass

	#print 'Time used: ', time() - start
        for mer_id in xrange(mer_count):
            if plus_mer_list[mer_id]:
                if plus[mer_id]:
                    plus[mer_id] += ';%s:%s' % (record.id, plus_mer_list[mer_id])
                else:
                    plus[mer_id] = '%s:%s' % (record.id, plus_mer_list[mer_id])

            if minus_mer_list[mer_id]:
                if minus[mer_id]:
                    minus[mer_id] += ';%s:%s' % (record.id, minus_mer_list[mer_id])
                else:
                    minus[mer_id] = '%s:%s' % (record.id, minus_mer_list[mer_id])

        memory_percent = get_memory_percent()
        if memory_percent > 50:
            if is_db_new:
                insert_db(conn, mer_count, plus, minus)
                is_db_new = False
            else:
                update_db(conn, mer_count, plus, minus)

            # Empty the container
            plus = ['']*mer_count
            minus = ['']*mer_count
            is_empty = True

            print 'Empty plus and minus due to the memory: %s.' % memory_percent


    if not is_empty:
        if is_db_new:
            insert_db(conn, mer_count, plus, minus)
        else:
            update_db(conn, mer_count, plus, minus)

    print "Time used: %s" % str(datetime.timedelta(seconds=(time() - start)))
    print 'Done.'

def main():
    '''main'''
    options = optget()
    index(options.filename, options.k)

if __name__ == "__main__":
    main()

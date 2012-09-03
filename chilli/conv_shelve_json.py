#!/usr/bin/env python

# By Wubin Qu
#
import sys
import os

def get_cache_json(filename):
    import json
    with open(filename) as fp:
        data = json.load(fp)

    return data

def get_cache_shelve(filename):
    '''Filename should be full path'''
    import os
    import shelve

    d = shelve.open(filename)
    data = d[os.path.basename(filename)]
    d.close()

    return data

def set_cache(data, filename):
    import json
    with open(filename, 'wb') as fp:
        json.dump(data, fp)

def get_cache(filename):
    '''Filename should be full path'''
    try:
	data = get_cache_json(filename)
    except:
	data = get_cache_shelve(filename)

    return data

def main():
    if len(sys.argv) != 2:
	print "Error: conv_shelve_json.py dir"
	exit()

    for root, dirs, files in os.walk(sys.argv[1]):
	for file in files:
	    if file.endswith('.uni'):
		filename = os.path.join(root, file)
		print filename
		try:
		    data = get_cache_json(filename)
		except:
		    data = get_cache_shelve(filename)
		    set_cache(data, filename)


if __name__ == '__main__':
    main()


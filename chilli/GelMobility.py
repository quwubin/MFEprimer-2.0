#!/usr/bin/env python
# -*- coding: utf-8 -*- #   
'''Mobility of the PCR amplicons in agarose gel electrophoresis
was calculated following the semi-log formula reported by Helling
et al [Ref. 1] in 1974.

    Y = a - b * ln (X + k)

  Where Y is the mobility and X is the size of the PCR amplicons.

  To get an reliable prediction of mobility of PCR amplicons, we
determined these parameters (a, b and k) by statistically analyzing
the commercially available DNA markers of known size (100 bp, 250 bp,
500 bp, 750 bp, 1000 bp and 2000 bp, from TaKaRa Inc.) in agarose gel
with different concentrations. The DNA fragment of 100 bp was taken as
a reference band in the agarose gel electrophoresis, and the relative
mobility was determined according to the 100 bp DNA band. When the 
mobility of the reference band reached 5 cm in running the agarose gel,
the electrophoresis was stoped and the mobility of other DNA bands
were measure to determine the relative mobility. Briefly, the 
parameters in the abovementioned formula were determined on the 
0.5% ~ 2% agarose gel electrophoresis with triplicate measurements
in out lab and shown below for producing the virtual electrophoretogram
in MPprimer program.

For 0.5% agarose gel concentration:

    a = 2.7094
    b = 0.2691
    k = 464.4412
    R^2 = 0.9989

For 1.0% agarose gel concentration:

    a = 2.3977
    b = 0.2700
    k = 73.9788
    R^2 = 0.9984

For 1.5% agarose gel concentration:

    a = 2.3221
    b = 0.2634
    k = 48.0873
    R^2 = 0.9977

For 2.0% agarose gel concentration:

    a = 2.1333
    b = 0.2561
    k = 18.5417
    R^2 = 0.9948

Reference

[1] Helling, R. B., Goodman, H. M., & Boyer, H. W. (1974). Analysis of 
endonuclease R-EcoRI fragments of DNA from lambdoid bacteriophages and 
other viruses by agarose-gel electrophoresis. Journal of Virology, 
14(5): 1235-1244.
'''

# Author: Wubin Qu [quwubin@gmail.com], BIRM, China
# Date: 2009-10-12
# License: GPL v3

import sys

def load_gel_para_dict(gel_conc=1.0, formula='Helling'):
    gel_para_dict = {}
    gel_para_dict['Helling'] = {}

    # Parameter for Helling's formula
    gel_para = {}
    gel_para[0.5] = {
	'a' : 2.7094,
	'b' : 0.2691,
	'k' : 464.4412,
    }

    gel_para[1.0] = {
	'a' : 2.3977,
	'b' : 0.2700,
	'k' : 73.9788,
    }

    gel_para[1.5] = {
	'a' : 2.3221,
	'b' : 0.2634,
	'k' : 48.0873,
    }

    gel_para[2.0] = {
	'a' : 2.1333,
	'b' : 0.2561,
	'k' : 18.5417,
    }

    gel_para_dict['Helling'] = gel_para

    err_msg = 'Gel concentration or formula is illegal, Currently, only 0.5, 1, 1.5, 2.0 for concentration and "Helling" for formula are allowed.'

    try: 
        gel_conc = float(gel_conc)
        a = gel_para_dict[formula][gel_conc]['a']
        b = gel_para_dict[formula][gel_conc]['b']
        k = gel_para_dict[formula][gel_conc]['k']
    except:
	print >> sys.stderr, err_msg
	exit()

    return gel_para_dict, a, b, k

def get_size_range(size, gel_conc=1.0, ref_mobility=50, offset=2, formula='Helling'):
    Y = cal_mobility(size, gel_conc)
    offset = int(offset)
    # Set 2 mm as the distance which the bands can be \
            #seperated by naked eyes
    Xmin = cal_size(Y + offset, gel_conc=gel_conc, ref_mobility=ref_mobility, formula='Helling')
    Xmax = cal_size(Y - offset, gel_conc=gel_conc, ref_mobility=ref_mobility, formula='Helling')

    return Xmin, Xmax

def cal_mobility(X, gel_conc=1.0, ref_mobility=50, formula='Helling'):
    '''Cal mobility based on size'''
    import math
    gel_para_dict, a, b, k = load_gel_para_dict(gel_conc=gel_conc, formula=formula)

    X = float(X)
    gel_conc = float(gel_conc)

    # X: size (bp)
    # ref_mobility: the mobility distance of the fastest DNA segment

    if formula == 'Helling':
        Y = a - b * math.log(X + k)
    else:
        pass
        #Y = math.exp(a - b * math.log(X + k))

    # Y: the relative mobility = mobility distance / ref_mobility
    Y = Y * ref_mobility
    # Y: the mobility distance
    return round(Y, 1)

def cal_size(Y, gel_conc=1.0, ref_mobility=50, formula='Helling'):
    '''Predict size based on the relative mobility'''
    import math

    gel_para_dict, a, b, k = load_gel_para_dict(gel_conc=gel_conc, formula=formula)

    # Y: the mobility distance
    Y = Y / ref_mobility
    # ref_mobility: the mobility distance of the fastest DNA segment
    if formula == 'Helling':
        #Y = a - b * math.log(X + k)
        X = math.exp((a - Y) / b) - k
    else:
        pass

    return int(round(X, 0))

def main ():
    '''Test'''
    import sys
    #X = 100
    #gel_conc = sys.argv[1]
    #mobility = cal_mobility(X, gel_conc, ref_mobility=50, formula='Helling')
    #print mobility
    #mobility = cal_mobility(X, gel_conc, ref_mobility=50)
    #print mobility
    min, max = get_size_range(sys.argv[2], gel_conc=1.0, ref_mobility=50, offset=sys.argv[1], formula='Helling')
    print min, max


if __name__ == '__main__':
    main()


#!/usr/bin/env python
# -*- coding: utf-8 -*- #
'''Thermodynamics Parameters for calculating Tm and Gibbs free energy

References

 Matched base-pair:

 [1] SantaLucia JR (1998) "A unified view of polymer, dumbbell
 and oligonucleotide DNA nearest-neighbor thermodynamics", Proc Natl
 Acad Sci 95:1460-65 http://dx.doi.org/10.1073/pnas.95.4.1460.

 Mismatched base-pair:

 [1]. Allawi, H.T. and SantaLucia, J. (1997) Thermodynamics and NMR 
 of internal GT mismatches in DNA, Biochemistry, 36, 10581-10594.
 [2]. Allawi, H.T. and SantaLucia, J. (1998) Nearest-neighbor 
 thermodynamics of internal A center dot C mismatches in 
 DNA: Sequence dependence and pH effects, Biochemistry, 37, 9435-9444.
 [3]. Allawi, H.T. and SantaLucia, J. (1998) Nearest neighbor 
 thermodynamic parameters for internal G center dot A mismatches 
 in DNA, Biochemistry, 37, 2170-2179.
 [4]. Allawi, H.T. and Santalucia, J. (1998) Thermodynamics of 
 internal C center dot T mismatches in DNA, Nucleic Acids 
 Research, 26, 2694-2701.
 
 NN model:

 [1] SantaLucia, J. (1998) A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics, Proceedings of the National Academy of Sci-ences of the United States of America, 95, 1460-1465.
    
 [2] von Ahsen, N., Wittwer, C.T. and Schutz, E. (2001) Oligonucleotide melting tempera-tures under PCR conditions: Nearest-neighbor corrections for Mg2+, deoxynu-cleotide triphosphate, and dimethyl sulfoxide concentrations with comparison to alternative empirical formulas, Clinical Chemistry, 47, 1956-1961.

by Wubin Qu <quwubin@gmail.com>
Copyright @ 2010, All Rights Reserved.
'''

Author = 'Wubin Qu <quwubin@gmail.com>, BIRM, China'
Date = '2010-3-2'
License = 'GPL v3'
Version = '1.0'

dH_full={
        'AATT' : -7.9, 'TTAA' : -7.9, 
        'ATTA' : -7.2, 'TAAT' : -7.2, 
        'CAGT' : -8.5, 'TGAC' : -8.5, 
        'GTCA' : -8.4, 'ACTG' : -8.4, 
        'CTGA' : -7.8, 'AGTC' : -7.8, 
        'GACT' : -8.2, 'TCAG' : -8.2, 
        'CGGC' : -10.6, 'GCCG' : -9.8, 
        'GGCC' : -8.0, 'CCGG' : -8.0, 
        'initCG' : 0.1, 'initGC' : 0.1, 
        'initAT' : 2.3, 'initTA' : 2.3, 
        # Like pair mismatches 
        'AATA' : 1.2, 'ATAA' : 1.2, 
        'CAGA' : -0.9, 'AGAC' : -0.9, 
        'GACA' : -2.9, 'ACAG' : -2.9, 
        'TAAA' : 4.7, 'AAAT' : 4.7, 
        'ACTC' : 0.0, 'CTCA' : 0.0, 
        'CCGC' : -1.5, 'CGCC' : -1.5, 
        'GCCC' : 3.6, 'CCCG' : 3.6, 
        'TCAC' : 6.1, 'CACT' : 6.1, 
        'AGTG' : -3.1, 'GTGA' : -3.1, 
        'CGGG' : -4.9, 'GGGC' : -4.9, 
        'GGCG' : -6.0, 'GCGG' : -6.0, 
        'TGAG' : 1.6, 'GAGT' : 1.6, 
        'ATTT' : -2.7, 'TTTA' : -2.7, 
        'CTGT' : -5.0, 'TGTC' : -5.0, 
        'GTCT' : -2.2, 'TCTG' : -2.2, 
        'TTAT' : 0.2, 'TATT' : 0.2, 
        # G.T mismatches 
        'AGTT' : 1.0, 'TTGA' : 1.0, 
        'ATTG' : -2.5, 'GTTA' : -2.5, 
        'CGGT' : -4.1, 'TGGC' : -4.1, 
        'CTGG' : -2.8, 'GGTC' : -2.8, 
        'GGCT' : 3.3, 'TCGG' : 3.3, 
        'GGTT' : 5.8, 'TTGG' : 5.8, 
        'GTCG' : -4.4, 'GCTG' : -4.4, 
        'GTTG' : 4.1, 'GTTG' : 4.1, 
        'TGAT' : -0.1, 'TAGT' : -0.1, 
        'TGGT' : -1.4, 'TGGT' : -1.4, 
        'TTAG' : -1.3, 'GATT' : -1.3, 
        #G.A mismatches
        'AATG' : -0.6, 'GTAA' : -0.6, 
        'AGTA' : -0.7, 'ATGA' : -0.7, 
        'CAGG' : -0.7, 'GGAC' : -0.7, 
        'CGGA' : -4.0, 'AGGC' : -4.0, 
        'GACG' : -0.6, 'GCAG' : -0.6, 
        'GGCA' : 0.5, 'ACGG' : 0.5, 
        'TAAG' : 0.7, 'GAAT' : 0.7, 
        'TGAA' : 3.0, 'AAGT' : 3.0, 
        #C.T mismatches
        'ACTT' : 0.7, 'TTCA' : 0.7, 
        'ATTC' : -1.2, 'CTTA' : -1.2, 
        'CCGT' : -0.8, 'TGCC' : -0.8, 
        'CTGC' : -1.5, 'CGTC' : -1.5, 
        'GCCT' : 2.3, 'TCCG' : 2.3, 
        'GTCC' : 5.2, 'CCTG' : 5.2, 
        'TCAT' : 1.2, 'TACT' : 1.2, 
        'TTAC' : 1.0, 'CATT' : 1.0, 
        #A.C mismatches
        'AATC' : 2.3, 'CTAA':2.3, 
        'ACTA' : 5.3, 'ATCA':5.3, 
        'CAGC' : 1.9, 'CGAC':1.9, 
        'CCGA' : 0.6, 'AGCC':0.6, 
        'GACC' : 5.2, 'CCAG':5.2, 
        'GCCA' : -0.7, 'ACCG':-0.7, 
        'TAAC' : 3.4, 'CAAT':3.4, 
        'TCAA' : 7.6, 'AACT':7.6
}

    #--------------------#
    # deltaS (cal/K.mol) #
    #--------------------#
dS_full={
        'AATT' : -22.2, 'TTAA':-22.2, 
        'ATTA' : -20.4, 'TAAT':-21.3, 
        'CAGT' : -22.7, 'TGAC':-22.7, 
        'GTCA' : -22.4, 'ACTG':-22.4, 
        'CTGA' : -21.0, 'AGTC':-21.0, 
        'GACT' : -22.2, 'TCAG':-22.2, 
        'CGGC' : -27.2, 'GCCG':-24.4, 
        'GGCC' : -19.9, 'CCGG':-19.9, 
        'initCG' : -2.8, 'initGC':-2.8, 
        'initAT' : 4.1, 'initTA':4.1, 
        'sym' : -1.4, 
        # : Like:pair:mismatches
        'AATA' : 1.7, 'ATAA':1.7, 
        'CAGA' : -4.2, 'AGAC':-4.2, 
        'GACA' : -9.8, 'ACAG':-9.8, 
        'TAAA' : 12.9, 'AAAT':12.9, 
        'ACTC' : -4.4, 'CTCA':-4.4, 
        'CCGC' : -7.2, 'CGCC':-7.2, 
        'GCCC' : 8.9, 'CCCG':8.9, 
        'TCAC' : 16.4, 'CACT':16.4, 
        'AGTG' : -9.5, 'GTGA':-9.5, 
        'CGGG' : -15.3, 'GGGC':-15.3, 
        'GGCG' : -15.8, 'GCGG':-15.8, 
        'TGAG' : 3.6, 'GAGT':3.6, 
        'ATTT' : -10.8, 'TTTA':-10.8, 
        'CTGT' : -15.8, 'TGTC':-15.8, 
        'GTCT' : -8.4, 'TCTG':-8.4, 
        'TTAT' : -1.5, 'TATT':-1.5, 
        # : G.T:mismatches
        'AGTT' : 0.9, 'TTGA':0.9, 
        'ATTG' : -8.3, 'GTTA':-8.3, 
        'CGGT' : -11.7, 'TGGC':-11.7, 
        'CTGG' : -8.0, 'GGTC':-8.0, 
        'GGCT' : 10.4, 'TCGG':10.4, 
        'GGTT' : 16.3, 'TTGG':16.3, 
        'GTCG' : -12.3, 'GCTG':-12.3, 
        'GTTG' : 9.5, 'GTTG':9.5, 
        'TGAT' : -1.7, 'TAGT':-1.7, 
        'TGGT' : -6.2, 'TGGT':-6.2, 
        'TTAG' : -5.3, 'GATT':-5.3, 
        #  G.A mismatches
        'AATG' : -2.3, 'GTAA' : -2.3, 
        'AGTA' : -2.3, 'ATGA' : -2.3, 
        'CAGG' : -2.3, 'GGAC' : -2.3, 
        'CGGA' : -13.2, 'AGGC' : -13.2, 
        'GACG' : -1.0, 'GCAG' : -1.0, 
        'GGCA' : 3.2, 'ACGG' : 3.2, 
        'TAAG' : 0.7, 'GAAT' : 0.7, 
        'TGAA' : 7.4, 'AAGT' : 7.4, 
        # C.T mismatches
        'ACTT' : 0.2, 'TTCA' : 0.2, 
        'ATTC' : -6.2, 'CTTA' : -6.2, 
        'CCGT' : -4.5, 'TGCC' : -4.5, 
        'CTGC' : -6.1, 'CGTC' : -6.1, 
        'GCCT' : 5.4, 'TCCG' : 5.4, 
        'GTCC' : 13.5, 'CCTG' : 13.5, 
        'TCAT' : 0.7, 'TACT' : 0.7, 
        'TTAC' : 0.7, 'CATT' : 0.7, 
        # A.C mismatches
        'AATC' : 4.6, 'CTAA' : 4.6, 
        'ACTA' : 14.6, 'ATCA' : 14.6, 
        'CAGC' : 3.7, 'CGAC' : 3.7, 
        'CCGA' : -0.6, 'AGCC' : -0.6, 
        'GACC' : 14.2, 'CCAG' : 14.2, 
        'GCCA' : -3.8, 'ACCG' : -3.8, 
        'TAAC' : 8.0, 'CAAT' : 8.0, 
        'TCAA' : 20.2, 'AACT' : 20.2
}


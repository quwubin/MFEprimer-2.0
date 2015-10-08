# MFEprimer-2.0 

A fast thermodynamics-based program for checking PCR primer specificity 

## Online server
* http://mfeprimer.com/
* http://biocompute.bmi.ac.cn/CZlab/MFEprimer-2.0/

## Introduction

Evaluating the specificity of PCR primers is an essential step in PCR primer design. The MFEprimer-2.0 server allows users to check primer specificity against genomic DNA and mRNA/cDNA sequence databases quickly and easily. MFEprimer-2.0 uses a k-mer index algorithm to accelerate the search process for primer binding sites and uses thermodynamics to evaluate binding stability between each primer and its DNA template. Several important characteristics such as the sequence, melting temperature and size of each amplicon, either specific or non-specific, are reported on the results page. Based on these characteristics and the user-friendly output, users can readily draw conclusions about the specificity of PCR primers. Analyses for degenerate primers and multiple PCR primers are also supported in MFEprimer-2.0. In addition, the databases supported by MFEprimer-2.0 are comprehensive, and custom databases can also be supported on request. The MFEprimer-2.0 server does not require a login and is freely available at http://biocompute.bmi.ac.cn/CZlab/MFEprimer-2.0. 

## Who need a local version of MFEprimer?

MFEprimer-2.0 has a web server version and a command-line version. [The web server host in our server](http://biocompute.bmi.ac.cn/CZlab/MFEprimer-2.0/) may meet the need of most users. While in the following conditions, users may need a command-line version or a local mirror of MFEprimer-2.0 web server:

  1. Private or custom database;
  2. Limited Internet access;
  3. Use MFEprimer-2.0 as a part of primer design software;
  4. A lot of primers to be evaluated and the task is beyond our server capability. 

## Versions

  * Command-line: [download the latest release](https://github.com/quwubin/MFEprimer/zipball/master), 
usually named "quwubin-MFEprimer-XXXXXXX.zip".
  * Local server (same as http://biocompute.bmi.ac.cn/CZlab/MFEprimer-2.0/): [download the latest release](https://github.com/quwubin/MFEprimerWeb/zipball/master), usually named "quwubin-MFEprimerWeb-XXXXXXX.zip".

## Command-line installation

### System requirement

  * System: Linux or Mac (not test, you may contact me if you want MFEprimer to run on Mac)
  * Python (>= 2.7) or PyPy (http://pypy.org/). I recommend PyPy, because MFEprimer-2.0 is more than 2 times speed up using pypy versus plain python. [Thanks Daniel Struck for this suggestion, here is his GitHub page: https://github.com/dstruck]
  * psutil (>= 3.0.0): download from here (https://github.com/giampaolo/psutil)

### Installation in Linux

  1. `mv quwubin-MFEprimer-XXXXXXX.zip $HOME/local/`   # You can put it anywhere
  2. `cd $HOME/local/`  # Go to the place
  3. `unzip quwubin-MFEprimer-XXXXXXX.zip`  # Unzip the file
  4. `mv quwubin-MFEprimer-XXXXXXX MFEprimer`  # Rename to normal MFEprimer
  5. `cd MFEprimer/test/`  # get to the test directory 
  6. `../IndexDb.sh test.rna`   # Index the database, it will create three files with suffix: .2bit .uni and .sqlite3.db.
  7. `../MFEprimer.py -i p.fa -d test.rna`   # Run MFEprimer and you will get the results if not errors found.
  8. Done. Good Luck.

### Installation in MacOSX

  1. `mv quwubin-MFEprimer-XXXXXXX.zip $HOME/local/`   # You can put it anywhere
  2. `cd $HOME/local/`  # Go to the place
  3. `unzip quwubin-MFEprimer-XXXXXXX.zip`  # Unzip the file
  4. `mv quwubin-MFEprimer-XXXXXXX MFEprimer`  # Rename to normal MFEprimer
  1. Go to the site http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.i386/
  2. Download the binaries for faToTwoBit and twoBitToFa
  3. Replaced the files in bin/32bit/ and bin/64bit/ with the files that you downloaded
  4. changed the permissions with chmod +x
  5. `cd MFEprimer/test/`  # get to the test directory 
  6. `../IndexDb.sh test.rna`   # Index the database, it will create three files with suffix: .2bit .uni and .sqlite3.db.
  7. `../MFEprimer.py -i p.fa -d test.rna`   # Run MFEprimer and you will get the results if not errors found.
  8. Done. Good Luck.

## Create a local mirror of MFEprimer-2.0

### System requirement

  * System: Linux server, large memory and disk space
  * Python (>= 2.7)
  * psutil: download from here (http://code.google.com/p/psutil/)
  * Flask (>= 0.8): http://flask.pocoo.org/
  * Apache with CGI support: other deployment options, please reading [here](http://flask.pocoo.org/docs/deploying/#deployment).

### Installation and test

  1. `mv quwubin-MFEprimerWeb-XXXXXXX.zip $HOME/local/`   # You can put it anywhere
  2. `cd $HOME/local/`  # Go to the place
  3. `unzip quwubin-MFEprimerWeb-XXXXXXX.zip`  # Unzip the file
  4. `mv quwubin-MFEprimerWeb-XXXXXXX MFEprimerWeb`  # Rename to normal MFEprimerWeb or MFEprimer-2.0 as I do
  1. `cd MFEprimerWeb` # get to the directory
  1. `./FirstRun.sh`  # Configuration 
  1. `./MFEprimerWeb.py` # Start the MFEprimer with Flask built-in web server. To stop it, enter "Ctrl-C"
  1. Use your web browser (IE9, Firefox, Google-Chrome, Safari, Opera), and enter "http://127.0.0.1:5000/"
  1. You will see the interface as shown in Fig. 1.
  1. Next, preparing your custom database.

### System configuration
  
  1. The Flask built-in web server is not suitable for production, you should deploy the MFEprimerWeb version with apache or other web servers. To do so, please read [Flask Deployment Options](http://flask.pocoo.org/docs/deploying/#deployment-options).
  1. For "batch mode", a system cron job is needed to check and run the batch jobs. Type `crontab -e` and add the below line to your crontab file:
```
*/1 * * * * /path/to/python /path/to/MFEprimerWeb/MFEprimer_Run_Jobs.py /path/to/MFEprimerWeb/batch_jobs
```
  2. When the MFEprimer web version is running, it will create several temporary files in "session" and "show_img" directories. We have to add two cron jobs to clean the temporary files regularly. Type `crontab -e` again and add the following lines to your crontab file:
```
  */1 * * * * /path/to/MFEprimerWeb/cron/clean_session.sh 
  */1 * * * * /path/to/MFEprimerWeb/cron/clean_img.sh   
```
  3. You may change the clean period depend on the server disk space status.

![Local server first run](https://github.com/quwubin/image/raw/master/MFEprimer/local_server_first_run.png)
>Fig. 1 MFEprimer-2.0 local server


## Preparing the database

  0. It usually needs large memory and disk space when indexing a database. For example, it will need about 12 GB memory and 80 GB disk space when indexing a human genome database with size of 3 GB in a 64bit Linux server. But for the custom database, which usually in small size (MB level), a personal computer with 2 GB memory may work well. Anyway, I recommend users to choose our server (http://biocompute.bmi.ac.cn/CZlab/MFEprimer-2.0/) first when checking the specificity of primers against public databases, such as human, mouse etc. 
  1. Preparing your custom database in FASTA-format and named it like "viruses.genomic" or "human.rna".
The name convention is "species.type".
  2. If your database is "viruses.genomic" and located in "$HOME/db" directory, then type `$HOME/local/MFEprimer/IndexDb.sh $HOME/db/viruses.genomic`. Please be patient because the indexing process may take several minutes or even hours. If you have downloaded the MFEprimerWeb version, the index command should be `$HOME/local/MFEprimerWeb/mfeprimer/IndexDb.sh $HOME/db/viruses.genomic`.
  3. Type `$HOME/local/MFEprimer/MFEprimer.py -i YourPrimer.fasta -d $HOME/db/viruses.genomic` for checking the specificity of primers against the custom (here is viruses.genomic) database. If you downloaded the MFEprimerWeb version, the test command should be `$HOME/local/MFEprimerWeb/mfeprimer/MFEprimer.py -i YourPrimer.fasta -d $HOME/db/viruses.genomic`.
  4. For the MFEprimerWeb version, a soft link should be created to let MFEprimer-2.0 server to find them. Type
```
    cd $HOME/local/MFEprimerWeb/MFEprimerDB/
    ln -s /$HOME/db/viruses.* .
```
  5. In most cases, the default k-value = 9 is sufficient for primer specificity checking. However, if users need more strict evaluation, the k-value can be adjusted to 8, 7 or even 6, however the running time will increase dramatically. So, to keep the server running well, we don't provide this option for choosing the k-value in our server. As an alternative solution, the command-line version of MFEprimer-2.0 can use any k-value easily with an option, for example, `./MFEprimer.py -i primer.fasta -d test.db -k 7 -o test.mfe`. Local server of MFEprimerWeb version can also use the custom k-value (by adjust the config.py parameter). To be noticed that, the user should re-index the database with the k-value (such as 7) first before run the MFEprimer. Here is the example of index command with custom-defined k-value:
```
$HOME/local/MFEprimerWeb/mfeprimer/IndexDb.sh $HOME/db/viruses.genomic 7
```

## More about "index"
   
Unlike MFEprimer 1.x versions, which use BLAST for primer binding sites search, MFEprimer-2.0 uses the k-mer index algorithm to speed up the primer binding sites search process. This is the speed problem I have to solve, while, the other question force me **MUST** to replace the BLAST. It's the "ACCURACY" problem. As we know that, BLAST is a famous program to find the homology sequence from a database by sequence similarity. However, the annealing process of primer and its target sequence is thermodynamics. They bind to each other just because they are stable in thermodynamics, not because they are matched in base pairs. For example, the mismatch "G-G" contributes as much as the Gibbs free energy of -2.2 kcal/mol to the duplex stability _[SantaLucia 2004]_. So the first step we have to do is to find all the possible binding sites with the k-mer index algorithm], and then to evaluate the binding stability using the Nearest-Neighbor model. 

Someone may think that what if let BLAST do the first step (find all the possible binding sites)? Actually, MFEprimer 1.x version uses this strategy. But we found that BLAST may lose many less significant hits. The reasons are 1)BLAST only reports the significant hits in statistically; 2) the BLAST output was controlled by "-E" and other options. Even we carefully handled these options, we still found many less significant hits were missed.

So, what is the "k-mer index algorithm"? And why this algorithm can't miss possible binding sites? Because MFEprimer-2.0 pre-stored all the positions of all k-mers. I will explain in details.

First of all, I have to explain the "k-mer". A “k-mer" is defined as a short DNA sequence with a length of k nucleotides. For example, "AAATTTCCC" is a mer with k=9.

The “index process” is to store all the positions of all k-mers which appears in all the sequences from a FASTA format database. We use the following Python-like seudo-code and Fig. 2 to illustrate the index process. For advanced users, they may look at our [Python code](https://github.com/quwubin/MFEprimer/blob/master/chilli/mfe_index_db.py) for details.

```
mer_pos_hash = {}
k = 9
for seq in fasta_db:
    mer_pos_hash[seq.id] = {}
    for i in range(seq.length - k):
        mer = seq[ i : i+k ]
        pos = i + k
        if mer not in mer_pos_hash[seq.id]:
            mer_pos_hash[seq.id][mer] = []

        mer_pos_hash[seq.id][mer].append( pos)
```

![Index algorithm](https://github.com/quwubin/image/raw/master/MFEprimer/IndexAlgorithm.png)
> Fig. 2 The k-mer index process in MFEprimer-2.0. Here k = 9 and the green lines show the mers.

We store the positions in SQLite3 database. The database schema is very simple. There are three fields:
  1. mer_id [Integer, Primary key]: the mer id. We don’t store the raw mer string into the database. Instead, we convert the mer string into a unique integers and store this unique integer as the mer_id.
  2. plus [text]: seq_id_1:pos_1,pos_2,pos_3,…,pos_n;seq_id_2:pos_1,pos_2…
  3. minus [text]: Same like plus but the sequence is the reverse complement sequence of the plus strand.

According to these explanations, we expect the users know why the k-mer index algorithm is more accurate than BLAST.

## MFEprimer 1.x vs. MFEprimer-2.0 (Running speed comparison)

We did a benchmark test to compare the running speed of MFEprimer 1.x and MFEprimer-2.0.

### Parameters

Here are the machine parameters:

  * System: Linux 2.6.18-194.el5xen #1 SMP x86_64 GNU/Linux
  * CPU: Intel(R) Xeon(R) CPU E5430  @ 2.66GHz  # For all the experiments, only one CPU or one processor was used
  * Memory: 16 GB

Here are the program running parameters:

  * MFEprimer 1.x: -W 9, -e 10000, -B 2000, -b 1000 -s 0.3
  * MFEprimer-2.0: -k 9, --size_stop=2000, --size_start=100, --ppc=30

### Datasets

Here are the data sets with size:

  * _C. elegans_: mRNA (24377 sequences, 37 M), Genome (7 sequences, 97 M bases)
  * Chiken: mRNA (19839 sequences, 41 M bases), Genome (32 sequences, 997 M bases)
  * Human: mRNA (46150 sequences, 125 M bases), Genome (25 sequences, 3134 M bases)

The primers are obtained from [_C. elegans_ RNAi library](http://cmgm.stanford.edu/~kimlab/primers.12-22-99.html). We used the 1 pair, 2 pairs, ... 30 pairs of primers for test.

The running time (in seconds) in Table 1 was counted by Linux command "/usr/bin/time" with option "-f '%e" to count only the elapsed time by the program.

### Results

From Table 1, we can see that for mRNA databases, which have thousands of short sequences (compared with the chromosome sequences), MFEprimer-2.0 have absolutely advantages. Even for 30 pair of primers, MFEprimer-2.0 takes less 10 seconds to finish the work. While for genome databases, which have fewer but large chromosome sequences, the running time increases when the number of primer pair increases. When using 20 pair of primers, it takes about 10 minutes. But for one pair of primers, which is the most normal case, MFEprimer-2.0 can finish the job within 1 seconds. So in the "single mode", MFEprimer-2.0 can quickly (usually less than 10 seconds) return the results for the tasks of one or two pair of primers. But for batch primers with large genome database, it's better to use the "batch mode" to examine the specificity of the PCR primers. 

In general, **the size of database sequence** and the **number of primers** have significant effect on the performance of MFEprimer-2.0. 

### Discussion

Obviously, MFEprimer 1.x has to run BLAST program for searching the whole database every time, however MFEprimer-2.0 only retrieves the position data from the SQL database. So as expected, MFEprimer-2.0 wins in almost every way. It seems that MFEprimer-2.0 lost in the aspects indicated by the yellow region. However, that's because MFEprimer-1.x missed many less significantly hits as I described in the previous section. 


![Benchmark result data](https://github.com/quwubin/image/raw/master/MFEprimer/benchmark_data.png)
> Table 1 Benchmark result data. MFEprimer-2.0 wins in almost every way. It seems that MFEprimer-2.0 lost in the aspects indicated by the yellow region. However, that's because MFEprimer-1.x missed many less significantly hits as I described in the previous section. 

## Bug tracker

Have a bug? Please create an issue here on GitHub!

https://github.com/quwubin/MFEprimer/issues


## Getting help

Email to Wubin Qu (quwubin@gmail.com).

## Citation

>Wubin Qu, Yang Zhou, Yanchun Zhang, Yiming Lu, Xiaolei Wang, Dongsheng Zhao, Yi Yang and Chenggang Zhang\*. MFEprimer-2.0: A fast thermodynamics-based program for checking PCR primer specificity. **_Nucleic Acids Res_**. 2012 (accepted).

>Wubin Qu, Zhiyong Shen, Dongsheng Zhao, Yi Yang and Chenggang Zhang. (2009) 
MFEprimer: multiple factor evaluation of the specificity of PCR primers, 
**_Bioinformatics_**, 25(2), 276-278.

## Authors

**Wubin Qu**

+ http://quwubin.sinaapp.com
+ http://github.com/quwubin

**Chenggang Zhang**

+ zhangcg@bmi.ac.cn

## Copyright and license

Copyright (c) 2008-2012. Wubin Qu (quwubin@gmail.com) and 
Chenggang Zhang (zhangcg@bmi.ac.cn, zcgweb@gmail.com), Beijing Institute of Radiation Medicine.

MFEprimer (all the versions) source and executables are freely available for academic, 
nonprofit and personal use. Commercial licensing information please contact 
Dr. Chenggang Zhang (zhangcg@bmi.ac.cn, zcgweb@gmail.com).

MFEprimer source may be downloaded from "https://github.com/quwubin/MFEprimer".

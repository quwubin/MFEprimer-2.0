MFEprimer
=================

A fast and thermodynamics-based PCR primer specificity checking program

Introduction
-----------

MFEprimer [v2.0] is a program to help the biologist to check the primer specificity 
against the whole genomic DNA and mRNA/cDNA database easily and quickly. It uses the 
k-mer index algorithm to speed up the primer binding sites searching process, and 
uses the thermodynamics to evaluate the binding stability between the primer and DNA template.
Several important characteristics such as the sequence, melting temperature and 
size of each amplicon, either specific or non-specific, were reported in the result page. 
Based on these characteristics and the user-friendly output, users can easily draw a conclusion 
about the specificity of the PCR primers. Degenerate primers and multiple PCR primers analysis 
were also supported in MFEprimer-2.0. In addition, the supported database MFEprimer-2.0 is 
comprehensive, while the custom database is also supported upon user’s request or use the
command-line version of MFEprimer-2.0. The MFEprimer-2.0 server is no login request and 
freely available at: http://biocompute.bmi.ac.cn/CZlab/MFEprimer-2.0/.

To get started, visit http://biocompute.bmi.ac.cn/CZlab/MFEprimer-2.0/!



Quick start
-----------

Web server: http://biocompute.bmi.ac.cn/CZlab/MFEprimer-2.0/

Command-line: [download the latest release](https://github.com/quwubin/MFEprimer/zipball/master), 
usually named "quwubin-MFEprimer-XXXXXXX.zip".


System requirement
-----------

  * Linux, Mac (not test, but it should support)

  * Python (>= 2.7)

  * psutil: download from here (http://code.google.com/p/psutil/)

Installation and test
-----------

  1. `mv quwubin-MFEprimer-XXXXXXX.zip $HOME/local/`   # You can put it anywhere
  2. `cd $HOME/local/`  # Go the place
  3. `unzip quwubin-MFEprimer-XXXXXXX.zip`  # Unzip the file
  4. `mv quwubin-MFEprimer-XXXXXXX MFEprimer`  # Rename to normal MFEprimer
  5. `cd MFEprimer/test/`  # change to the test directory 
  6. `../IndexDb.sh test.rna`   # Index the database, it will create three files with suffix: .2bit .uni and .sqlite3.db.
  7. `../MFEprimer.py -i p.fa -d test.rna`   # Run MFEprimer and you will get the results if not errors found.
  8. Done. Good Luck.

Preparing the database
----------

For transparency and insight into our release cycle, and for striving to maintain backward compatibility, Bootstrap will be maintained under the Semantic Versioning guidelines as much as possible.

Releases will be numbered with the follow format:

`<major>.<minor>.<patch>`

And constructed with the following guidelines:

* Breaking backward compatibility bumps the major (and resets the minor and patch)
* New additions without breaking backward compatibility bumps the minor (and resets the patch)
* Bug fixes and misc changes bumps the patch

For more information on SemVer, please visit http://semver.org/.



Bug tracker
-----------

Have a bug? Please create an issue here on GitHub!

https://github.com/twitter/bootstrap/issues



Twitter account
---------------

Keep up to date on announcements and more by following Bootstrap on Twitter, [@TwBootstrap](http://twitter.com/TwBootstrap).



Blog
----

Read more detailed announcements, discussions, and more on [The Official Twitter Bootstrap Blog](http://blog.getbootstrap.com).



Mailing list
------------

Have a question? Ask on our mailing list!

twitter-bootstrap@googlegroups.com

http://groups.google.com/group/twitter-bootstrap



IRC
---

Server: irc.freenode.net

Channel: ##twitter-bootstrap (the double ## is not a typo)



Developers
----------

We have included a makefile with convenience methods for working with the Bootstrap library.

+ **dependencies**
Our makefile depends on you having recess, uglify.js, and jshint installed. To install, just run the following command in npm:

```
$ npm install recess uglify-js jshint -g
```

+ **build** - `make`
Runs the recess compiler to rebuild the `/less` files and compiles the docs pages. Requires recess and uglify-js. <a href="http://twitter.github.com/bootstrap/less.html#compiling">Read more in our docs &raquo;</a>

+ **test** - `make test`
Runs jshint and qunit tests headlessly in phantom js (used for ci). Depends on having phatomjs installed.

+ **watch** - `make watch`
This is a convenience method for watching just Less files and automatically building them whenever you save. Requires the Watchr gem.



Authors
-------

**Mark Otto**

+ http://twitter.com/mdo
+ http://github.com/markdotto

**Jacob Thornton**

+ http://twitter.com/fat
+ http://github.com/fat



Copyright and license
---------------------

Copyright (c) 2008-2012. Wubin Qu <quwubin@gmail.com> and 
Chenggang Zhang <zhangcg@bmi.ac.cn>, Beijing Institute of Radiation Medicine.

MFEprimer (all the versions) source and executables are freely available for academic, 
nonprofit and personal use. Commercial licensing information please contact 
Dr. Chenggang Zhang <zhangcg@bmi.ac.cn>.

MFEprimer source may be downloaded from "https://github.com/quwubin/MFEprimer".

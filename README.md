# PremDPPI
## About
<font size=4> 
  
PremDPPI can identify the mutations that disrupt protein-protein interactions. It can be used for finding mutations that affect protein function. The 3D structure of a protein-protein complex is required for performing the prediction.
  
</font>


## Source code releases
<font size=4> 
  
You can download [releases](https://github.com/minghuilab/PremDPPI/releases) on github.

</font>

## Installation

#### I. PREREQUISITES

<font size=4>
 
PremDPPI requires the following software and packages.

1. PROVEAN

   This is available at the PROVEAN website.

   http://provean.jcvi.org/index.php/

2. NCBI BLAST 2.4.0

   This is available at the NCBI ftp site.

   ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.4.0/

3. DSSP

   This is available at the DSSP website.

   https://swift.cmbi.umcn.nl/gv/dssp/

4. VMD

   This is available at the VMD website.

   https://www.ks.uiuc.edu/Research/vmd/

6. CHARMM

   This is available at the CHARMM website.

   https://www.charmm.org/

7. NAMD

   This is available at the NAMD website.

   https://www.ks.uiuc.edu/Research/namd/

8. iPot

   This is available at the iPot website.

   https://github.com/gjoni/iPot/

9. Python packages: pandas and rpy2

   To install these packages you can use the following command:
</font>

<font size=4>

	$ conda install -c conda-forge pandas
	$ conda install -c r rpy2

</font> 

<font size=4>

10. R packages: randomForest

</font>

<font size=4>

	$ install.packages('randomForest')

</font> 

#### II. INSTALLATION INSTRUCTIONS

<font size=4>

1. Download and/or install prerequisites described above.

2. Download and unpack the distribution:

</font>

<font size=4>

	$ wget https://github.com/minghuilab/PremDPPI/archive/v1.0.tar.gz
	$ tar -zxvf PremDPPI-1.0.tar.gz

</font> 

<font size=4>

3. Change to the source directory:

</font>

<font size=4>

	$ cd PremDPPI-1.0/version1.0

</font> 

<font size=4>

4. Change the path parameters in PremDPPI.py:

</font>

<font size=4>

	workdir = Your working directory
	pathvmd = path for running VMD software  # /usr/local/bin/vmd
	pathcharmm = path for running CHARMM software  # /usr/local/bin/charmm
	pathnamd2 = path for running NAMD software  # /usr/local/bin/NAMD_2.12_Source/Linux-x86_64-g++/namd2
	pathmkdssp = path for running DSSP software  # /usr/local/bin/mkdssp
	pathpsiblast = path for running PSI-BLAST software  # /usr/local/bin/blast/psiblast
	pathblastdb = path for blastdb  # /usr/local/bin/blastdb/nr
	pathprovean = path for PROVEAN software  # /usr/bin/provean.sh
	pathipot = path for iPot software  # /usr/local/bin/iPot/
	
</font>

&nbsp; &nbsp; The FoldX software needs to be installed in the working directory.

#### III. RUNNING PremDPPI

<font size=4>

	$ python PremDPPI.py -i 2021031508311331364218043

</font> 

## Platform

<font size=4>

PremDPPI is only intended to run on *linux* operating systems.

</font>

## Issues

<font size=4>

You will need to have Python 2 (or 3) and R 3.4.0 (or higher) installed.

</font>

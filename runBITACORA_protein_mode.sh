#!/bin/bash

##########################################################
##                                                      ##
##                       BITACORA                       ##
##                                                      ##
##          Bioinformatics tool to assist the           ##
##      comprehensive annotation of gene families       ##
##                                                      ##
##                                                      ##
## Developed by Joel Vizueta                            ##
## Contact: via github or jvizueta@ub.edu               ##
##                                                      ##
##########################################################

VERSION=1.1

##########################################################
##              EXPORT EXECUTABLES TO PATH              ##
##########################################################

# Perl needs to be installed in the system. 
# HMMER and BLAST executables need to be included in $PATH

export PATH=$PATH:/path/to/blast/bin
export PATH=$PATH:/path/to/hmmer/bin

# PATH to BITACORA Scripts folder. 
SCRIPTDIR=/path/to/Scripts

##########################################################
##                   PREPARE THE DATA                   ##
##########################################################

# Include here the name of your species. i.e. NAME=Dmel
NAME=OrganismName

# Write all files with its full or relative PATH

# Include the fasta file containing the annotated proteins
PROTFILE=/path/to/protein.fasta

# Include the folder containing the FPDB databases (Including a fasta and HMM file named as YOURFPDB_db.fasta and YOURFPDB_db.hmm); Multiple FPDB can be included in the folder to be searched
QUERYDIR=/path/to/query_folder


##########################################################
##                 EDITABLE PARAMETERS                  ##
##########################################################

# Set CLEAN=T if you want to clean the output folder. Intermediate files will not be erased but saved in the Intermediate_files folder. Otherwise, set CLEAN=F to keep all files in the output folder
CLEAN=T

# You can modify the E-value used to filter BLAST and HMMER. Default is 1e-5
EVALUE=1e-5

# Number of threads to be used in blast searches
THREADS=1


##########################################################
##                      HOW TO RUN                      ##
##########################################################

# Once you have included all of the above variables, you can run BITACORA as in:
#$ bash runBITACORA_protein_mode.sh


##########################################################
##                   PIPELINE - CODE                    ##
##########################################################

echo -e "\n#######################  Running BITACORA Protein mode  #######################";
echo "BITACORA protein-mode version $VERSION";
date

# Checking if provided data is ok

if [[ ! -f $SCRIPTDIR/check_data_protein_mode.pl ]] ; then
	echo -e "BITACORA can't find Scripts folder in $SCRIPTDIR. Be sure to add also Scripts at the end of the path as /path/Scripts";
	echo -e "BITACORA died with error\n";
	exit 1;
fi

perl $SCRIPTDIR/check_data_protein_mode.pl $PROTFILE $QUERYDIR 2>BITACORAstd.err

ERRORCHECK="$(grep -c 'ERROR' BITACORAstd.err)"

if [ $ERRORCHECK != 0 ]; then
	cat BITACORAstd.err;
	echo -e "BITACORA died with error\n";
	exit 1;
fi

# Run step 1

perl $SCRIPTDIR/runanalysis_protein_mode.pl $NAME $PROTFILE $QUERYDIR $EVALUE $THREADS 2>>BITACORAstd.err

ERRORCHECK="$(grep -c 'ERROR' BITACORAstd.err)"

if [ $ERRORCHECK != 0 ]; then
	cat BITACORAstd.err;
	echo -e "BITACORA died with error\n";
	exit 1;
fi

ERRORCHECK="$(grep -c 'Segmentation' BITACORAstd.err)"

if [ $ERRORCHECK != 0 ]; then
	cat BITACORAstd.err;
	echo -e "BITACORA died with error\n";
	exit 1;
fi


# Cleaning 

if [ $CLEAN = "T" ]; then
	perl $SCRIPTDIR/runcleaning_protein_mode.pl $NAME $QUERYDIR
	echo -e "Cleaning output folders\n";
fi


rm BITACORAstd.err

echo -e "BITACORA completed without errors :)";
date

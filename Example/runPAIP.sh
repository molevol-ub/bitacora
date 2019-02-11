#!/bin/bash

##########################################################
##                                                      ##
## PAIP: Protein annotation and identification pipeline ##
##                                                      ##
## Developed by Joel Vizueta                            ##
## Contact: via github or jvizueta@ub.edu               ##
##                                                      ##
##########################################################

VERSION=1.0

##########################################################
##    		   EXPORT EXECUTABLES TO PATH 				##
##########################################################

# Perl needs to be installed. 
# HMMER and BLAST executables need to be included in $PATH

export PATH=$PATH:/path/to/blast/bin
export PATH=$PATH:/path/to/hmmer/bin

# PATH to PAIP Scripts folder. 
SCRIPTDIR=../Scripts

##########################################################
##                  PREPARE THE DATA 					##
##########################################################

# Include here the name of your species. i.e. NAME=Dmel
NAME=Dmel

# Write all files with its full or relative PATH

# Include the fasta file containing the genome sequences 
GENOME=Files/Drosophila_melanogaster.BDGP6.dna.chromosome.2R.fa

# Include the GFF3 file. Different formats as GTF can be used, however, check the Manual to adapt PAIP to use them
#GFFFILE=../Files/Cdip_ah2p.gff3
GFFFILE=Files/Drosophila_melanogaster.BDGP6.95.chromosome.2R.reformatted.gff3

# Include the fasta file containing the annotated proteins
PROTFILE=Files/Drosophila_melanogaster.BDGP6.95.chromosome.2R.pep.fasta

# Include the folder containing the QUERY database (Including a fasta and HMM file named as QUERY_db.fasta and QUERY_db.hmm); Multiple query DB can be included in the folder to be searched
QUERYDIR=DB



##########################################################
##         		    EDITABLE PARAMETERS					##
##########################################################

# Set CLEAN=T if you wish to clean the output folder. Intermediate files will not be erased but saved in the Intermediate_files folder. Otherwise, set CLEAN=F to keep all files in the output folder
CLEAN=T

# You can modify the E-value used to filter BLAST and HMMER. Default is 1e-5
EVALUE=1e-5

# Maximum length of an intron used to join putative exons of a gene. Default value is conservative and can also join exons from different genes. See the manual for more details.
MAXINTRON=16000

# Number of threads to be used in blast searches
THREADS=2

##########################################################
##         		       HOW TO RUN						##
##########################################################

# Once you have included all of the above variables, you can run PAIP as in:
#$ bash runPAIP.sh



##########################################################
##			   	        PIPELINE 						##
##########################################################

echo -e "\n#######################  Running PAIP  #######################";
echo "PAIP version $VERSION";
date

# Checking if provided data is ok

if [[ ! -f $SCRIPTDIR/check_data.pl ]] ; then
	echo -e "PAIP can't find Scripts folder in $SCRIPTDIR. Be sure to add also Scripts at the end of the path as /path/Scripts";
	echo -e "PAIP died with error\n";
	exit 1;
fi

perl $SCRIPTDIR/check_data.pl $GFFFILE $GENOME $PROTFILE $QUERYDIR 2>PAIPstd.err

ERRORCHECK="$(grep -c 'ERROR' PAIPstd.err)"

if [ $ERRORCHECK != 0 ]; then
	cat PAIPstd.err;
	echo -e "PAIP died with error\n";
	exit 1;
fi

# Run phase 1

perl $SCRIPTDIR/runanalysis.pl $NAME $PROTFILE $QUERYDIR $GFFFILE $GENOME $EVALUE $THREADS 2>>PAIPstd.err

ERRORCHECK="$(grep -c 'ERROR' PAIPstd.err)"

if [ $ERRORCHECK != 0 ]; then
	cat PAIPstd.err;
	echo -e "PAIP died with error\n";
	exit 1;
fi

# Run phase 2

perl $SCRIPTDIR/runanalysis_2ndround_genomic_withgff.pl $NAME $PROTFILE $QUERYDIR $GENOME $GFFFILE $EVALUE $MAXINTRON $THREADS 2>>PAIPstd.err

ERRORCHECK="$(grep -c 'ERROR' PAIPstd.err)"

if [ $ERRORCHECK != 0 ]; then
	cat PAIPstd.err;
	echo -e "PAIP died with error\n";
	exit 1;
fi


# Cleaning 

if [ $CLEAN = "T" ]; then
	perl $SCRIPTDIR/runcleaning.pl $NAME $QUERYDIR
	echo -e "Cleaning output folders\n";
fi


#rm PAIPstd.err

echo -e "PAIP completed without errors :)";
date

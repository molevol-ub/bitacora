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

VERSION=1.4

##########################################################
##              EXPORT EXECUTABLES TO PATH              ##
##########################################################

# Perl needs to be installed in the system. 
# HMMER and BLAST executables need to be included in $PATH

export PATH=$PATH:/path/to/blast/bin
export PATH=$PATH:/path/to/hmmer/bin

# PATH to BITACORA Scripts folder. 
SCRIPTDIR=/path/to/Scripts

# In case of using GeMoMa (set as default), specify the PATH to jar file. Otherwise, set GEMOMA=F in editable parameters section to use the close-proximity method, which does not require any external software
GEMOMAP=/path/to/GeMoMa.jar

##########################################################
##                   PREPARE THE DATA                   ##
##########################################################

# Include here the name of your species. i.e. NAME=Dmel
NAME=OrganismName

# Write all files with its full or relative PATH

# Include the fasta file containing the genome sequences 
GENOME=/path/to/genome.fasta

# Include the folder containing the FPDB databases (Including a fasta and HMM file named as YOURFPDB_db.fasta and YOURFPDB_db.hmm); Multiple FPDB can be included in the folder to be searched
QUERYDIR=/path/to/query_folder


##########################################################
##                 EDITABLE PARAMETERS                  ##
##########################################################

# Set CLEAN=T if you want to clean the output folder. Intermediate files will not be erased but saved in the Intermediate_files folder. Otherwise, set CLEAN=F to keep all files in the same output folder
CLEAN=T

# You can modify the E-value used to filter BLAST and HMMER. Default is 1e-5
EVALUE=1e-3

# Number of threads to be used in blast searches
THREADS=1

# (Default) GEMOMA=T (with upper case) will use GeMoMa software to predict novel genes from TBLASTN alignments (PATH to jar file need to be specified in GEMOMAP variable) 
# Otherwise, set GEMOMA=F to predict new genes by exon proximity (close-proximity method)
GEMOMA=T

# (Used when GEMOMA=F; close-proximity method) Maximum length of an intron used to join putative exons of a gene. Default value is conservative and can also join exons from different genes (labeled in output files with _Ndom) 
# The provided script in Scripts/Tools/get_intron_size_fromgff.pl can estimate intron length statistics for a specific GFF. See the manual for more details
MAXINTRON=15000

# Set GENOMICBLASTP=T in order to conduct both BLASTP and HMMER to curate novel annotated genes (Note that this option is the most sensitive but greatly depends on the database quality and could result in false positives) 
# Otherwise, BITACORA will only use the protein domain (HMMER) to validate new annotated genes (In this case, the probability of detecting all copies is lower, but it will avoid to identify unrelated genes)
GENOMICBLASTP=F

# An additional validation and filtering of the resulting annotations can be conducted using the option ADDFILTER. 
# If ADDFILTER=T, BITACORA will cluster highly similar sequences (with 98% identity; being isoforms or resulting from putative assembly artifacts), and will discard all annotations with a length lower than the specified in FILTERLENGTH parameter.
ADDFILTER=T
FILTERLENGTH=30

# Alternatively, BITACORA can report all annotated genes, without any clustering of identical copies. Set RETAINNONFILTER=T in this case. 
RETAINNONFILTER=F


##########################################################
##                      HOW TO RUN                      ##
##########################################################

# Once you have included all of the above variables, you can run BITACORA as in:
#$ bash runBITACORA_genome_mode.sh


##########################################################
##                   PIPELINE - CODE                    ##
##########################################################

echo -e "\n#######################  Running BITACORA Genome mode  #######################";
echo "BITACORA genome-mode version $VERSION";
date

# Checking if provided data is ok

if [[ ! -f $SCRIPTDIR/check_data_genome_mode.pl ]] ; then
	echo -e "BITACORA can't find Scripts folder in $SCRIPTDIR. Be sure to add also Scripts at the end of the path as /path/Scripts";
	echo -e "BITACORA died with error\n";
	exit 1;
fi

if [ $GEMOMA == "T" ] ; then
	perl $SCRIPTDIR/check_data_genome_mode.pl $GENOME $QUERYDIR $GEMOMA $GEMOMAP 2>BITACORAstd.err
fi

if [ $GEMOMA != "T" ] ; then
	perl $SCRIPTDIR/check_data_genome_mode.pl $GENOME $QUERYDIR $GEMOMA 2>BITACORAstd.err
fi

ERRORCHECK="$(grep -c 'ERROR' BITACORAstd.err)"

if [ $ERRORCHECK != 0 ]; then
	cat BITACORAstd.err;
	echo -e "BITACORA died with error\n";
	exit 1;
fi


# Run step 2

if [ $GEMOMA == "T" ] ; then
	if [ $GENOMICBLASTP == "T" ] ; then
		perl $SCRIPTDIR/runanalysis_2ndround_v2_genomic_nogff_gemoma.pl $NAME $QUERYDIR $GENOME $EVALUE $MAXINTRON $THREADS $GEMOMAP 2>>BITACORAstd.err 2>BITACORAstd.err
	fi

	if [ $GENOMICBLASTP != "T" ] ; then
		perl $SCRIPTDIR/runanalysis_2ndround_genomic_nogff_gemoma.pl $NAME $QUERYDIR $GENOME $EVALUE $MAXINTRON $THREADS $GEMOMAP 2>>BITACORAstd.err 2>BITACORAstd.err
	fi
fi

if [ $GEMOMA != "T" ] ; then
	if [ $GENOMICBLASTP == "T" ] ; then
		perl $SCRIPTDIR/runanalysis_2ndround_v2_genomic_nogff.pl $NAME $QUERYDIR $GENOME $EVALUE $MAXINTRON $THREADS 2>>BITACORAstd.err 2>BITACORAstd.err
	fi
	
	if [ $GENOMICBLASTP != "T" ] ; then
		perl $SCRIPTDIR/runanalysis_2ndround_genomic_nogff.pl $NAME $QUERYDIR $GENOME $EVALUE $MAXINTRON $THREADS 2>>BITACORAstd.err 2>BITACORAstd.err
	fi
fi

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


# Run additional filtering and clustering
	
if [ $ADDFILTER == "T" ] ; then
	perl $SCRIPTDIR/runfiltering_genome_mode.pl $NAME $QUERYDIR $FILTERLENGTH 2>>BITACORAstd.err 2>BITACORAstd.err
fi

ERRORCHECK="$(grep -c 'ERROR' BITACORAstd.err)"

if [ $ERRORCHECK != 0 ]; then
	cat BITACORAstd.err;
	echo -e "BITACORA died with error\n";
	exit 1;
fi


# Cleaning 

if [ $RETAINNONFILTER == "T" ]; then
	perl $SCRIPTDIR/runcleaning_genome_mode_allcopies.pl $NAME $QUERYDIR
	echo -e "Cleaning output folders\n";
fi

if [ $RETAINNONFILTER != "T" ]; then
	if [ $CLEAN == "T" ]; then
	perl $SCRIPTDIR/runcleaning_genome_mode.pl $NAME $QUERYDIR
	echo -e "Cleaning output folders\n";
	fi
fi


rm BITACORAstd.err

echo -e "BITACORA completed without errors :)";
date

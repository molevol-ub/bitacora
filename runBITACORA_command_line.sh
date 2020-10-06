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

VERSION=1.3

# Default values for editable parameters
CLEAN=T
EVALUE=1e-3
THREADS=1
ALGORITHM=gemoma
MAXINTRON=15000
GENOMICBLASTP=F
ADDFILTER=T
FILTERLENGTH=30

BITMODE=full
NAME=Out

# Add here the path to BLAST, HMMER, BITACORA script folder and GeMoMa jar file to avoid specifying them in the command line argument
BPATH=''
HPATH=''
SCRIPTDIR=''
GEMOMAP=''

function usage {
	echo
	echo "Usage:"
	echo "  $(basename $0)"
	echo "      -m BITACORA runnning mode. Specify 'full', 'genome' or 'protein' (Default=$BITMODE)"
	echo "      -q Folder containing the query database or multiple databases, named as Example1_db.fasta and Example1_db.hmm (Mandatory)"
	echo "      -g Genome fasta file (Mandatory in full and genome mode)"
	echo "      -f GFF file (Mandatory in full mode)"
	echo "      -p Protein fasta file (Mandatory in full and protein mode)"
	echo "      -n Organism Name, i.e. 'Dmel' (Default=$NAME)"
	echo "      -sp PATH to BITACORA Scripts folder (Mandatory)"
	echo "      -gp PATH to GeMoMa jar file (Mandatory in case of using this algorithm to predict novel genes)"
	echo "      -bp PATH to BLAST binaries (Optional if it is already included in PATH)"
	echo "      -hp PATH to HMMER executable (Optional if it is already included in PATH)"
	echo "      -c Clean output files if '-c T'. Specify 'T' or 'F' (Default=$CLEAN)"
	echo "      -e E-value used to filter BLAST and HMMER output (Default=$EVALUE)"
	echo "      -t Number of threads (Default=$THREADS)"
	echo "      -a Algorithm used to predict novel genes. Specify 'gemoma' or 'proximity' (Default=$ALGORITHM)"
	echo "      -i Maximum intron length to join putative exons in the close-proximity algorithm (Default=$MAXINTRON)"
	echo "      -b Specify '-b T' to conduct an additional BLASTP search in addition to HMMER to validate novel genes (Default=$GENOMICBLASTP)"	
	echo "      -r Conduct an additional filtering of the annotations if -r T. Specify 'T' or 'F' (Default=$ADDFILTER)"		
	echo "      -l Minimum length to retain identified genes (Default=$FILTERLENGTH)"											
	echo "      -h Show this help. See BITACORA documentation for further details about each parameter"
	echo ""
	echo "  Example for running BITACORA full mode using GeMoMa algorithm and 4 threads:"
	echo "      ./$(basename $0) -m full -sp /path/to/Scripts -gp /path/to/GeMoMa.jar -n Dmel -g /path/to/genome.fasta -f /path/to/GFF.gff3 -p /path/to/protein.fasta -q /path/to/query_folder -t 4"
	echo ""	
	echo "  Example for running BITACORA genome mode using GeMoMa algorithm:"
	echo "      ./$(basename $0) -m genome -sp /path/to/Scripts -gp /path/to/GeMoMa.jar -n Dmel -g /path/to/genome.fasta -q /path/to/query_folder"
	echo ""	
	echo "  Example for running BITACORA protein mode and 2 threads:"
	echo "      ./$(basename $0) -m protein -sp /path/to/Scripts -n Dmel -p /path/to/protein.fasta -q /path/to/query_folder -t 2"
	echo ""	
	echo ""
	exit 0
}


# Read options

if [ "$#" -lt "1" ]; # at least 1 argument
then
	usage
fi

while [ $# -gt 0 ]; do
	case "$1" in
		-h|-help) usage
				;;
		-m) shift
			BITMODE=$1
			;;				
		-bp) shift
			BPATH=$1
			;;
		-hp) shift
			HPATH=$1
			;;
		-sp) shift
			SCRIPTDIR=$1
			;;
		-gp) shift
			GEMOMAP=$1
			;;
		-n)	shift
			NAME=$1
			;;	
		-g)	shift
			GENOME=$1
			;;	
		-f)	shift
			GFFFILE=$1
			;;	
		-p)	shift
			PROTFILE=$1
			;;	
		-q)	shift
			QUERYDIR=$1
			;;		
		-c)	shift
			CLEAN=$1
			;;							
		-e)	shift
			EVALUE=$1
			;;	
		-t)	shift
			THREADS=$1
			;;	
		-a)	shift
			ALGORITHM=$1
			;;	
		-i)	shift
			MAXINTRON=$1
			;;	
		-b)	shift
			GENOMICBLASTP=$1
			;;	
		-r)	shift
			ADDFILTER=$1
			;;	
		-l)	shift
			FILTERLENGTH=$1
			;;							
		*)	echo 
			echo "ERROR - Invalid option: $1"
			echo
			usage
			;;
	esac
	shift
done


GEMOMA=''
if [ $ALGORITHM == "gemoma" ] ; then
	GEMOMA=T
elif [ $ALGORITHM == "proximity" ] ; then
	GEMOMA=F
else
	echo -e "\nERROR, no recognized algorithm was detected. Please specify 'gemoma' or 'proximity' in -a option\n"
	exit 1;
fi

export PATH=$BPATH:$PATH
export PATH=$HPATH:$PATH

if [ $BITMODE == "full" ] ; then
	##########################################################
	##                   PIPELINE - CODE                    ##
	##########################################################

	echo -e "\n#######################  Running BITACORA  #######################";
	echo "BITACORA version $VERSION";
	date

	# Checking if provided data is ok

	if [[ ! -f $SCRIPTDIR/check_data.pl ]] ; then
		echo -e "BITACORA can't find Scripts folder in $SCRIPTDIR. Be sure to add also Scripts at the end of the path as /path/Scripts";
		echo -e "BITACORA died with error\n";
		usage
		exit 1;
	fi

	if [ $GEMOMA == "T" ] ; then
		perl $SCRIPTDIR/check_data.pl $GFFFILE $GENOME $PROTFILE $QUERYDIR $GEMOMA $GEMOMAP 2>BITACORAstd.err
	fi

	if [ $GEMOMA != "T" ] ; then
		perl $SCRIPTDIR/check_data.pl $GFFFILE $GENOME $PROTFILE $QUERYDIR $GEMOMA 2>BITACORAstd.err
	fi

	ERRORCHECK="$(grep -c 'ERROR' BITACORAstd.err)"

	if [ $ERRORCHECK != 0 ]; then
		cat BITACORAstd.err;
		echo -e "BITACORA died with error\n";
		exit 1;
	fi


	# Run step 1

	perl $SCRIPTDIR/runanalysis.pl $NAME $PROTFILE $QUERYDIR $GFFFILE $GENOME $EVALUE $THREADS 2>>BITACORAstd.err

	ERRORCHECK="$(grep -c 'ERROR' BITACORAstd.err)"

	if [ $ERRORCHECK != 0 ]; then
		cat BITACORAstd.err;
		echo -e "BITACORA died with error\n";
		exit 1;
	fi


	# Run step 2

	if [ $GEMOMA == "T" ] ; then
		if [ $GENOMICBLASTP == "T" ] ; then
			perl $SCRIPTDIR/runanalysis_2ndround_v2_genomic_withgff_gemoma.pl $NAME $PROTFILE $QUERYDIR $GENOME $GFFFILE $EVALUE $MAXINTRON $THREADS $GEMOMAP 2>>BITACORAstd.err 2>BITACORAstd.err
		fi

		if [ $GENOMICBLASTP != "T" ] ; then
			perl $SCRIPTDIR/runanalysis_2ndround_genomic_withgff_gemoma.pl $NAME $PROTFILE $QUERYDIR $GENOME $GFFFILE $EVALUE $MAXINTRON $THREADS $GEMOMAP 2>>BITACORAstd.err 2>BITACORAstd.err
		fi
	fi

	if [ $GEMOMA != "T" ] ; then
		if [ $GENOMICBLASTP == "T" ] ; then
			perl $SCRIPTDIR/runanalysis_2ndround_v2_genomic_withgff.pl $NAME $PROTFILE $QUERYDIR $GENOME $GFFFILE $EVALUE $MAXINTRON $THREADS 2>>BITACORAstd.err 2>BITACORAstd.err
		fi

		if [ $GENOMICBLASTP != "T" ] ; then
			perl $SCRIPTDIR/runanalysis_2ndround_genomic_withgff.pl $NAME $PROTFILE $QUERYDIR $GENOME $GFFFILE $EVALUE $MAXINTRON $THREADS 2>>BITACORAstd.err 2>BITACORAstd.err
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
		perl $SCRIPTDIR/runfiltering.pl $NAME $QUERYDIR $FILTERLENGTH 2>>BITACORAstd.err 2>BITACORAstd.err
	fi

	ERRORCHECK="$(grep -c 'ERROR' BITACORAstd.err)"

	if [ $ERRORCHECK != 0 ]; then
		cat BITACORAstd.err;
		echo -e "BITACORA died with error\n";
		exit 1;
	fi


	# Cleaning 

	if [ $CLEAN = "T" ]; then
		perl $SCRIPTDIR/runcleaning.pl $NAME $QUERYDIR
		echo -e "Cleaning output folders\n";
	fi


	rm BITACORAstd.err

	echo -e "BITACORA completed without errors :)";
	date

	exit

elif [ $BITMODE == "genome" ] ; then
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
		usage
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

	if [ $CLEAN = "T" ]; then
		perl $SCRIPTDIR/runcleaning_genome_mode.pl $NAME $QUERYDIR
		echo -e "Cleaning output folders\n";
	fi


	rm BITACORAstd.err

	echo -e "BITACORA completed without errors :)";
	date

	exit

elif [ $BITMODE == "protein" ] ; then
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
		usage
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


	# Run additional filtering and clustering
		
	if [ $ADDFILTER == "T" ] ; then
		perl $SCRIPTDIR/runfiltering_protein_mode.pl $NAME $QUERYDIR $FILTERLENGTH 2>>BITACORAstd.err 2>BITACORAstd.err
	fi

	ERRORCHECK="$(grep -c 'ERROR' BITACORAstd.err)"

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

	exit

else
	echo -e "\nERROR, please specify a valid BITACORA running mode: 'full', 'genome' or 'protein' \n";
	exit
fi

exit


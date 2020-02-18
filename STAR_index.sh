#!/bin/bash
# This script start the pipeline from .fastq.gz files assuming 
# all single-end.fastq.gz files are located in fastqs/ relative to $WORKDIR and
# all paired-end.fastq.gz files are located in paired_fastqs/ relative to $WORKDIR.


########## Step 0: check command line args and make sure files exist ##########
## Initialize variables with default values:
WORKDIR=""
GENOME=""
STAR=""
fastqc=""
skip_qc=false
Genomeversion=""
n_cpus=8
## Parse command line args
while [[ $# -gt 0 ]]; do
	key="$1"
	case $key in
		-S|--STAR)
		STAR="$2"
		shift #pass argument
		;;
		-f|--fasqc)
		fastqc="$2"
		shift #pass argument
		;;
		-g|--genome)
		GENOME="$2"
		shift # past argument
		;;
		-w|--workdir)
		WORKDIR="$2"
		shift # past argument
		;;
		--hg)
		Genomeversion="$2"
		shift # past argument
		;;
		-t|--cpus)
		n_cpus="$2"
		shift # past argument
		;;        
		--skip-qc)
		skip_qc=true
		shift # past argument
		;;		
		-h|--help)
		echo "Usage: ./STAR_INDEX.sh -g <GENOME> -w <WORKDIR>"
		echo "Options:"

		echo "  -S|--STAR: Path for STAR executable"
		echo "  -h|--hg: genome version (19/38)"
		echo "  -t, --cpus: Number of CPUs to use, default to 8"
		exit
		;;
		*)
		# unknown option
		echo "Unknown option: $key, exiting."
		echo "Usage: ./STAR_INDEX -g <GENOME> -w <WORKDIR>"
		echo "Options:"
		echo "  --skip-qc: Skip fastQC steps for the fastq files"
		echo "  -h|--hg: genome version (19/38)"
		echo "  -t, --cpus: Number of CPUs to use, default to 8"
		exit
		;;
	esac
	shift # past argument or value  
done


## Detect number of CPUs and use min(N_CPUS, n_cpus) for jobs
N_CPUS=$(nproc)
N_CPUS=$(($N_CPUS>n_cpus?n_cpus:$N_CPUS))
echo "Number of CPUs to be used: $N_CPUS"


## Check genome versiom
if  [ "$Genomeversion" != 19 ] && [ "$Genomeversion" != 38 ]  ; then
	echo "invalid genome version $hg. Please enter 19 or 38"
	exit 1
else
	echo "hg=hg$Genomeversion"
fi


## Check $WORKDIR
if [[ ! -d $WORKDIR ]]; then
	echo "Could not find working directory: $WORKDIR, exiting. Please make sure the working directory exists"
	exit 1
else
	## Convert to absolute paths
	GENOME=$(readlink -e $GENOME)
	WORKDIR=$(readlink -e $WORKDIR)
	echo "GENOME=$GENOME"
	echo "WORKDIR=$WORKDIR"
fi

if [[ ! -d $STAR ]]; then
	echo "Could not find STAR executalbe: $STAR, exiting. Please make sure the working directory exists"
	exit 1

## Check $GENOME
if [[ ! -d $GENOME ]]; then
	echo "Could not find reference genome: $GENOME, exiting. Please make sure the working directory exists"
	exit 1
else
	GENOME_GTF="$GENOME/hg$Genomeversion.gtf"
	GENOME_FA="$GENOME/hg$Genomeversion.fa"
	if [[ ! -f $GENOME_GTF ]]; then
		echo "$GENOME_GTF not found, exiting"
		exit 1
	fi
	if [[ ! -f $GENOME_FA ]]; then
		echo "$GENOME_FA not found, exiting"
		exit 1
	fi
	STAR_INDEX="$GENOME/"
fi

## Make STAR index if not exists
echo $STAR_INDEX
if [ ! -d $STAR_INDEX ]; then
	echo "STAR index does not exist, building STAR index"
	$STAR\
		--runThreadN $N_CPUS \
		--runMode genomeGenerate \
		--genomeDir $STAR_INDEX \
		--genomeFastaFiles $GENOME_FA \
		--sjdbGTFfile $GENOME_GTF \
		--sjdbOverhang 100
fi


#!/bin/bash
# This script start the pipeline from .fastq.gz files assuming 
# all single-end.fastq.gz files are located in fastqs/ relative to $WORKDIR and
# all paired-end.fastq.gz files are located in paired_fastqs/ relative to $WORKDIR.


########## Step 0: check command line args and make sure files exist ##########
## Initialize variables with default values:
WORKDIR=""
GENOME=""

STAR="${PWD}/STAR"

echo $STAR
fastqc="${PWD}/FastQC/fastqc"
skip_qc=false
Genomeversion=""
n_cpus=8
featureCounts="${PWD}/featureCounts"
echo $featureCounts

## Parse command line args
while [[ $# -gt 0 ]]; do
	key="$1"
	case $key in
		#	-S|--STAR)
		#STAR="$2"
		#shift #pass argument
		#;;
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

	#	echo "  -S|--STAR: Path for STAR executable"
		echo "  --skip-qc: Skip fastQC steps for the fastq files"
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

## Check $fastqc
if (  [ ! -d $fastqc ]  && "$skip_qc" = false ) ; then
	echo "Could not find fastqc executalbe: $fastqc, exiting. Please make sure the working directory exists"
	exit 1

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
fi
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


cd $WORKDIR

## create dirs if not exists
mkdir -p fastQC_output
mkdir -p star_output
mkdir -p featureCount_output


########## Step 2: QC, align and assemble sequencing reads ##########
echo "Started to align reads to the genome and assemble transcriptome"
if [[ -d fastqs ]]; then
	## Align and assemble single-end sequencing reads
	cd fastqs
	for fq in $(ls); do
		basename=$(echo $fq | cut -f1 -d '.')
		if [ "$skip_qc" = false ]; then
			echo "Performing FastQC for $basename"
			$fastqc $fq -o ../fastQC_output
		fi

		echo "Aligning reads from $basename to the reference genome"
		$STAR \
			--genomeDir $STAR_INDEX \
			--sjdbGTFfile $GENOME_GTF \
			--runThreadN $N_CPUS \
			--outSAMstrandField intronMotif \
			--outFilterIntronMotifs RemoveNoncanonical \
			--outFileNamePrefix ../star_output/$basename \
			--readFilesIn $fq \
			--readFilesCommand zcat \
			--outSAMtype BAM Unsorted \
			--outReadsUnmapped Fastx \
			--outSAMmode Full

		suffix="Aligned.out.bam"
		outname="$basename.count.txt"
		bam="../star_output/$basename$suffix"
		$featureCounts \
			-T $N_CPUS \
			-t exon \
			-g gene_id \
			-a $GENOME_GTF \
			-o ../featureCount_output/$outname \
			$bam
	done
	cd ..
fi	

if [[ -d paired_fastqs ]]; then
	## Align and assemble paired-end sequencing reads
	cd paired_fastqs
	for basename in $(ls | cut -f1 -d '_' | sort | uniq); do
		echo $basename
		fq1="_1.fastq"
		fq2="_2.fastq"
		fq1=$basename$fq1
		fq2=$basename$fq2
		if [ "$skip_qc" = false ]; then
			echo "Performing FastQC for $basename"
			$fastqc $fq1 -o ../fastQC_output
			$fastqc $fq2 -o ../fastQC_output
		fi
		echo "Aligning reads from $basename to the reference genome"
		$STAR \
			--genomeDir $STAR_INDEX \
			--sjdbGTFfile $GENOME_GTF \
			--runThreadN $N_CPUS \
			--outSAMstrandField intronMotif \
			--outFilterIntronMotifs RemoveNoncanonical \
			--outFileNamePrefix ../star_output/$basename \
			--readFilesIn $fq1 $fq2 \
			--readFilesCommand zcat \
			--outSAMtype BAM Unsorted \
			--outReadsUnmapped Fastx \
			--outSAMmode Full

		suffix="Aligned.out.bam"
		outname="$basename.count.txt"
		bam="../star_output/$basename$suffix"
		$featureCounts \
			-T $N_CPUS \
			-t exon \
			-g gene_id \
			-a $GENOME_GTF \
			-o ../featureCount_output/$outname \
			$bam
	done
	cd ..
fi

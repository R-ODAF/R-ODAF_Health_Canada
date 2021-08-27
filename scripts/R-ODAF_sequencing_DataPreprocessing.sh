#!/bin/bash
####################################################################################
### Sequencing R-ODAF, Omics Data Analysis Framework for Regulatory application  ###
####################################################################################

####################################################
#### Settings which need to be adapted by user #####
####################################################
source ${HOME}/miniconda3/etc/profile.d/conda.sh
project="2021_Robaire_OPE_TempOSeq" # $1 # Call script using project name as the only argument
# Specify the directory for the output
OUTPUT_DIR="${HOME}/shared/projects/${project}/data/output/"
# Specify location of input fastq files. ALL FILES IN THE FOLDER WILL BE PROCESSED 
RAW_SAMPLE_DIR="${HOME}/shared/projects/${project}/data/raw/"
# Specify extention of input files (".fastq" or ".fastq.gz") 
SUFFIX_INPUTFILES='.fastq.gz'
# Specify the sequencing type (RNASeq or TempOSeq)
SEQTYPE="TempOSeq"
# Specify the sequencing mode used to obtain the data
SEQMODE="single" #specify "paired" or "single" end mode
# Specify the read suffix (e.g. "_R1_001")
PAIRED_END_SUFFIX_FORWARD="_R1_001"
# *IF* paired end mode was used, specify the reverse suffix as well (e.g. "_R2")
PAIRED_END_SUFFIX_REVERSE="_R2_001"

# Choose the main organism for genome alignment (e.g "Rat_6.0.97"). {NOTE: This ID is a label specific for this script and is made for the user to identify which genome version was used. It can contain any text}.
# hg38 | Rat_6.0.84 | S1500 | HumanWT
ORGANISM_GENOME_ID="S1500"
# PATH/Directory in which the genome files are located
# ${HOME}/shared/dbs/human/hg38/ | ${HOME}/shared/dbs/rat/ensembl/rnor6_0/v84/genome
# ${HOME}/shared/dbs/biospyder/TSQR_Scripts_Human_Surrogate_2.0/reference/humansurrogate2_0
# ${HOME}/shared/dbs/biospyder/R-Scripts/Human_Whole_Transcriptome/TSQR_Scripts_HumanWT_1.1/reference/humanwt1_1
GENOME_FILES_DIR="${HOME}/shared/dbs/biospyder/TSQR_Scripts_Human_Surrogate_2.0/reference/humansurrogate2_0"
# Filename of genome fasta file (without path)
# Homo_sapiens_assembly38.fasta | Rnor_6.0.fa 
# S1500: humansurrogate2_0.fa
# Human WT: humanwt1_1.fa
GENOME_FILE_NAME="humansurrogate2_0.fa"
# Filename of GTF file (without path)
# hg38.ensGene.gtf | Rattus_norvegicus.Rnor_6.0.84.andERCC.gtf
# S1500: humansurrogate2_0.gtf
# Human WT: humanwt1_1.gtf
GTF_FILE_NAME="humansurrogate2_0.gtf"
# Whether the genome indexing has already been done. When "Yes" is specified, the indexing will be skipped. If "No" The index will be made
GENOME_INDEX_ALREADY_AVAILABLE="Yes" #Specify "Yes" or "No"
RSEM_INDEX_ALREADY_AVAILABLE="Yes" #Specify "Yes" or "No"
#Specify whether you are working with a large (=human) genome. Specify "Yes" when working with human or "No"
LARGE_GENOME="Yes"

# System parameters
# Specify amount of CPUs to use for alignment step (use 20 or 30)
CPU_FOR_ALIGNMENT=40 
# Specify amount of CPUs to use (use 6 or higher)
CPU_FOR_OTHER=40

### No other input required ###

#################################
# Running the sequencing R-ODAF #
#################################

declare BASEDIR=${OUTPUT_DIR}   #specify working directory
#declare SEQMODE="single"   #specify "paired" or "single" end mode
#declare RAW_SAMPLE_DIR=${BASEDIR}   #specify location of fastq files
if [[ "${SUFFIX_INPUTFILES}" =~ ^\. ]]
 then echo "Suffix has a leading period, continuing."
 else 
 echo "Suffix did not have a leading period, adding one now..."
 SUFFIX_INPUTFILES=${SUFFIX_INPUTFILES/#/.}
 echo ${SUFFIX_INPUTFILES}
fi
declare SUFFIX_IN=${SUFFIX_INPUTFILES} # specify extension of input files (".fastq", ".fastq.gz", etc.)

#For paired end mode: specify the SUFFIXES used for read1 and read2 
declare PAIR1=${PAIRED_END_SUFFIX_FORWARD}
declare PAIR2=${PAIRED_END_SUFFIX_REVERSE}

#Specify reference genome files
declare GENOMEDIR=${GENOME_FILES_DIR}       #location of genome files
declare GENOME="${GENOMEDIR}/${GENOME_FILE_NAME}"  #genome fasta file
declare	GTF="${GENOMEDIR}/${GTF_FILE_NAME}"   #genome GTF file
declare GenomeID=${ORGANISM_GENOME_ID}  #Specify the genome name (e.g. Species+GenomeVersion or Species+DownloadDate) to prevent overwriting other indexed genomes

declare LargeGenome="No"     #Specify "Yes" or "No" to indicate if Human=Large Genome is used
#Specify "Yes" or "No" to indicate if GenomeIndexing has already been done for STAR (if no is specified, the index will be made)
GenomeIndexDone=${GENOME_INDEX_ALREADY_AVAILABLE}
#Specify "Yes" or "No" to indicate if GenomeIndexing has already been done for RSEM (if no is specified, the index will be made)
GenomeIndexRSEMDone=${RSEM_INDEX_ALREADY_AVAILABLE}

declare CPUs=${CPU_FOR_OTHER}   #Specify amount of CPUs to use
declare CPUs_align=${CPU_FOR_ALIGNMENT}   #Specify amount of CPUs to use for alignment step

#######################################################################
### Defining parameters for script to run (no user input necessary) ### 
#######################################################################

declare SOURCEDIR=${RAW_SAMPLE_DIR}
declare TRIMM_DIR="${BASEDIR}/Trimmed_reads/"
declare OUTPUTDIR=${TRIMM_DIR}
declare QC_DIR_fastp="${OUTPUTDIR}/fastpQCoutput/"
declare QC_DIR_multiQC="${OUTPUTDIR}/MultiQC/"
declare align_DIR="${OUTPUTDIR}/STAR/"
declare Quant_DIR="${OUTPUTDIR}/RSEM/"
declare RSEM_GENOMEDIR="${GENOMEDIR}/RSEM/"
declare TEMPOSEQR="${HOME}/shared/projects/${project}/scripts/temposeq/TempO-SeqR_v3.1.R"
declare TEMPOSEQFIGS="${HOME}/shared/projects/${project}/scripts/temposeq/generate_figures.R"

declare SUFFIX1=${SUFFIX_IN}
declare SUFFIX_out="_trimmed${SUFFIX_IN}"

echo "Intializing directories..."
mkdir -p ${TRIMM_DIR}
mkdir -p ${QC_DIR_fastp}
mkdir -p ${QC_DIR_multiQC}
mkdir -p ${align_DIR}
mkdir -p ${Quant_DIR}
mkdir -p ${RSEM_GENOMEDIR}
echo "Done."

echo "Starting logging..."
exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
mydate="$(date +'%d.%m.%Y.%H-%M')"
exec 1>${OUTPUT_DIR}/log_${mydate}.out 2>&1

echo $SHELL

echo "Activating required software."
if [ ${SEQTYPE} == "TempOSeq" ]
then
  conda activate odaf
  #conda activate odaf-star2.7.1 # If genome index is a different version...
elif [ ${SEQTYPE} == "RNASeq" ]
then
  conda activate odaf
else
  echo "Sequencing type not recognized. Quitting."
  break
fi

##################################
### Trimming raw reads : Fastp ###
##################################

# Trimming single end reads 

if [ ${SEQMODE} == "single" ]; then
	declare FILES1="${SOURCEDIR}*${SUFFIX1}"
for FILENAME in ${FILES1[@]}; do
	echo -e "[TRIMMING] fastp: [${FILENAME:${#SOURCEDIR}:-${#SUFFIX1}}]" 
	#Single end
	#Prevent overwrite:
	if [ -e ${TRIMM_DIR}${FILENAME:${#SOURCEDIR}:-${#SUFFIX1}}${SUFFIX_out}  ]; then
		echo "File exists, continuing"
	else
		fastp --in1 ${FILENAME} \
		--out1 ${TRIMM_DIR}${FILENAME:${#SOURCEDIR}:-${#SUFFIX1}}${SUFFIX_out} \
		--json "${QC_DIR_fastp}${FILENAME:${#SOURCEDIR}:-${#SUFFIX1}}_fastp.json" \
		--html ${QC_DIR_fastp}"${FILENAME:${#SOURCEDIR}:-${#SUFFIX1}}fastp.html" \
		--cut_front --cut_front_window_size 1 \
		--cut_front_mean_quality 3 \
		--cut_tail \
		--cut_tail_window_size 1 \
		--cut_tail_mean_quality 3 \
		--cut_right --cut_right_window_size 4 \
		--cut_right_mean_quality 15 \
		--length_required 36 \
		--thread ${CPUs}
	fi
done; fi

# Trimming paired end reads  
if [ ${SEQMODE} == "paired" ]; then
	declare FILES1="${SOURCEDIR}*${PAIR1}${SUFFIX1}";
	for FILENAME in ${FILES1[@]}; do
	READ1=${FILENAME}	
	READ2=${FILENAME:0:-${#SUFFIX1}-${#PAIR1}}${PAIR2}${SUFFIX1}
	echo -e "[TRIMMING] fastp: [${READ1:${#SOURCEDIR}:-${#PAIR1}-${#SUFFIX1}}]"
	#Prevent overwrite:
	if [ -e ${TRIMM_DIR}${READ1:${#SOURCEDIR}:-${#SUFFIX1}}${SUFFIX_out}  ] && \
	   [ -e ${TRIMM_DIR}${READ2:${#SOURCEDIR}:-${#SUFFIX1}}${SUFFIX_out}  ]; then
		echo "Files exist, continuing"
	else
	fastp \
	--in1 ${READ1} \
	--in2 ${READ2} \
	--out1 ${TRIMM_DIR}${READ1:${#SOURCEDIR}:-${#SUFFIX1}}${SUFFIX_out} \
	--out2 ${TRIMM_DIR}${READ2:${#SOURCEDIR}:-${#SUFFIX1}}${SUFFIX_out} \
	--json "${QC_DIR_fastp}${READ1:${#SOURCEDIR}:-${#SUFFIX1}-${#PAIR1}}PE_fastp.json" \
	--html "${QC_DIR_fastp}${READ1:${#SOURCEDIR}:-${#SUFFIX1}-${#PAIR1}}PE_fastp.html" \
	--cut_front --cut_front_window_size 1 \
	--cut_front_mean_quality 3 \
	--cut_tail \
	--cut_tail_window_size 1 \
	--cut_tail_mean_quality 3 \
	--cut_right \
	--cut_right_window_size 4 \
	--cut_right_mean_quality 15 \
	--length_required 36 \
	--thread ${CPUs}
	fi
done; fi

################################
#INFORMATION ON TRIMMING PROCESS
#--cut_front --cut_front_window_size 1 --cut_front_mean_quality 3 == Trimmomatic  "LEADING:3"
#--cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 3 == Trimmomatic  "TRAILING:3"
#--cut_right --cut_right_window_size 4 --cut_right_mean_quality 15 == Trimmomatic  "SLIDINGWINDOW:4:15"

#### If additional trimming is needed (see multiQCreport): 
# add to R1: --trim_front1 {amount_bases} and/or --trim_tail1 {amount_bases}
# add to R2: --trim_front2 {amount_bases} and/or --trim_tail2 {amount_bases}

####################################################
### Alignment of reads (paired end & single end) ###
####################################################

if [ ${SEQTYPE} == "RNASeq" ]; then

#Indexing reference genome for STAR
cd ${GENOMEDIR}
if [ ${GenomeIndexDone} == "No" ]; then
	if [ ${LargeGenome} == "Yes" ]; then
		echo "Indexing Large(=Human) genome"
		# additional parameters to save RAM usage: --genomeSAsparseD 1 --genomeChrBinNbits 15
		STAR \
		--runMode genomeGenerate \
		--genomeDir ${GENOMEDIR} \
		--genomeFastaFiles ${GENOME} \
		--sjdbGTFfile ${GTF} \
		--sjdbOverhang 99 \
		--runThreadN ${CPUs} \
		--genomeSAsparseD 1 \
		--genomeChrBinNbits 15
	else
		echo "Indexing Small(=non-human) genome"
		STAR \
		--runMode genomeGenerate \
		--genomeDir ${GENOMEDIR} \
		--genomeFastaFiles ${GENOME} \
		--sjdbGTFfile ${GTF} \
		--sjdbOverhang 99 \
		--runThreadN ${CPUs}
	fi
fi
cd ${SOURCEDIR}

# Aligning reads single end
if [ ${SEQMODE} == "single" ]; then
	declare FILES1=${OUTPUTDIR}*${PAIR1}*${SUFFIX1}
  for FILENAME in ${FILES1[@]}; do
	READ1=${FILENAME:0:-${#SUFFIX1}}	
	#Prevent overwrite: 
	if [ -e ${align_DIR}${READ1:${#OUTPUTDIR}:-(${#PAIR1}+8)}Log.final.out  ]; then
		echo "File exists, continuing"
	else
		echo -e "[ALIGNING] STAR : [${READ1:${#OUTPUTDIR}:-${#PAIR1}}]"
		if [[ "${SUFFIX1}" == *".gz" ]]; then
			echo "gzipped file detected, using zcat to read"
			STAR \
			--runThreadN ${CPUs_align} \
			--genomeDir ${GENOMEDIR} \
			--readFilesIn "${READ1}${SUFFIX1}" \
			--quantMode TranscriptomeSAM \
			--readFilesCommand zcat \
			--outFileNamePrefix ${align_DIR}${READ1:${#OUTPUTDIR}:-(${#PAIR1}+8)}
		else
			echo "uncompressed FASTQ format detected"
			STAR \
			--runThreadN ${CPUs_align} \
			--genomeDir ${GENOMEDIR} \
			--readFilesIn "${READ1}${SUFFIX1}" \
			--quantMode TranscriptomeSAM \
			--outFileNamePrefix ${align_DIR}${READ1:${#OUTPUTDIR}:-(${#PAIR1}+8)}
		fi
	fi
done; fi

# Aligning reads paired end
if [ ${SEQMODE} == "paired" ]; then
	declare FILES1=${OUTPUTDIR}*${PAIR1}*${SUFFIX1}		
  for FILENAME in ${FILES1[@]}; do
	READ1=${FILENAME:0:-${#SUFFIX1}}	
	READ2=${FILENAME:0:-${#SUFFIX1}-(${#PAIR1}+8)}${PAIR2}"_trimmed"
	#Prevent overwrite: 
	if [ -e ${align_DIR}${READ1:${#OUTPUTDIR}:-(${#PAIR1}+8)}Log.final.out  ]; then
		echo "File exists, continuing"
	else
		echo -e "[ALIGNING] STAR : [${READ1:${#OUTPUTDIR}:-(${#PAIR1}+8)}]"
		if [[ "${SUFFIX1}" == *".gz" ]]; then
			echo "gzipped file detected, using zcat to read"
			STAR \
			--runThreadN ${CPUs_align} \
			--genomeDir ${GENOMEDIR} \
			--readFilesIn "${READ1}${SUFFIX1}" "${READ2}${SUFFIX1}" \
			--quantMode TranscriptomeSAM \
			--readFilesCommand zcat \
			--outFileNamePrefix ${align_DIR}${READ1:${#OUTPUTDIR}:-(${#PAIR1}+8)}
		else
			echo "uncompressed FASTQ format detected"
			STAR \
			--runThreadN ${CPUs_align} \
			--genomeDir ${GENOMEDIR} \
			--readFilesIn "${READ1}${SUFFIX1}" "${READ2}${SUFFIX1}" \
			--quantMode TranscriptomeSAM \
			--outFileNamePrefix ${align_DIR}${READ1:${#OUTPUTDIR}:-(${#PAIR1}+8)}
		fi
	fi
done; fi

#######################
# QUANTIFICATION RSEM #
#######################

#Indexing reference genome for RSEM
echo "Indexing genome"
cd ${GENOMEDIR}
if [ ${GenomeIndexRSEMDone} == "No" ]; then 
	rsem-prepare-reference --gtf ${GTF} ${GENOME} ${RSEM_GENOMEDIR}${GenomeID}
	GenomeIndexRSEMDone="Yes"
fi

echo "Quantifying single end"
cd ${SOURCEDIR}
# Quantify reads single end
if [ ${SEQMODE} == "single" ]; then
	declare FILES1=${OUTPUTDIR}*${PAIR1}*${SUFFIX1}		
for FILENAME in ${FILES1[@]}; do
	READ1=${FILENAME:0:-${#SUFFIX1}}	
	echo -e "[QUANTIFYING] RSEM : [${READ1:${#OUTPUTDIR}:-(${#PAIR1}+8)}]" 
	rsem-calculate-expression \
	-p ${CPUs} \
	--bam "${align_DIR}${READ1:${#OUTPUTDIR}:-(${#PAIR1}+8)}Aligned.toTranscriptome.out.bam" \
	--no-bam-output ${RSEM_GENOMEDIR}${GenomeID} ${Quant_DIR}${READ1:${#OUTPUTDIR}:-(${#PAIR1}+8)}
done;fi

echo "Quantifying paired end"
# Quantify reads paired end
if [ ${SEQMODE} == "paired" ]; then
	declare FILES1=${OUTPUTDIR}*${PAIR1}*${SUFFIX1}		
for FILENAME in ${FILES1[@]}; do
	READ1=${FILENAME:0:-${#SUFFIX1}}	
	READ2=${FILENAME:0:-${#SUFFIX1}-(${#PAIR1}+8)}${PAIR2}"_trimmed"
	echo -e "[QUANTIFYING] RSEM : [${READ1:${#OUTPUTDIR}:-(${#PAIR1}+8)}]" 
	rsem-calculate-expression \
	-p ${CPUs} \
	--paired-end \
	--bam "${align_DIR}${READ1:${#OUTPUTDIR}:-(${#PAIR1}+8)}Aligned.toTranscriptome.out.bam" \
	--no-bam-output ${RSEM_GENOMEDIR}${GenomeID} ${Quant_DIR}${READ1:${#OUTPUTDIR}:-(${#PAIR1}+8)}
done;fi

# Output summary file from RSEM results
declare FILELIST=$(find ${Quant_DIR} -name "*genes.results"  -printf "%f\t")
rsem-generate-data-matrix $FILELIST > ${Quant_DIR}/genes.data.tsv
sed -i 's/\.genes.results//g' ${Quant_DIR}/genes.data.tsv

declare FILELIST=$(find ${Quant_DIR} -name "*isoforms.results"  -printf "%f\t")
rsem-generate-data-matrix $FILELIST > ${Quant_DIR}/isoforms.data.tsv
sed -i 's/\.genes.results//g' ${Quant_DIR}/isoforms.data.tsv

else
if [ ${SEQTYPE} == "TempOSeq" ]; then
# Command Line Arguments:
# 1: FASTA Reference File
# 2: Directory of FASTQ files to align
# 3: Number of CPUs to use
cd ${BASEDIR}
echo "Parameters passed to Rscript for TempO-Seq Alignments..."
echo ${TEMPOSEQR}
echo ${GENOME}
echo ${TRIMM_DIR}
echo ${CPU_FOR_ALIGNMENT}
Rscript ${TEMPOSEQR} ${GENOME} ${TRIMM_DIR} ${CPU_FOR_ALIGNMENT}
declare READCOUNTS="${BASEDIR}/count_table.csv"
declare MAPPEDUNMAPPED="${BASEDIR}/mapped_unmapped.csv"
echo ${READCOUNTS}
echo ${MAPPEDUNMAPPED}
Rscript ${TEMPOSEQFIGS} ${READCOUNTS} ${MAPPEDUNMAPPED}
fi
fi

###################################################################################################
### Quality control raw reads: Fastp + RSEM + STAR MultiQC report ###
#####################################################################

# Run MultiQC on outputs
multiqc \
--cl_config "extra_fn_clean_exts: { '_fastp.json' }" \
--cl_config "sample_names_replace_exact: True" \
${BASEDIR} \
--filename MultiQC_Report.html \
--interactive \
--replace-names "${RAW_SAMPLE_DIR}/../../metadata/metadata.multiqc.txt" \
--outdir ${QC_DIR_multiQC}

###################################################################################################

conda env export > conda_environment.${mydate}.yml

conda deactivate
echo "Pre-processing of data complete."

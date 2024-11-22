#!/bin/bash -l

#SBATCH -A mandrup    # Names the account for tracking usage, will stay as your lab name
#SBATCH --output slurm-%x.%A.%a.log   # Names an output file for what would go to STDOUT, %x %A %a represent jobname jobid jobarrayindex
#SBATCH --mail-user victorg@bmb.sdu.dk   # Names an address to be emailed when the job finishes
#SBATCH --mail-type END,FAIL,ARRAY_TASKS  # Specifies when the job finishes (either correctly or failed)
#SBATCH --job-name this_job   # Gives the job a name, so you can find its log file & see it in the queue status, etc
#SBATCH --nodes 1         # How many nodes to be used, ie one compute nodes
#SBATCH --mem 90G        # The job can use up to __GB of ram. It is mutually exclusive with --mem-per-cpu.
#SBATCH -p CLOUD       # Names to use the serial partition
#SBATCH --cpus-per-task 16    # How many cores on that node
##SBATCH --mem-per-cpu 2500M   # he job can use up to 2.5 GB (non-integers are not allowed) of ram per cpu, i.e. 80 GB ram. NOT IN USE
#SBATCH -t 20:30:00       # Means to use up to hours:minutes:seconds of run time before it will get killed

# run in u1-standard-16 node

# Start runtime
START=$(date +%s)
echo -e "\nStarting processing"

#we set OMP_NUM_THREADS to the number of available cores
echo "Running on $SLURM_CPUS_ON_NODE CPU cores"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# If using #SBATCH --cpus-per-task then
#we set MEM_PER_THREADS to the max memory per CPU
#echo "Running on $SLURM_MEM_PER_CPU M max RAM per CPU"
#export MEM=$((${SLURM_CPUS_PER_TASK}*${SLURM_MEM_PER_CPU}))
#echo "Running on ${MEM} M max RAM"


# else using #SBATCH --mem then
echo "Running on $SLURM_MEM_PER_NODE M max RAM"

# One should only specify 80-90% of available memory
MEM=$((${SLURM_MEM_PER_NODE}/1024)) # Memory in GB
MEM=$((9*${MEM}/10)) # Leave 19 GB ram aside of the available RAM
MEM_PER_THREAD=$((${MEM}/${SLURM_CPUS_PER_TASK}))

echo "${MEM} G max total mem"
echo "${MEM_PER_THREAD} G max mem per thread"

## Read optional arguments
function usage {
  echo -e "\n Usage:$(basename $0) -g <genome>  <input_files>"
  echo "Options:"
  echo " -g <genome>   - Specify the genome (mm10, mm39, hg38)"
  echo " -h            - Display this help message"
  exit 1
}

while getopts g:h opt; do
    case "${opt}" in
      g) GENOME="${OPTARG}"
      ;;
      \?) ## Invalid option
      echo "Invalid option: -${OPTARG}" >&2
      usage
      exit 1
      ;;
      :) ## Argument required
      echo "Option -${OPTARG} requires an argument." >&2
      usage
      exit 1
      ;;
    esac
done
shift "$((OPTIND-1))"  #This tells getopts to move on to the next argument.

if [ -z "${GENOME}" ]; then
  echo "Genome option (-g) is required."
  usage
  exit 1
fi

echo "Using genome setup: ${GENOME}"

# PATH to the reference genomes
if [ "${GENOME}" == "mm10" ]; then
	REFERENCE=/work/References/Mouse/mm10/bowtie1.3.1/mm10
	MIRNA_PRECURSOR=/work/References/mirbase/RELEASE_22.1/premature_mmu.fa
	MIRNA_MATURE=/work/References/mirbase/RELEASE_22.1/mature_mmu.fa

elif [ "${GENOME}" == "mm39" ]; then
  REFERENCE=/work/References/Mouse/mm39/bowtie1.3.1/mm39
  MIRNA_PRECURSOR=/work/References/mirbase/RELEASE_22.1/premature_mmu39.fa
  MIRNA_MATURE=/work/References/mirbase/RELEASE_22.1/mature_mmu.fa

elif [ "${GENOME}" == "hg38" ]; then
  REFERENCE=/work/References/Human/hg38_analysisSet/bowtie1.3.1/hg38.analysisSet
  MIRNA_PRECURSOR=/work/References/mirbase/RELEASE_22.1/premature_hsa_hg38.fa
  MIRNA_MATURE=/work/References/mirbase/RELEASE_22.1/mature_hsa.fa

else
  echo "Invalid genome option: ${GENOME}"
  usage
  exit 1
fi

echo "Using reference: ${REFERENCE}"

# Input Files
if [ $# -eq 0 ]; then
  echo "No input files provided."
  usage
fi

INPUT1=${1?Missing input R1_001.fastq.gz file}
INPUT1=$(readlink -f "${INPUT1}")

echo "Input file: ${INPUT1}"

# ===============================================
# Create Main Output Directory
# ================================================

PREFIX=$(basename "${INPUT1}" _R1_001.fastq.gz)

if [ -f "${PREFIX}" ]; then
  >&2 echo "Error: Output location (${PREFIX}) already exists as a file"
  exit 1
fi

if [ -d "${PREFIX}" ]; then
  echo "Warning: Output location (${PREFIX}) is already a directory, reusing, could overwrite"
  # If you don't want to reuse, you could make it exit 1 here to kill the script if
  # the folder already exists
else
  mkdir "${PREFIX}"
fi

cd "${PREFIX}"

# Output Files
TRIMMED1=$(basename "${INPUT1}" _R1_001.fastq.gz)_trimmed_R1_001.fq # Output AdapterRemoval
FINAL_COUNTS="${PREFIX}_miRNA_RawCounts.txt"

# Commands
echo "AdapterRemoval..."
START_SUBPROCESS=$(date +%s)

module load AdapterRemoval

ILLUMINA=$(echo 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA')  #Default, no need to specify
SMALL_RNA=$(echo 'TGGAATTCTCGGGTGCCAAGG')
Nextera=$(echo 'CTGTCTCTTATACACATCT')

# Command
# For Nextflex V3
# --trim3p 4
# --minlength default 15
# AdapterRemoval --threads "${OMP_NUM_THREADS}" --trim3p 4 --trim5p 4 --minlength 10 --file1 "${INPUT1}" --basename "${PREFIX}" --adapter1 "${SMALL_RNA}" --output1 "${TRIMMED1}"

# For Nextflex V4
AdapterRemoval --threads "${OMP_NUM_THREADS}" --file1 "${INPUT1}" --basename "${PREFIX}" --minlength 16 --adapter1 "${SMALL_RNA}" --output1 "${TRIMMED1}"

END_SUBPROCESS=$(date +%s)
RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

# Alignment
echo "miRDeep2...mapper.pl"
START_SUBPROCESS=$(date +%s)

module load miRDeep2

# -o 			# number of threads to use for bowtie
# -e            # input file is fastq format
# -l 18			# discard reads shorter than int nts, default = 18
# -h            # parse to fasta format
# -m            # collapse reads
# -s file       # print processed reads to this file [output file]
# -t file       # print read mappings to this file [output file]
mapper.pl "${PREFIX}"_trimmed_R1_001.fq -o "${OMP_NUM_THREADS}" -e -j -h -m -v -p "${REFERENCE}" -s ${PREFIX}.processed.fa -t ${PREFIX}.mapped.arf

END_SUBPROCESS=$(date +%s)
RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

# Quantifying
echo "miRDeep2...quantifier.pl"
START_SUBPROCESS=$(date +%s)
# -p precursor.fa  # miRNA precursor sequences from miRBase
# -m mature.fa     # miRNA sequences from miRBase
# -c c [file]    config.txt file with different sample ids...or just the one sample id
# -d to leave out pdfs
quantifier.pl -d -p "${MIRNA_PRECURSOR}" -m "${MIRNA_MATURE}" -r ${PREFIX}.processed.fa

mv ./miRNAs_expressed_all_samples_*.csv ./"${PREFIX}_miRNAs_expressed.txt"
mv ./expression_*.html ./"${PREFIX}_expression.html"
mv ./expression_analyses/expression_analyses_*/* ./expression_analyses/
mv ./expression_analyses ./"${PREFIX}_expression_analyses"

cat "${PREFIX}"_miRNAs_expressed.txt | cut -f1-3 | awk '{ print $1 "\t" $3 "\t" $2}' > "${FINAL_COUNTS}" # Three tab separated Columns miRNA \t precursor \t read_count

END_SUBPROCESS=$(date +%s)
RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

chmod 770 ./*
chmod 770 ./*/*
chmod 770 ./*/*/*

# Finalize
END=$(date +%s)
RUNTIME=$((END-START))
H=$((RUNTIME / 3600 ))  # Calculate hours
M=$(( (RUNTIME / 60 ) % 60 ))  # Calculate minutes
S=$(( RUNTIME % 60 ))  # Calculate seconds
echo -e "\tProcessing completed. Total run time: ${H} hours, ${M} minutes, and ${S} seconds."

echo "Run multiqc with option -e bowtie1"
echo "...Done!"

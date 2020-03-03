#!/bin/bash

# ------------------------------------------------------------------------------------------------ #
# depletePhiX.sh
#     deplete PhiX spike-in reads from FASTQ files that were added during library prep
# ------------------------------------------------------------------------------------------------ #
echo -e "Welcome to 'depletePhiX.sh' !"

# If one step fails, stop the script and exit
set -eo pipefail

# ------------------------------------------------------------------------------------------------ #
# Ensure the script is called correctly
# ------------------------------------------------------------------------------------------------ #
# Set up a usage statement in case this program is called incorrectly
usage() { echo -e "\nERROR: Missing sample name (needed for naming files) or input files. \n" \
                  "Proper usage for trimming adapters from paired-end reads: \n\n" \
                  "$0 -s mysample -1 sample-reads_R1.fq -2 sample-reads_R2.fq \n\n" \
                  "Optional parameter: \n" \
                  "-o (output directory for saving trimmed files; [default = './' (current directory)]) \n\n" \
                  "Example of a complex run: \n" \
                  "$0 -s my_sample -1 sample-reads_R1.fq -2 sample-reads_R2.fq -o trimmed_reads/ \n\n" \
                  "Exiting program. Please retry with corrected parameters..." >&2; exit 1;
        }

# Make sure the pipeline is invoked correctly, with project and sample names
while getopts "s:o:p:1:2:" arg; do
        case ${arg} in
                s ) # Take in the sample name for naming
                  sampleID=${OPTARG}
                        ;;
                p) # path to phiX genome
                  phiX_genome=${OPTARG}
                        ;;
                1 ) # path to forward reads fastq
                  trimmed_fwd_reads=${OPTARG}
                        ;;
                2 ) # path to reverse reads fastq
                  trimmed_rev_reads=${OPTARG}
                        ;;
            o ) # set the output directory
                  output_directory=${OPTARG}
                        ;;
                * ) # Display help
                  usage
                        ;;
        esac
done
shift $(( OPTIND-1 ))

# Check that required parameters are provided
if [[ -z "${sampleID}" ]] || [[ -z "${phiX_genome}" ]] || \
   [[ -z "${trimmed_fwd_reads}" ]] || [[ -z "${trimmed_rev_reads}" ]]; then
    usage
fi

# Set up an empty log file
cat /dev/null > ${sampleID}.phiX_depletion.log

# Create an output directory to store the output files
## If user didn't provide one, just use a subdirectory in working directory
if [[ -z ${output_directory} ]]; then
    output_directory="./"
    mkdir -p ${output_directory}

else
    # If user provided a desired output directory: check to make sure output directory doesn't exist;
    # then create output directory; if error, just default to a results subdirectory within current dir
    if [[ ! -d ${output_directory} ]]; then
        mkdir -p ${output_directory} || \
        {
            echo "Cannot create user-provided output directory. Defaulting to current working directory './'" | \
            tee ${sampleID}.phiX_depletion.log

            output_directory="./"
            mkdir ${output_directory}
        }
    fi
fi
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #

function deplete_phiX() {

    # --FUNCTION---------------------------------------------------------------------------------- #
    # Name:        deplete_phiX
    # Description: Remove PhiX spike-in reads from sample
    # -------------------------------------------------------------------------------------------- #

    # -------------------------------------------------------------------------------------------- #
    # Ensure that the necessary software is installed
    command -v bwa > /dev/null || {
        echo -e "ERROR: This script requires 'bwa' but it could not found. \n" \
                "Please install this application. \n" \
                "Exiting with error code 6..." >&2; exit 6
        }

    command -v samtools > /dev/null || {
        echo -e "ERROR: This script requires 'samtools' but it could not found. \n" \
                "Please install this application. \n" \
                "Exiting with error code 6..." >&2; exit 6
        }
    # -------------------------------------------------------------------------------------------- #

    # -------------------------------------------------------------------------------------------- #
    # Adapter trimming log info
    echo "Began adapter trimming at:    $(date)" | \
    tee -a ${SAMPLE}.phiX_de.log
    # -------------------------------------------------------------------------------------------- #

    # -------------------------------------------------------------------------------------------- #
    # Deplete phiX
    # -------------------------------------------------------------------------------------------- #
    # Build BWA index out of the reference
    bwa index \
    -p phiX_index \
    $phiX_genome

    # Perform the mapping
    bwa mem \
    phiX_index \
    $trimmed_fwd_reads $trimmed_rev_reads > ${sampleID}.phiX_mapped.sam

    # Get summary stats of the mapping
    samtools flagstat ${sampleID}.phiX_mapped.sam > \
    ${sampleID}.phiX.stats.txt

    # Retrieve phiX-depleted reads & sort output
    samtools view -f 4 -bh ${sampleID}.phiX_mapped.sam | \
    samtools sort > \
    ${sampleID}.phiX_depleted.sorted.bam

    # Save how many phiX reads were in the sequencing data
    samtools idxstats ${sampleID}.phiX_depleted.sorted.bam > \
        ${sampleID}.phiX.counts.txt

    # Print that phiX count to the screen
    echo "name   length   mapped_reads   unmapped_reads"
    cat ${sampleID}.phiX.counts.txt

    # Convert the phiX-depleted BAM back to individual FASTQ files
    samtools fastq \
    -1 ${sampleID}.phiX_depleted_R1.fq \
    -2 ${sampleID}.phiX_depleted_R2.fq \
    -s ${sampleID}.phiX_depleted.singletons.fq \
    ${sampleID}.phiX_depleted.sorted.bam

    # Compress the phiX-depleted fastqs
    gzip ${sampleID}.phiX_depleted*fq
    # -------------------------------------------------------------------------------------------- #

    # -------------------------------------------------------------------------------------------- #
    # Adapter trimming log info
    echo "Finished adapter trimming at:    $(date)" | \
    tee -a ${SAMPLE}.phiX_depletion.log
   # -------------------------------------------------------------------------------------------- #

}

deplete_phiX

#!/bin/bash

# ------------------------------------------------------------------------------------------------ #
# mapToTVV.sh
#     Take adapter-trimmed, phiX-depleted dsRNA-seq reads and map to to the five TVV species + satellites
# ------------------------------------------------------------------------------------------------ #
echo -e "Welcome to 'mapToTVV.sh' !"

# If one step fails, stop the script and exit
set -eo pipefail

# ------------------------------------------------------------------------------------------------ #
# Ensure the script is called correctly
# ------------------------------------------------------------------------------------------------ #
# Set up a usage statement in case this program is called incorrectly
usage() { echo -e "\nERROR: Missing sample name (needed for naming files) or input files. \n" \
                  "Proper usage for trimming adapters from paired-end reads: \n\n" \
                  "$0 -s mysample -1 sample-reads_R1.fq -2 sample-reads_R2.fq -t tvv_directory/ \n\n" \
                  "Optional parameter: \n" \
                  "-o (output directory for saving trimmed files; [default = './' (current directory)]) \n\n" \
                  "Example of a complex run: \n" \
                  "$0 -s my_sample -1 sample-reads_R1.fq -2 sample-reads_R2.fq -t tvv_genomes/ -o './' \n\n" \
                  "Exiting program. Please retry with corrected parameters..." >&2; exit 1;
        }

# Make sure the pipeline is invoked correctly, with project and sample names
while getopts "s:o:t:1:2:" arg; do
        case ${arg} in
                s ) # Take in the sample name for naming
                  sample=${OPTARG}
                        ;;
                1 ) # path to forward reads fastq
                  forward_reads=${OPTARG}
                        ;;
                2 ) # path to reverse reads fastq
                  reverse_reads=${OPTARG}
                        ;;
                t ) # path to directories containing the TVV reference genomes
                  tvv_directory=${OPTARG}
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
if [[ -z "${sample}" ]] || [[ -z "${forward_reads}" ]] || [[ -z "${reverse_reads}" ]]; then
    usage
fi

# Set up an empty log file
cat /dev/null > ${sample}.mapToTVV.log

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
            tee ${sample}.mapToTVV.log

            output_directory="./"
            mkdir ${output_directory}
        }
    fi
fi
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #

function mapping () {

    # --FUNCTION---------------------------------------------------------------------------------- #
    # Name:        mapping
    # Description: Map to dsRNA-seq reads to all TVVs and satellites
    # -------------------------------------------------------------------------------------------- #

    # -------------------------------------------------------------------------------------------- #
    # Ensure that the necessary software is installed
    command -v bbsplit.sh > /dev/null || {
        echo -e "ERROR: This script requires 'bbsplit.sh' but it could not found. \n" \
                "Please install this application. \n" \
                "Exiting with error code 6..." >&2; exit 2
        }
    # -------------------------------------------------------------------------------------------- #

    # -------------------------------------------------------------------------------------------- #
    # Adapter trimming log info
    echo "Began mapping to TVVs at:    $(date)" | \
    tee -a ${sample}.mapToTVV.log
    # -------------------------------------------------------------------------------------------- #

    # -------------------------------------------------------------------------------------------- #
    # Set up TVV genomes for mapping
    tvv1="${tvv_directory}/tvv1.fasta"
    tvv2="${tvv_directory}/tvv2.fasta"
    tvv3="${tvv_directory}/tvv3.fasta"
    tvv4="${tvv_directory}/tvv4.fasta"
    tvv5="${tvv_directory}/tvv5.fasta"
    tvv_satellites="${tvv_directory}/tvv-dsRNA-satellites.fasta"

    # Make sure the provided TVV genomes are present
    if [[ ! -f $tvv1 || ! -f $tvv2 || ! -f  $tvv3 || ! -f  $tvv4 || ! -f $tvv5 || ! -f  $tvv_satellites ]];
        then echo "Missing one or more TVV genomes files. " \
                  "Ensure all the files are present for TVV1-TVV5 + TVV-dsRNA-satellites and " \
                  "in the indiciated directory: ${tvv_directory}" | tee -a ${sample}.mapToTVV.log
             exit 3
     fi
    # -------------------------------------------------------------------------------------------- #

    # -------------------------------------------------------------------------------------------- #
    # Perform the mapping
    # -------------------------------------------------------------------------------------------- #
    bbsplit.sh \
        in1=$forward_reads \
        in2=$reverse_reads \
        ref="${tvv1},${tvv2},${tvv3},${tvv4},${tvv5},${tvv_satellites}" \
        basename="${sample}_%.R#.fq" \
        outu="${sample}_unmapped-reads.R#.fq" \
        ambig2=best \
        1>> "${sample}.mapToTVV.log" \
        2>> "${sample}.mapToTVV.log"
    # -------------------------------------------------------------------------------------------- #

    # -------------------------------------------------------------------------------------------- #
    # Adapter trimming log info
    echo "Finished adapter trimming at:    $(date)" | \
    tee -a ${sample}.mapToTVV.log
   # -------------------------------------------------------------------------------------------- #

}
# ------------------------------------------------------------------------------------------------ #

mapping

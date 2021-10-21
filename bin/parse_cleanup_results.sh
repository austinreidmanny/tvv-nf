#!/bin/bash

# ----------------------------------------------------
# parse_cleanup_results.sh
# ----------------------------------------------------
# This script reads the results of cleanup_blast.sh
# and formats the results into a nice summary table
# ----------------------------------------------------


# ----------------------------------------------------
# Ensure this script is called correctly
# ----------------------------------------------------
set -eo pipefail


# ----------------------------------------------------
# Create a usage statement
# ----------------------------------------------------
usage() { echo -e "# ---------------------------------------------------- # \n" \
				  "This script will read the results of cleanup_blast.sh and format the results. \n\n" \
	              "Proper usage: $0 -f cleanup_results.txt \n\n" \
                  "Please retry with correct parameters. Exiting with error code 1... \n" \
				  "# ---------------------------------------------------- # \n"
		  exit 1 
	  }


# ----------------------------------------------------
# Take in cleanup_blast results.txt file
# ----------------------------------------------------
# Read user-provided parameters
while getopts "f:*:" arg; do
	case ${arg} in
		f ) # Take in the file
			input=${OPTARG} ;;

		* ) # if anything else is given, just print usage statement
			usage ;;

	esac
done

shift $(( OPTIND - 1 ))

# If no input is provided, tell that to the user and exit
if [[ -z $input ]] ; then
   	echo -e "\nERROR: No input detected. \n"
	usage
fi



# ----------------------------------------------------
# Define the main code
# ----------------------------------------------------
main() {

	# Name an output file, based on the input file, saved in the current directory
	output="$(basename $input '.txt').parsed.txt"

	# Create a sample name, so we know which input file was used
	sample="$(basename $input | cut -d "." -f 1)"

	# Count the number of each TVV sequence in the file and store those values

		# (note: if grep does not find the pattern (e.g., any TVV1 reads), it will return an error 1 code
		# which kills the program; the last command checks if there is an error code of 1, an ignores it
		# this will still allow grep to give a true error code, which will be >1 (according to documentation)

 	tvv1_count=$(cat $input | grep -c "TVV1\|Trichomonas_vaginalis_virus_1" || [[ $? == 1 ]] )
    tvv2_count=$(cat $input | grep -c "TVV2\|Trichomonas_vaginalis_virus_2\|Trichomonas_vaginalis_virus_II" || [[ $? == 1 ]] )
	tvv3_count=$(cat $input | grep -c "TVV3\|Trichomonas_vaginalis_virus_3" || [[ $? == 1 ]] )
	tvv4_count=$(cat $input | grep -c "TVV4\|Trichomonas_vaginalis_virus_4" || [[ $? == 1 ]] )
	tvv5_count=$(cat $input | grep -c "TVV5\|Trichomonas_vaginalis_virus_5" || [[ $? == 1 ]] )

	# Print each species with each count
	echo -e "${sample}\tTVV1\t${tvv1_count}" >  $output
	echo -e "${sample}\tTVV2\t${tvv2_count}" >> $output
	echo -e "${sample}\tTVV3\t${tvv3_count}" >> $output
	echo -e "${sample}\tTVV4\t${tvv4_count}" >> $output
	echo -e "${sample}\tTVV5\t${tvv5_count}" >> $output

	}	


# ----------------------------------------------------
# Run the main code
# ----------------------------------------------------
main


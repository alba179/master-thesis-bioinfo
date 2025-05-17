#!/bin/bash
# 09/10/2024
# Author: Alba Méndez Alejandre
# Description: This script is going to perform vep annotation for a input release (as 102, for example) a input species, and a defined file paths using .fofn
# We are going to use the offline cache, and output a vcf annotated with "_vep" suffix
# Usage: 
# ./annotate_vep.sh 102 mus_musculus /dir/to/.fofn

release=$1 # input only the release number, as 102
species=$2 # species name, such as mus_musculus or homo_sapiens
file_of_filenames=$3 # the third input is the .fofn file, that are in our workstation

# check to ensure correct input arguments
if [ $# -ne 3 ]; then
    echo "Usage: $0 <release> <species> <file_of_filenames>"
    exit 1
fi

# check to ensure correct full file paths)
if [ ! -f "$file_of_filenames" ]; then
    echo "The file $file_of_filenames does not exist. Exiting."
    exit 1
fi

# check to ensure conda env exists
if ! conda info --envs | grep "vep_$release"; then
    echo "Conda environment vep_$release does not exist. Exiting."
    exit 1
fi

source ~/miniforge3/etc/profile.d/conda.sh # we always have to add this in bash scripts for conda activate
conda activate vep_$release # activate conda env with that specific release


while read -r full_path; do # this reads the input list
	
	vcf_file=$(basename $full_path) # we extract the basename of the fullpath for our file we want to annotate
	dir_workst=/media/storage/mcGinn_2021/scomatic/vep_annotation_final

	printf "\n[ $(date +'%Y/%m/%d %T.%3N') ] Processing files in $vcf_file_workst \n\n"
        
	dir_folia=/home/albax/mcGinn_2021/vep_annotation/
	mkdir -p dir_folia

	input_file_folia=$dir_folia/$vcf_file

	rsync --progress -hh -av matterhorn:${full_path} $dir_folia

	# we want the annotated vcf file to be in the same directory as the original vcf file
	output_file_folia=${vcf_file%.vcf}.vep.vcf # remove the vcf and append .vep to the output file
	
	if 
		vep -i $input_file_folia -o ${dir_folia}/${output_file_folia} --fork 14 --cache --dir_cache ~/.vep --species $species --cache_version $release --offline --everything --vcf; then
		printf "[ $(date +'%Y/%m/%d %T.%3N') ] Successfully annotated $vcf_file_workst"
	else
		printf "[ $(date +'%Y/%m/%d %T.%3N') ] Failed to annotate $vcf_file_workst"
	fi
	
	rm -f ${dir_folia}/${vcf_file}

	rsync --progress -hh -av ${dir_folia}/*.* matterhorn:$dir_workst

	rm -rf $dir_folia

done < $file_of_filenames # pass the .fofn file to read, with the list of full paths to process and run vep into


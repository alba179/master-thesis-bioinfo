#!/bin/bash
# [folia]
# This script must be executed in the server where we want to run STAR

# activate env where we have installed STAR
source /home/albax/miniforge3/bin/activate STAR

FOLIA_BASE="/home/albax/mcGinn_2021"
MATTERHORN_BASE="/media/storage/mcGinn_2021"
CONCATENATED_RUNS_DIR="${MATTERHORN_BASE}/concatenated_runs/sample_concatenated"


SAMPLE_NAMES=("SIGAD8" "SIGAE8" "SIGAF6" "SIGAF8" "SIGAG6" "SIGAG8" "SIGAH6" "SIGAH8") # list of samples to proccess

i=0
for SAMPLE in "${SAMPLE_NAMES[@]}"; 
do
	i=$((i + 1))

	SAMPLE_DIR="${CONCATENATED_RUNS_DIR}/${SAMPLE}"
	LOCAL_SAMPLE_DIR="${FOLIA_BASE}${SAMPLE_DIR}"
	echo "Processing sample: ${SAMPLE}"

	R2_FILE="${SAMPLE_DIR}/${SAMPLE}*R2*run1.fastq.gz"
	R1_FILE="${SAMPLE_DIR}/${SAMPLE}*R1*run1.fastq.gz"

	# Files' pull from matterhorn to folia
	rsync -avR --progress -hh "alba@matterhorn:${R2_FILE}" "${FOLIA_BASE}"
	rsync -avR --progress -hh "alba@matterhorn:${R1_FILE}" "${FOLIA_BASE}"

	LOCAL_R2_FILE="${FOLIA_BASE}${R2_FILE}"
	LOCAL_R1_FILE="${FOLIA_BASE}${R1_FILE}"

	echo "R2 file (folia): ${LOCAL_R2_FILE}"
    echo "R1 file (folia): ${LOCAL_R1_FILE}"
	
	cd ${FOLIA_BASE}
	mkdir -p "${FOLIA_BASE}/STAR_out/${SAMPLE}"

  # Process files
#	~/miniforge3/envs/STAR/bin/STAR \
	STAR --runMode alignReads --genomeDir /home/albax/reference_genomes/STAR_indexes/Mus_musculus/mm10 --outFileNamePrefix STAR_out/SLX-17937_${SAMPLE} --readFilesCommand zcat --readFilesIn ${LOCAL_R2_FILE} ${LOCAL_R1_FILE} --soloType CB_UMI_Simple --soloCBwhitelist /home/albax/mcGinn_2021/STAR_whitelist/3M-february-2018.txt --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 --soloStrand Forward --soloFeatures Gene GeneFull Velocyto --outSAMattributes NH HI AS nM CR UR CB UB GX GN sM sS sQ --outSAMtype BAM SortedByCoordinate --runThreadN 4

	if [[ $? -ne 0 ]]; then
		echo "STAR alignment failed for sample ${SAMPLE}."
		continue
	fi

        # Transfer back to matterhorn
	ssh alba@matterhorn "mkdir -p /media/storage/mcGinn_2021/STARalignment/${SAMPLE}/"

	rsync -av --progress -hh ${FOLIA_BASE}/STAR_out/*${SAMPLE}* alba@matterhorn:${MATTERHORN_BASE}/STARalignment/${SAMPLE}/

	# Remove to save space
	echo -e "\n[ $(date +'%Y/%m/%d %T.%3N') ] Finished processing ${SAMPLE}"
	rm ${LOCAL_R1_FILE} ${LOCAL_R2_FILE}
	rm -r ${FOLIA_BASE}/STAR_out/*${SAMPLE}*/

done

# clean up base directory
rm -r ${FOLIA_BASE}/STAR_out

conda deactivate
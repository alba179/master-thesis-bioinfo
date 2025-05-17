#!/bin/bash
# [matterhorn]
# AUTHOR: Alba Mendez Alejandre
# DESCRIPTION: Download the raw fastQ via ftp from ENA. The metadata is stored in the downloaded csv. 
# DATE: 10/06/2024

seq_ids=("run1" "run2")

for seq_number in "${seq_ids[@]}"; do

    # Path to our CSV file
    CSV_FILE="/mnt/D/mcGinn_2021/data_info/${seq_number}_download.csv"

    # Directory where we are going to download the data (bam and fastq)
    DOWNLOAD_DIR="./${seq_number}_run1/${seq_number}_adult_P70"

    # Create the directory if it doesn't exist
    mkdir -p "$DOWNLOAD_DIR"

    # Iterate through each line in the CSV file
    { 
        read # skip header
        while IFS=, read -r name source_name bam_uri R1_fastq_uri R2_fastq_uri I1_fastq_uri # these are the columns of the file
        do
            # Skip header
            if [[ $name == "name" ]]; then
            continue
            fi

            # Ensure the URI is not empty
            if [[ -z $R1_fastq_uri ]]; then
            echo "Empty URI for name: $name"
            continue
            fi

            # Download R1_fastq_uri
            curl -L -o "${DOWNLOAD_DIR}/${name}_$(basename "$R1_fastq_uri")" "$R1_fastq_uri"
            
            if [[ -z $R2_fastq_uri ]]; then
            echo "Empty URI for name: $name"
            continue
            fi

            # Download R2_fastq_uri
            curl -L -o "${DOWNLOAD_DIR}/${name}_$(basename "$R2_fastq_uri")" "$R2_fastq_uri"
            
            if [[ -z $I1_fastq_uri ]]; then
            echo "Empty URI for name: $name"
            continue
            fi

            # Download I1_fastq_uri
            curl -L -o "${DOWNLOAD_DIR}/${name}_$(basename "$I1_fastq_uri")" "$I1_fastq_uri"
        
            if [[ -z $bam_uri ]]; then
            echo "Empty URI for name: $name"
            continue
            fi

        done < "$CSV_FILE"
    }
done

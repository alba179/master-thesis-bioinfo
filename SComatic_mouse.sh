#!/bin/bash
# Script with bams separated by run (SLX-17923 vs SLX-18123). Each run bams are analyzed in a scomatic pipeline. This means, the output of this script will give two final .tsv: one for the run of SLX-17923, and another one for the run of SLX-18123. Each of these runs/.tsv, will have 8 samples: sample1_Epi, samp2_Epi, samp3_Epi, samp4_Epi, samp5_Epi, samp6_Epi, samp7_Epi, samp8_Epi. 
# There are 6 sample which are control adult, + one sample is old CTL, and + one sample is old DEN. 
# Each of the sample has the key like "SIGAD8", "SIGAG6", etc. 
# DISCLAIMER: ALWAYS USE THE SAME REFERENCE GENOME!!! If you have aligned you rreads with mm10, scomatic will have to be ran with mm10. No coordinates' change with bedtools nor similar strange things. 
# SComatic nowadays has only panel of normales (PON) and editing sites (A-to-I events) for: GRCh38 and mm10. There aren't PON nor editing sites for mouse assembly mm39. 
# The PON for human has been made with GATK data 1000 genomes project. No info of how they did it.
# The PON for mouse mm10 was done in July 20, using mm10 Tabula muris data. No info of how exactly they did it. 

# Run this script in a environment with scomatic installed. 

# conda activate SComatic

SCOMATIC=~/bin/SComatic-main

run_names=("SLX-17937" "SLX-18123") # list of runs

echo "Starting Scomatic analysis."


i=0
for run in "${run_names[@]}";
do
	i=$((i + 1)) 

  # iterate though each bam in the directory where they are located
  bam_dir="/media/storage/mcGinn_2021/STARalignment/bam/${run}"
  output_dir="/media/storage/mcGinn_2021/scomatic/output/${run}"

  mkdir -p $output_dir

  logfile="$output_dir/scomatic_log.txt"
  mkdir -p "$(dirname "$logfile")"

  log_message() {
      local log_time
      log_time=$(date +"%Y-%m-%d %H:%M:%S")
      local message="[$log_time] $1"
      echo "$message" >> "$logfile"
      echo "$message"
  }

  # Initialize log file
  echo "### Scomatic Pipeline Log ###" > "$logfile"
  log_message "Running scomatic pipeline for run: ${run}"


  # Step 1: Splitting alignment file in cell type specific bams
  log_message "Step 1: Splitting alignment file in cell type specific bams..."

  output_dir1=$output_dir/Step1_BamCellTypes
  mkdir -p $output_dir1
  meta_dir="/media/storage/mcGinn_2021/scomatic/markers"

  # Iterate through each BAM file in the bam_dir
  for bam_file in "$bam_dir"/*.bam
  do
      # Extract the sample name from the file name
      base_name=$(basename "$bam_file" Aligned.sortedByCoord.out.bam) # SLX-17937_SIGAD8
      sample=${base_name#*_} # remove the prefix up to '_', so that we get only SIGAD8

      echo "Started sample: $sample"

      meta_file="$meta_dir/${run}/esoph_markers_scomatic_${base_name}.tsv"
      
      # Run the Python script with the appropriate arguments
      python $SCOMATIC/scripts/SplitBam/SplitBamCellTypes.py --bam "$bam_file" \
          --meta "$meta_file" \
          --id "$sample" \
          --n_trim 5 \
          --max_nM 5 \
          --max_NH 1 \
          --outdir "$output_dir1"

      log_message "Finished sample: $sample"

      sample_list+=("$sample")
  done

  # Now, we have each .bam file for each sample located in =$output_dir/Step1_BamCellTypes. In total, 16 .tsv
  # The bam file will be like: SLX-17937_SIGAD8.SIGAD8_Epi.bam


  #######################################
  # Step 2: Collecting base count information
  
  REF=/media/storage/reference_genomes/Mus_musculus.GRCm38.dna.primary_assembly.fa

  output_dir2=$output_dir/Step2_BaseCellCounts
  mkdir -p $output_dir2

  for bam in "$output_dir1"/*.bam
  do
    
    # Cell type, which will be SIGAD8_Epi in this case
    cell_type=$(basename $bam | awk -F'.' '{print $(NF-1)}') # this takes the cell_type from the bam name (the second field to last by separating by '.')

    # Temp folder
    temp=$output_dir2/temp_${cell_type}
    mkdir -p $temp

    log_message "Processing base counts for cell type: $cell_type"

    # Command line to submit to cluster
    python $SCOMATIC/scripts/BaseCellCounter/BaseCellCounter.py --bam $bam \
      --ref $REF \
      --chrom all \
      --out_folder $output_dir2 \
      --min_bq 30 \
      --tmp_dir $temp \
      --nprocs 10

    rm -rf $temp

    log_message "Finished base count processing for cell type: $cell_type"
  done

  log_message "Listing files in $output_dir2 before renaming:"
  ls "$output_dir2"
  # Now, we have each .tsv for each sample in output_dir2=$output_dir/Step2_BaseCellCounts
  # We are going to move both .tsv from the same sample to a new folder inside output_dir2, named after the sample name.

  # echo "Organizing .tsv files into sample-specific folders..."
  log_message "Renaming .tsv files..."


  for tsv in "$output_dir2"/*.tsv
  do
      # Create a new directory for the sample
      # Move the .tsv files corresponding to the sample into the sample directory
 
      base_name=$(basename "$tsv")
      
      # SIGAD8.SIGAD8_Epi.tsv ; we want SIGAD8_Epi.tsv -> this will be our unique "cell types"
    
      # Extract run name, and spec cell
      sample_run=$(echo "$base_name" | awk -F'.' '{print $1}') # extracts SIGAD8
      spec_cell=$(echo "$base_name" | awk -F'_' '{print $2}') # extracts Epi.tsv

      new_filename="${output_dir2}/${sample_run}_${spec_cell}" # SIGAD8_Epi.tsv

      # Change each name of the .tsv files by adding its run
      # mv "$output_dir2/${sample_name}"*.tsv "$sample_dir/${run_name}${sample_name}"*.tsv
      mv "$tsv" "$new_filename"

      # echo "Moved .tsv files for sample: $sample to $sample_dir"
      log_message "Renamed $tsv to $new_filename"
  done

  log_message "Listing files in $output_dir2 after renaming:"
  ls "$output_dir2"

  # echo "Moved the files successfully!"
  log_message "Renamed the files successfully!"

  #######################################
  # Step 3: Merging base count matrices
  log_message "Step 3: Merging base count matrices..."

  output_dir3=$output_dir/Step3_BaseCellCountsMerged
  mkdir -p $output_dir3

  log_message "Merging base count matrices into a single .tsv file..."

  python $SCOMATIC/scripts/MergeCounts/MergeBaseCellCounts.py \
    --tsv_folder ${output_dir2} \
    --outfile ${output_dir3}/${run}.BaseCellCounts.AllCellTypes.tsv


  #######################################
  # Step 4: Detection of somatic mutations
  log_message "Step 4: Detection of somatic mutations..."

  # Step 4.1
  log_message "Step 4.1: Variant calling..."
  output_dir4=$output_dir/Step4_VariantCalling
  mkdir -p $output_dir4

  python $SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step1.py \
            --infile ${output_dir3}/${run}.BaseCellCounts.AllCellTypes.tsv \
            --outfile ${output_dir4}/${run} \
	    --max_cell_types 1 \
	    --min_cell_types 2 \
            --ref $REF


  # Step 4.2
  log_message "Step 4.2: Somatic mutation detection..."
  editing=$SCOMATIC/RNAediting/AllEditingSites.mm10.txt
  PON=$SCOMATIC/PoNs/PoN.scRNAseq.mm10.tsv

  python $SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step2.py \
            --infile ${output_dir4}/${run}.calling.step1.tsv \
            --outfile ${output_dir4}/${run}.calling.step2.tsv \
            --editing $editing \
            --pon $PON

  # extra step: Intersection with bed file
  log_message "Extra step: Intersection with bed file..."
  bedtools intersect -header -a ${output_dir4}/${run}.calling.step2.tsv -b $SCOMATIC/bed_files_of_interest/UCSC.mm10.without.repeatmasker.bed | awk '$1 ~ /^#/ || $6 == "PASS"' > ${output_dir4}/${run}.calling.step2.pass.tsv

  log_message "Finished scomatic pipeline for ${run} successfully!"

done

log_message "Pipeline completed successfully for both runs!"

# Display log file
echo "Log file saved to: $logfile"
cat "$logfile"


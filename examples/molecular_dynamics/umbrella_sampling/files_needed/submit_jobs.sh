#!/bin/bash

main_folder="/home/users/w/wega/MD_BPe_ionpair/MD/umbrella"

for subfolder in "$main_folder"/umbrella_*; do

    if [ -d "$subfolder" ]; then

        job_file="$subfolder/umbrella_job.sh"

        if [ -f "$job_file" ]; then

            echo "Submitting: $subfolder"

            (
                cd "$subfolder" || exit 1
                sbatch ./umbrella_job.sh
            )

        else
            echo "No umbrella_job.sh found in $subfolder"
        fi

    fi

done

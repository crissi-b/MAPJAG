#!/bin/bash
#SBATCH -n 70
#SBATCH -N 1
#SBATCH --mem 299G
#SBATCH --time 48:0:0
#SBATCH --qos castles
#SBATCH --mail-type ALL
#SBATCH --account=croftap-stia-atac


export MRO_DISK_SPACE_CHECK=disable

set -e
module purge; module load bluebear

cd /rds/projects/c/croftap-celldive01/xenium/analysis/tx_counts_FOV/
for f in *.csv; do /rds/projects/c/croftap-celldive01/CosmX_analysis/baysor0.6.0/baysor/bin/baysor run -s 10 -x x_location -y y_location -g feature_name -m 10 --save-polygons=GeoJSON --prior-segmentation-confidence=0.7 --n-clusters=10 -o /rds/projects/c/croftap-celldive01/xenium/analysis/baysor_outs/"$f" "$f" :cell_id_new; done


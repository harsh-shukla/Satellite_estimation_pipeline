## LOAD MODULES
## module load python/3.10.2
####

HIFI_DATA=/dfs7/jje/hgshukla/HiFi_Subreads/ISO1/CCS_Gen_6_2/QV_0_99/ISO1_hifi_099.fastq.gz


python3 trf_parser.Disk.py -i $HIFI_DATA -o satellite.bed -e /data/homezvol1/hgshukla/Softwares/TRF-mod/trf-mod -t 32 --trf_options "-p 15 -s 100"

#python3 resolve_overlap_rearrange.py -i satellite.bed -o satellite.FILTNORM.bed --min-alignment-length 50

#python3 summarize_trf_content.py -i satellite.FILTNORM.bed -o satellite_summary.tsv --total-bp-sequenced satellite.bed.total_bp

#python3 generate_html_report.py -i satellite_summary.tsv -o satellite_report.html


# UNLOAD MODULES
# module unload python/3.10.2

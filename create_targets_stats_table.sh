K_mer_len=35
allowed_mismatches=2
genome_reference=resources/hg19_ion/hg19.fasta
genome_nickname=hg19_ion

targets=resources/AmpliseqExome/AmpliSeqExome.20141113.designed.bed
targets_nickname=ampliseq_exome

cov_files_regex="resources/cov.xls/*OKQC*"


[ ! -d "results/${genome_nickname}-index" ] && echo "Creating index for ${genome_nickname}; fasta=${genome_reference}"
[ ! -d "results/${genome_nickname}-index" ] && genmap index -F "${genome_reference}" -I "results/${genome_nickname}-index"


map_output=results/"${genome_nickname}"-K"${K_mer_len}"-E"${allowed_mismatches}"
[ ! -f "${map_output}.bedgraph" ] && echo "Computing mappability; index=results/${genome_nickname}-index"
[ ! -f "${map_output}.bedgraph" ] && genmap map -K "${K_mer_len}" -E "${allowed_mismatches}" \
-I "results/${genome_nickname}-index" \
-S "${targets}" \
-O "${map_output}" -bg


overall_map_output="results/${targets_nickname}-overall-mappability.bed"
[ ! -f "${overall_map_output}" ] && echo "Computing overall mappability; bedgraph=${map_output}.bedgraph"
[ ! -f "${overall_map_output}" ] && python scripts/overall_mappability.py \
--targets "${targets}" \
--mappability "${map_output}.bedgraph" \
--k-mer-length "${K_mer_len}" \
--output "${overall_map_output}"


stats_table_output="results/${targets_nickname}-stats-table.tsv"
[ ! -f "${stats_table_output}" ] && echo "Computing stats for ${targets_nickname}; cov_files=${map_output}.bedgraph"
[ ! -f "${stats_table_output}" ] && python scripts/target_stats_table.py \
--overall-mappability "${overall_map_output}" \
--cov-files-regex "${cov_files_regex}" \
--output "${stats_table_output}"

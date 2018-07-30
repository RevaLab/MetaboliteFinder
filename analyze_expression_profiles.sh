#!/bin/bash
#module load python/3.6.2
#module load py_packages/3.6
#
#source /sc/orga/projects/Signatures/GeneExpressionVulnerability/venv3.6/bin/activate

INPUT_FILE=$1
BASE_DIR="/Users/calina01/PycharmProjects/MetaboliteFinder"

if [ $# -eq 0 ]
  then
    echo "USAGE: sh analyze_expression_profiles.sh path/to/input"
    exit
fi

if [ ! -f ${INPUT_FILE} ]
then
    echo "Could not find $INPUT_FILE."
    exit
fi

read -p 'Database to use (k for KEGG, r for Reactome, s for HMDB/SMPDB): ' DB_CHOICE
read -p 'Enrichment (p for positive or n for negative): ' ENRICHMENT
read -p 'P-value threshold: ' P_THRESHOLD
#read -p 'Create summary with samples (y for yes, n for no): ' SHOW_SAMPLES
#read -p 'Create summary with metabolites (y for yes, n for no): ' SHOW_METABOLITES
read -p 'Genes of interest (u for upregulated, d for downregulated): ' GENE_EXPRESSION_CHOICE
read -p 'Percent cutoff (x/100): ' PERCENT_CUTOFF


if [ ${DB_CHOICE} == 'r' ] || [ ${DB_CHOICE} == 'R' ]; then
    DATABASE='reactome'
elif [ ${DB_CHOICE} == 'k' ] || [ ${DB_CHOICE} == 'K' ]; then
    DATABASE='KEGG'
elif [ ${DB_CHOICE} == 's' ] || [ ${DB_CHOICE} == 'S' ]; then
    DATABASE='HMDB_SMPDB'
else
    DATABASE='KEGG'
    echo 'Database choice not recognized. Defaulting to KEGG.'
fi

if [ ${GENE_EXPRESSION_CHOICE} == 'd' ] || [ ${GENE_EXPRESSION_CHOICE} == 'D' ]; then
    ASCENDING='y'
elif [ ${GENE_EXPRESSION_CHOICE} == 'u' ] || [ ${GENE_EXPRESSION_CHOICE} == 'U' ]; then
    ASCENDING='n'
else
    ASCENDING='y'
    echo 'Gene expression choice not recognized. Defaulting to downregulated.'
fi

## RUN SCRIPT
#read -p 'Run script as background process (y for yes, n for no): ' IS_BACKGROUND
#
#if [ ${IS_BACKGROUND} == 'y' ] || [ ${IS_BACKGROUND} == 'Y' ]; then
#    nohup python run_metabolite_analysis.py ${INPUT_FILE} ${DATABASE} ${ENRICHMENT} ${P_THRESHOLD} ${PERCENT_CUTOFF} ${ASCENDING} > GEV.out 2> GEV.err &
#
#    echo "Your program is running. You may close this screen and return later."
#elif [ ${IS_BACKGROUND} == 'n' ] || [ ${IS_BACKGROUND} == 'N' ]; then
#    python run_metabolite_analysis.py ${INPUT_FILE} ${DATABASE} ${ENRICHMENT} ${P_THRESHOLD} ${ASCENDING} ${PERCENT_CUTOFF}
#else
#    echo "Defaulting to foreground."
#    python run_metabolite_analysis.py ${INPUT_FILE} ${DATABASE} ${ENRICHMENT} ${P_THRESHOLD} ${ASCENDING} ${PERCENT_CUTOFF}
#fi

echo "#!/bin/bash
#module load python/3.6.2
#module load py_packages/3.6

python ${BASE_DIR}/run_metabolite_analysis.py ${INPUT_FILE} ${DATABASE} ${ENRICHMENT} ${P_THRESHOLD} ${ASCENDING} ${PERCENT_CUTOFF} ${BASE_DIR}" > "${BASE_DIR}/run_job.sh"

sh "${BASE_DIR}/run_job.sh"
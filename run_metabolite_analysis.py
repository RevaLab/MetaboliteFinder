import sys
import time

import numpy as np

from tumor import Tumor
from sample import Sample


t0 = time.time()
input_file = sys.argv[1]
db_name = sys.argv[2]
enrichment = sys.argv[3]
p_threshold = sys.argv[4]
# show_samples = sys.argv[5]
# show_metabolites = sys.argv[6]
ascending = sys.argv[5]
percent_cutoff = sys.argv[6]
base_dir = sys.argv[7]

exclude_unique_pw = False

tumor = Tumor(input_file,
              db_name,
              enrichment,
              p_threshold,
              # show_samples,
              # show_metabolites,
              ascending,
              percent_cutoff,
              base_dir)

sample_count = 0
all_sample_count = len(tumor.gene_expression_table.columns)

for sample_id in tumor.gene_expression_table:
    sample_count += 1

    sample = Sample(sample_id,
                    tumor.output_dir,
                    tumor.gene_expression_table[sample_id],
                    tumor.db_name,
                    tumor.db,
                    tumor.enrichment,
                    tumor.p_threshold,
                    tumor.ascending)

    percents = np.arange(tumor.percent_start, tumor.percent_cutoff, 0.1)

    for percent in percents:
        # t0_p_s = time.time()

        percent = round(percent, 1)
        print("Analyzing: {}% -- sample {}/{} -- {}".format(percent, sample_count, all_sample_count, tumor.input_file))
        percent_genes = sample.find_percent_genes(percent)

        sample.perform_enrichment_analysis_for_percent_genes(percent_genes, percent, exclude_unique_pw)
        tumor.add_enriched_pathways_to_final_summary(sample.percent_enrichments[percent],
                                                     percent)

        tumor.add_percent_genes_to_gene_summary(percent,
                                                percent_genes)
        # t1_p_s = time.time()
        # print("Time per sample/percent: {}".format(t1_p_s - t0_p_s))
    sample.write_percent_enrichment_table()

tumor.write_final_summary_table()
tumor.write_gene_summary()
tumor.write_metabolite_summary()

t1 = time.time()
print("Time elapsed: {}".format(t1 - t0))
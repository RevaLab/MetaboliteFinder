import csv
import os
import pickle

import pandas as pd
import numpy as np


class Tumor:
    def __init__(self,
                 input_file,
                 db_name,
                 enrichment,
                 p_threshold,
                 # show_samples,
                 # show_metabolites,
                 ascending,
                 percent_cutoff,
                 base_dir):

        yes_or_no = {'y': True, 'n': False}
        """
            pathway: {percent: [samples]}
        """
        self.final_summary = {}

        """
            percent: [genes]
        """
        self.gene_summary = {}

        self.ascending = yes_or_no[ascending]
        self.input_file = input_file
        self.db_name = db_name
        self.base_dir = base_dir
        self.db = pickle.load(open('{}/databases/{}.pkl'.format(self.base_dir, self.db_name), 'rb'))
        self.enrichment = enrichment
        self.p_threshold = float(p_threshold)
        self.p_threshold_str = p_threshold
        # self.show_samples = yes_or_no[show_samples]
        # self.show_metabolites = yes_or_no[show_metabolites]
        self.percent_cutoff_str = percent_cutoff
        self.percent_cutoff = float(percent_cutoff) + 0.1

        print("Creating output directory")
        self.output_dir = self.create_output_dir()

        print("Reading input file")
        self.gene_expression_table = self.read_file()

        if self.ascending:
            self.percent_start = 0.0
        else:
            self.percent_start = 0.1

    def create_output_dir(self):
        input_file_basename = os.path.basename(self.input_file)
        basename_wo_ext = '.'.join(input_file_basename.split(".")[:-1])
        if self.ascending:
            output_dir = "{}/output_dir/{}_{}_{}_bottom_{}".format(self.base_dir,
                                                                   basename_wo_ext,
                                                                   self.p_threshold_str,
                                                                   self.db_name,
                                                                   self.percent_cutoff_str)
        else:
            output_dir = "{}/output_dir/{}_{}_{}_top_{}".format(self.base_dir,
                                                                basename_wo_ext,
                                                                self.p_threshold_str,
                                                                self.db_name,
                                                                self.percent_cutoff_str)
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)

        return output_dir

    def read_file(self):
        gene_expression_table = pd.read_csv(self.input_file, header=0, index_col=0, sep="\t")
        print("Finished reading input file")
        return gene_expression_table

    def add_enriched_pathways_to_final_summary(self, enrichment_results, percent):
        if not len(enrichment_results):
            return {}

        if percent not in self.final_summary:
            self.final_summary[percent] = {}

        for pathway in enrichment_results:
            if pathway not in self.final_summary[percent]:
                self.final_summary[percent][pathway] = []

            self.final_summary[percent][pathway].append(enrichment_results[pathway]['pval'])

    def add_percent_genes_to_gene_summary(self, percent, genes):
        if percent not in self.gene_summary:
            self.gene_summary[percent] = {}

        for gene in genes:
            if gene not in self.gene_summary[percent]:
                self.gene_summary[percent][gene] = 0

            self.gene_summary[percent][gene] += 1

    def write_final_summary_table(self):
        with open('{}/final_summary.txt'.format(self.output_dir), 'w') as tsvout:
            tsvout = csv.writer(tsvout, delimiter="\t")

            for percent in self.final_summary:
                for pathway in self.final_summary[percent]:
                    geom_mean = self.geo_mean_overflow(self.final_summary[percent][pathway])
                    row = [percent,
                           pathway,
                           self.db_name,
                           'Occurrences:',
                           len(self.final_summary[percent][pathway]),
                           'Geom_mean:',
                           geom_mean,
                           'Nmtb:',
                           len(self.db['dict'][pathway]['metabolites'])
                           ]
                    tsvout.writerow(row)
        pickle.dump(self.final_summary, open('final_summary.pkl', 'wb'))
        print('finished writing pathway summary')

    def write_gene_summary(self):
        with open('{}/gene_summary.txt'.format(self.output_dir), 'w') as tsvout:
            tsvout = csv.writer(tsvout, delimiter='\t')
            for percent in self.gene_summary:
                for gene in self.gene_summary[percent]:
                    row = [
                        percent,
                        gene,
                        self.gene_summary[percent][gene]
                    ]
                    tsvout.writerow(row)
        print('finished writing gene summary')

    def write_metabolite_summary(self):
        metabolite_summary = {}
        with open('{}/metabolite_summary.txt'.format(self.output_dir), 'w') as tsvout:
            tsvout = csv.writer(tsvout, delimiter='\t')

            for percent in self.final_summary:
                for pathway in self.final_summary[percent]:
                    for metabolite in self.db['dict'][pathway]['metabolites']:

                        if metabolite not in metabolite_summary:
                            metabolite_summary[metabolite] = {
                                'pathways': set(),
                                'samples': set()
                            }

                        metabolite_summary[metabolite]['pathways'].add(pathway)
                        updated_samples_list = metabolite_summary[metabolite]['samples'].union(
                            self.final_summary[percent][pathway]
                        )
                        metabolite_summary[metabolite]['samples'] = updated_samples_list

            for metabolite in metabolite_summary:
                row = [metabolite,
                       self.db_name,
                       'Samples:',
                       len(metabolite_summary[metabolite]['samples']),
                       'Pathways:',
                       len(metabolite_summary[metabolite]['pathways'])
                       ]
                row += list(metabolite_summary[metabolite]['pathways'])
                tsvout.writerow(row)
        print("finished writing metabolite summary")

    @staticmethod
    def geo_mean_overflow(iterable):
        a = np.log(iterable)
        return np.exp(a.sum() / len(a))

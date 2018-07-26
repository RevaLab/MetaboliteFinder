import csv

from pathway_enrichment import pathway_enrichment_analysis_pw_bg


class Sample:
    def __init__(self,
                 sample_id,
                 output_dir,
                 gene_expressions,
                 db_name,
                 db,
                 enrichment,
                 p_threshold,
                 ascending):

        self.db_name = db_name
        self.db = db
        self.sample_id = sample_id
        self.ascending = ascending
        self.output_dir = output_dir
        self.enrichment = enrichment
        self.p_threshold = p_threshold

        self.output_file = '{}/{}'.format(self.output_dir, self.sample_id)

        self.percent_enrichments = {}
        self.percent_enrichments_pw_inclusive = {}

        self.gene_expressions = gene_expressions.sort_values(ascending=self.ascending)
        self.expression_bins = self.bin_genes_by_expression()
        self.zero_genes = self.gene_expressions.loc[self.gene_expressions == 0.0].axes[0].tolist()
        self.all_genes = list(self.gene_expressions.index)
        self.all_genes_with_expression = self.gene_expressions[len(self.zero_genes):]

    def bin_genes_by_expression(self):
        expression_bins = {}

        for gene, expression in self.gene_expressions.iteritems():
            expression = round(expression, 3)
            if expression not in expression_bins:
                expression_bins[expression] = []
            expression_bins[expression].append(gene)
        return expression_bins

    def find_percent_genes(self, percent):
        percent_genes = []

        if self.ascending:
            percent_index = int(len(self.all_genes_with_expression) * (percent/100)) - 1
            expression_cutoff_at_percent_index = self.all_genes_with_expression[percent_index]
        else:
            percent_index = int(len(self.gene_expressions) * (percent/100)) - 1
            expression_cutoff_at_percent_index = self.gene_expressions[percent_index]

        if percent_index < 0:
            expression_cutoff_at_percent_index = 0.0

        for expression in self.expression_bins:
            if self.ascending and expression <= expression_cutoff_at_percent_index:
                percent_genes.extend(self.expression_bins[expression])

            if not self.ascending and expression >= expression_cutoff_at_percent_index:
                percent_genes.extend(self.expression_bins[expression])

        return percent_genes

    def perform_enrichment_analysis_for_percent_genes(self, percent_genes, percent, exclude_unique_pw):
        enriched_pathways_in_percent = pathway_enrichment_analysis_pw_bg(percent_genes,
                                                                         self.all_genes,
                                                                         self.db,
                                                                         self.enrichment,
                                                                         self.p_threshold,
                                                                         exclude_unique_pw)
        self.percent_enrichments[percent] = enriched_pathways_in_percent

    def write_percent_enrichment_table(self):
        with open('{}.txt'.format(self.output_file), 'w') as tsvout:
            tsvout = csv.writer(tsvout, delimiter="\t")

            for percent in self.percent_enrichments:
                enriched_pathways = self.percent_enrichments[percent]
                for pathway in enriched_pathways:
                    row = self.write_output_row(percent,
                                                pathway,
                                                enriched_pathways[pathway])
                    tsvout.writerow(row)
        print('finished writing: {}'.format(self.output_file))

    @staticmethod
    def write_output_row(percent, pathway, enriched_pathway):
        row = [
            'Percent:',
            percent,
            'Neff_genes:',
            enriched_pathway['input_genes'],
            'Pathway:',
            pathway,
            enriched_pathway['db'],
            'length:',
            enriched_pathway['length'],
            'overlap:',
            enriched_pathway['overlap'],
            'pval:',
            enriched_pathway['pval'],
            'enrichment:',
            enriched_pathway['enrichment'],
            'bg:',
            enriched_pathway['db_u_input'],
            'Nmtb:',
            len(enriched_pathway['metabolites'])
        ]
        return row

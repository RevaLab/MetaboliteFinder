import scipy.stats as stats


# bg = input_bg_u_pw | input_bg
def pathway_enrichment_analysis_pw_bg(input_genes, input_bg, db, enrichment, p_threshold, exclude_unique_pw):
    enriched_pathways = {}

    input_genes = set(input_genes)
    input_count = len(input_genes)
    pathways = db['dict']
    bg = set(input_bg)

    for pathway in pathways:
        if exclude_unique_pw:
            # Consider only pathway genes annotated in input bg
            pathway_genes = pathways[pathway]['genes'].intersection(bg)
        else:
            pathway_genes = pathways[pathway]['genes']
            bg = bg.union(pathway_genes)

        pathway_count = len(pathway_genes)
        bg_count = len(bg)

        overlap = input_genes.intersection(pathway_genes)
        overlap_count = len(overlap)

        if not overlap_count:
            continue

        pathway_only_count = pathway_count - overlap_count
        input_only_count = input_count - overlap_count

        bg_only_count = bg_count - pathway_count - input_count + overlap_count

        enrichment_coefficient = find_enrichment_coefficient(overlap_count,
                                                             pathway_count,
                                                             input_count,
                                                             bg_count)

        if enrichment_coefficient <= 0 and enrichment == 'p':
            continue

        if enrichment_coefficient > 0 and enrichment == 'n':
            continue

        tail = find_tail(overlap_count, pathway_count, input_count, bg_count)

        try:
            odds_ratio, pvalue = stats.fisher_exact(
                [
                    [overlap_count, pathway_only_count],
                    [input_only_count, bg_only_count]
                ],
                alternative=tail
            )
        except ValueError:
            print('Something is wrong with this 2x2 table: \n')
            print([overlap_count, pathway_only_count], [input_only_count, bg_only_count])
            break

        if pvalue > float(p_threshold):
            continue

        enriched_pathways[pathway] = {
            'input_genes': input_count,
            'db': pathways[pathway]['db'],
            'length': pathway_count,
            'overlap': overlap_count,
            'pval': pvalue,
            'enrichment': enrichment_coefficient,
            'db_u_input': bg_count,
            'metabolites': pathways[pathway]['metabolites'],
            'overlap_genes': overlap
        }

    return enriched_pathways


def find_enrichment_coefficient(overlap_count, pathway_count, input_count, bg_count):
    enrichment_threshold = (pathway_count * input_count) / bg_count

    coefficient = overlap_count / enrichment_threshold

    if coefficient < 1:
        return 1 / coefficient * -1

    return coefficient


def find_tail(overlap_count, pw_count, input_count, bg_count):
    enrichment_threshold = (pw_count * input_count) / bg_count

    if overlap_count >= enrichment_threshold:
        tail = 'greater'
    else:
        tail = 'less'

    return tail


def intersect(a, b):
    """ return the intersection of two lists """
    return list(set(a) & set(b))


def union(a, b):
    """ return the union of two lists """
    return list(set(a) | set(b))



# bg = db_u_input
# def pathway_enrichment_analysis(input_genes, db, enrichment, p_threshold):
#     enriched_pathways = {}
#
#     input_count = len(input_genes)
#     db_all = db['all']
#     pathways = db['dict']
#     db_u_input = union(db_all, input_genes)
#     bg_count = len(db_u_input)
#
#     for pathway in pathways:
#         pathway_genes = pathways[pathway]['genes']
#         pathway_count = len(pathway_genes)
#
#         overlap = intersect(input_genes, pathway_genes)
#         overlap_count = len(overlap)
#         if not overlap_count:
#             continue
#
#         pathway_only_count = pathway_count - overlap_count
#         input_only_count = input_count - overlap_count
#
#         bg_only_count = bg_count - pathway_count - input_count + overlap_count
#
#         enrichment_coefficient = find_enrichment_coefficient(overlap_count,
#                                                              pathway_count,
#                                                              input_count,
#                                                              bg_count)
#
#         if enrichment_coefficient <= 0 and enrichment == 'p':
#             continue
#
#         if enrichment_coefficient > 0 and enrichment == 'n':
#             continue
#
#         tail = find_tail(overlap_count, pathway_count, input_count, bg_count)
#
#         try:
#             odds_ratio, pvalue = stats.fisher_exact(
#                 [
#                     [overlap_count, pathway_only_count],
#                     [input_only_count, bg_only_count]
#                 ],
#                 alternative=tail
#             )
#         except ValueError:
#             print('Something is wrong with this 2x2 table: \n')
#             print([overlap_count, pathway_only_count], [input_only_count, bg_only_count])
#             break
#
#         if pvalue > float(p_threshold):
#             continue
#
#         enriched_pathways[pathway] = {
#             'input_genes': input_count,
#             'db': pathways[pathway]['db'],
#             'length': pathway_count,
#             'overlap': overlap_count,
#             'pval': pvalue,
#             'enrichment': enrichment_coefficient,
#             'db_u_input': bg_count,
#             'metabolites': pathways[pathway]['metabolites'],
#             'overlap_genes': overlap
#         }
#
#     return enriched_pathways
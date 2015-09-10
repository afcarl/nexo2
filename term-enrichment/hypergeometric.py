from term_mapping import TermMapping
import qvalue
import numpy as np
import math
from scipy import stats

# p-value cutoff
p_val_threshold = 0.01

# At least this number of genes should be assigned
gene_threshold = 2


class HypergeometricTest:

    _mapper = None

    def __init__(self):
        self._mapper = TermMapping()

    def perform_test(self, ontology_type, genes_of_interest, cutoff):

        if cutoff is None:
            cutoff = p_val_threshold

        total_num_genes = self._mapper.get_gene_count(ontology_type)
        term2genes = self._mapper.get_term_mapping(ontology_type)
        gene2terms = self._mapper.get_gene_mapping(ontology_type)

        filtered = self.__filter_input(genes_of_interest, gene2terms)

        # print(filtered)
        # print(len(filtered))

        sample_map = self.__calculateTermFrequency(filtered, gene2terms)

        n = len(filtered)

        # Number of tests performed: will be used for correction.
        num_tests = len(sample_map)

        results = [{}] * num_tests

        pvals = np.zeros(num_tests)

        idx = 0
        for term in sample_map:
            # Calculate p-value
            sampled = sample_map[term]
            assigned_genes = term2genes[term]
            k = len(sampled)
            m = len(assigned_genes)

            p = stats.hypergeom.pmf(k,total_num_genes,m,n)
            p_corrected = p * num_tests
            pvals[idx] = p

            result = {}
            result['id'] = term
            result['p-value'] = p_corrected
            result['background'] = m
            result['genes'] = list(sampled)
            results[idx] = result
            idx = idx + 1
            # print(term + " = " + str(p) + ", k = " + str(k) + ",
            # m = " + str(m) + ", total = " + str(total_num_genes) + ", n= " + str(n))

        # Correct border values (library does not accept 0 & 1)
        i = 0
        for p in pvals:
            if p >= 1.0 or math.isnan(p):
                pvals[i] = 0.9999999999999999999999
            elif p <= 0:
                pvals[i] = 0.0000000000000000000001
            i += 1
        qvals = qvalue.estimate(pvals)
        filtered_results = []

        idx = 0

        for term in sample_map:
            qv = qvals[idx]
            res = results[idx]

            pv = res['p-value']
            k = len(res['genes'])
            res['q-value'] = qv
            if pv < cutoff and k >= gene_threshold:
                filtered_results.append(res)

            idx = idx + 1

        return {
            'results': filtered_results,
            'total_genes': total_num_genes
        }

    # remove junks & duplicates
    def __filter_input(self, input_list, gene2terms):
        filtered = set([])

        for gene in input_list:
            if gene in gene2terms:
                filtered.add(gene)

        return list(filtered)

    def __calculateTermFrequency(self, genes_of_interest, gene2terms):
        sample_term_map = {}

        for gene in genes_of_interest:
            terms = []
            if gene in gene2terms:
                terms = gene2terms[gene]
            else:
                continue

            for term in terms:
                if term in sample_term_map:
                    assigned_genes = sample_term_map[term]
                    assigned_genes.add(gene)
                    sample_term_map[term] = assigned_genes
                else:
                    assigned_genes = set([])
                    assigned_genes.add(gene)
                    sample_term_map[term] = assigned_genes
        return sample_term_map

import scipy.stats
from numpy import array, empty
import numpy as np
from genebench.datatypes import GeneDiffInput, GeneMethodResult
from genebench.diffmethods.base.diffmethod import DiffMethod

class TTest(DiffMethod):

    def run(self, input: GeneDiffInput) -> GeneMethodResult:
        genes, pvalues = self._get_pvalues(input.control,
                                           input.perturbed,
                                           input.genes)
        pvalues = self._correct_pvalues(pvalues)
        values = np.ones(len(genes)) - pvalues
        values = np.where(values >= 0.95, values, np.clip(values - 0.5, 0, 1))
        genes, values = self.sort(genes, values)
        return GeneMethodResult.from_separate_lists(genes, values)

    def _get_pvalues(self, A, B, genes):
        pvalues = []
        for i in range(len(A)):
            # equal_var=False specifies to use Welch's t-test.
            ttest_results = scipy.stats.ttest_ind(A[i], B[i], equal_var=False)
            pvalue = ttest_results.pvalue
            pvalues.append(pvalue)

        return genes, pvalues

    def _correct_pvalues(self, pvalues, correction_type="FDR"):
        pvalues = array(pvalues)
        sample_size = pvalues.shape[0]
        qvalues = empty(sample_size)
        if correction_type == "Bonferroni":
            # Bonferroni correction
            qvalues = sample_size * pvalues
        elif correction_type == "Bonferroni-Holm":
            # Bonferroni-Holm correction
            values = [(pvalue, i) for i, pvalue in enumerate(pvalues)]
            values.sort()
            for rank, vals in enumerate(values):
                pvalue, i = vals
                qvalues[i] = (sample_size-rank) * pvalue
        elif correction_type == "FDR":
            # Benjamini-Hochberg, AKA - FDR test
            values = [(pvalue, i) for i, pvalue in enumerate(pvalues)]
            values.sort()
            values.reverse()
            new_values = []
            for i, vals in enumerate(values):
                rank = sample_size - i
                pvalue, index = vals
                new_values.append((sample_size/rank) * pvalue)
            for i in range(0, int(sample_size)-1):
                if new_values[i] < new_values[i+1]:
                    new_values[i+1] = new_values[i]
            for i, vals in enumerate(values):
                pvalue, index = vals
                qvalues[index] = new_values[i]
        return qvalues

from genebench.datatypes import GeneDiffInput
from genebench.datatypes import GeneMethodResult
import numpy as np


class DiffMethod:

    def setup(self, config):
        pass

    def train(self, input: GeneDiffInput):
        pass

    def run(self, input: GeneDiffInput) -> GeneMethodResult:
        pass

    def sort(self, genes, values, is_reversed=True):
        genes = np.array(list(genes))
        values = np.array(values)
        sort_indices = np.argsort(values)
        values = values[sort_indices]
        genes = genes[sort_indices]
        genes = [str(x) for x in genes]
        if is_reversed:
            genes = list(reversed(genes))
            values = list(reversed(values))
        return genes, values

import numpy as np
from genebench.datatypes import GeneDiffInput, GeneMethodResult
from genebench.diffmethods.base.diffmethod import DiffMethod


class Random(DiffMethod):
    def run(self, input: GeneDiffInput) -> GeneMethodResult:
        rankings = np.random.rand(len(input.genes))
        genes, scores = self.sort(input.genes, rankings)
        return GeneMethodResult.from_separate_lists(genes, scores)

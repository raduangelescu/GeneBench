from datatypes import GeneDiffInput, GeneMethodResult
from diffmethods.base.diffmethod import DiffMethod
import numpy as np


class Random(DiffMethod):
    def run(self, input: GeneDiffInput) -> GeneMethodResult:
        rankings = np.random.rand(len(input.genes))
        genes, scores = self.sort(input.genes, rankings)
        return GeneMethodResult.from_separate_lists(genes, scores)

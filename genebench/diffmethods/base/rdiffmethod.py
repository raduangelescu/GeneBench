import os
import rpy2.robjects as ro
from rpy2.robjects import (pandas2ri, numpy2ri)
from rpy2.robjects.conversion import localconverter
from genebench.datatypes import GeneDiffInput, GeneMethodResult
from genebench.utils import Utils
from genebench.diffmethods.base.diffmethod import DiffMethod

class RDiffMethodConfig:
    def __init__(self,
                 file_name,
                 method_name):
        self.file_name = file_name
        self.method_name = method_name


class RDiffMethod(DiffMethod):
    def setup(self, config: RDiffMethodConfig):
        self.config = config
        self.logger = Utils.get_logger(self.config.method_name)
        r_file_name = self.config.file_name
        current_file_path = os.path.dirname(os.path.abspath(__file__))
        self.abs_r_file_path = os.path.join(current_file_path,
                                            'R',
                                            r_file_name)

    def train(self, A, B, genes):
        self.logger.info(f'no training for method {self.config.method_name}')

    def prepare_data(self, A, B, genes):
        pass

    def run_custom_method(r, data_df, mask, gene_df):
        pass

    def get_converters(self):
        return ro.default_converter + pandas2ri.converter

    def post_proces_results(self, genes, values):
        return genes, values

    def run(self, input: GeneDiffInput) -> GeneMethodResult:
        numpy2ri.activate()
        pandas2ri.activate()
        r = ro.r
        r.source(self.abs_r_file_path)
        A = input.control
        B = input.perturbed
        genes = input.genes
        data_df, mask, gene_df = self.prepare_data(A, B, genes)
        converters = self.get_converters()
        with localconverter(converters):
            self.logger.info(f'running method {self.config.method_name}')
            genes, values = self.run_custom_method(r, data_df, mask, gene_df)
            genes, values = self.post_proces_results(genes, values)
            return GeneMethodResult.from_separate_lists(genes, values)

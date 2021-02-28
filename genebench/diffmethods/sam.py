import rpy2.robjects as ro
import pandas as pd
from rpy2.robjects import (pandas2ri, numpy2ri)
import numpy as np
import scipy
from genebench.diffmethods.base.rdiffmethod import RDiffMethod
from genebench.diffmethods.base.rdiffmethod import RDiffMethodConfig

class SAM(RDiffMethod):

    def setup(self, config):
        super().setup(RDiffMethodConfig(file_name='sam.r',
                                        method_name='sam'))

    def prepare_data(self, A, B, genes):
        pdA = pd.DataFrame(A)
        pdB = pd.DataFrame(B)
        data_df = pd.concat([pdA, pdB], axis=1, ignore_index=True)
        mask = ([1.0] * pdA.shape[1])
        mask.extend([2.0] * pdB.shape[1])
        gene_df = pd.DataFrame(genes)
        mask_df = pd.DataFrame(mask).transpose()
        return [data_df, mask_df, gene_df]

    def get_converters(self):
        ro.conversion.py2ri = numpy2ri
        convert = ro.default_converter + pandas2ri.converter + numpy2ri.converter
        return convert

    def post_proces_results(self, genes, values):
        n = len(genes)
        pvalues = scipy.stats.t.sf(np.abs(values), n-1) * 2
        values = np.ones(n) - pvalues
        values = np.where(values >= 0.95, values, np.clip(values - 0.5, 0, 1))
        return self.sort(genes, values)

    def run_custom_method(self, r, data, mask, genes):
        np_data = np.array(data.to_numpy())
        np_class = np.array(mask)
        np_gene_names = np.array(genes.to_numpy())
        result = r.calculate_sam(np_data, np_class, np_gene_names)
        genes_ids = result['Gene ID']
        genes_ids = np.transpose(genes_ids)[0]
        genes_ids = genes_ids.tolist()
        scores = result['Score(d)']
        scores = np.transpose(scores)[0]
        scores = scores.tolist()
        return genes_ids, scores

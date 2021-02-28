import rpy2.robjects as ro
import pandas as pd
import numpy as np
from genebench.diffmethods.base.rdiffmethod import RDiffMethod
from genebench.diffmethods.base.rdiffmethod import RDiffMethodConfig


class LIMMA(RDiffMethod):

    def setup(self, config):
        super().setup(RDiffMethodConfig(file_name='limma.r',
                                        method_name='limma'))

    def prepare_data(self, A, B, genes):
        pdA = pd.DataFrame(A)
        pdB = pd.DataFrame(B)
        data_df = pd.concat([pdA, pdB],
                            axis=1,
                            ignore_index=True)
        mask = ([0] * pdA.shape[1])
        mask.extend([1] * pdB.shape[1])
        mask = pd.DataFrame(mask)
        return [data_df, mask, genes]

    def post_proces_results(self, genes, values):
        values = np.ones(len(genes)) - values
        return self.sort(genes, values)

    def run_custom_method(self, r, data, mask, genes):
        r_data = ro.conversion.py2rpy(data)
        r_class = ro.conversion.py2rpy(mask)
        result = r.calculate_limma(r_data, r_class)
        result_py = ro.conversion.rpy2py(result)
        p_values = list(result_py['adj.P.Val'])
        return genes, p_values

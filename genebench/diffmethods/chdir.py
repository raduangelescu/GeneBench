import numpy as np
import pandas as pd

from genebench.diffmethods.base.rdiffmethod import RDiffMethod
from genebench.diffmethods.base.rdiffmethod import RDiffMethodConfig
# Note that we don't sort results as they should be already sorted
# from the R algorithm: https://rdrr.io/cran/GeoDE/man/chdirAnalysis.html


class ChDir(RDiffMethod):
    def setup(self, config):
        config = RDiffMethodConfig(file_name='characteristic_direction.r',
                                   method_name='CharacteristicDirection')
        super().setup(config)

    def prepare_data(self, A, B, genes):
        pdA = pd.DataFrame(A)
        pdA.insert(0, "gene_names", genes)
        pdB = pd.DataFrame(B)
        data_df = pd.concat([pdA, pdB], axis=1, ignore_index=True)
        mask = ([1.0] * (pdA.shape[1]-1))
        mask.extend([2.0] * pdB.shape[1])
        mask = pd.DataFrame(mask)
        gene_df = pd.DataFrame(genes)
        return [data_df, mask, gene_df]

    def run_custom_method(self, r, data, mask, genes):
        result = r.characteristic_direction(data, mask, genes)
        genes = list(result.index.values)
        values = result.to_numpy().tolist()
        return genes, values

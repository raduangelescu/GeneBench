from numpy.lib.function_base import diff
from genebench.silico.generators.base import SilicoDataGenerator
from genebench.datatypes import GeneData
from genebench.datatypes import GeneDiffValidation
from genebench.datatypes import GeoData
from genebench.utils import Utils
import numpy as np
import random


class LinearDataGeneratorParam():
    def __init__(self,
                 diff_factor,
                 noise_factor,
                 num_replicates):
        self.diff_factor = diff_factor
        self.noise_factor = noise_factor
        self.num_replicates = num_replicates
        pass


class LinearDataGenerator(SilicoDataGenerator):
    def __init__(self, config):
        super().__init__(config)
        self.param = LinearDataGeneratorParam(**self.config.params)
        self.logger = Utils.get_logger('LinearDataGenerator')
        pass

    def generate_single(self, validation_data, id, num_genes) -> GeoData:
        all_tfs = list(validation_data.data.keys())
        picked_tf = random.choice(all_tfs)
        perturbed_genes = set(validation_data.data[picked_tf])
        genes = Utils.get_random_gene_names(num_genes)
        mask = []
        for gene in genes:
            if gene in perturbed_genes:
                mask.append(1)
            else:
                mask.append(0)
        mask = np.array(mask)
        df_factor = self.param.diff_factor
        validation = []
        for index, mask_value in enumerate(mask.tolist()):
            if mask_value == 1:
                validation.append(genes[index])
        num_replicates = self.param.num_replicates
        mask = np.array([mask])
        mask = np.repeat(mask, num_replicates, axis=0).T
        control = np.random.rand(num_genes, num_replicates)
        effect = np.random.rand(num_genes, num_replicates) * df_factor
        perturbation = control + np.multiply(mask, effect)

        gene_data = GeoData({
            "name": f"SIL_{id}",
            "perturbed_series_names": ['fakeseries'],
            "control_series_names": ['fakeseries'],
            "extra_info": {"none": "none"},
            "perturbed_array": perturbation.tolist(),
            "control_array": control.tolist(),
            "source": self.config.source,
            "genes": genes,
            "pf": picked_tf
        })

        return gene_data

    def generate_experiments(self,
                             validation_data: GeneDiffValidation,
                             num_experiments,
                             num_genes):
        experiments = []
        for id in range(0, num_experiments):
            experiments.append(self.generate_single(validation_data,
                                                    id,
                                                    num_genes))
        return experiments

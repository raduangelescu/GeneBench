from genebench.evaluationmetrics.base import Metric
from genebench.datatypes import GeneDiffValidation, GeneMethodResult
from genebench.utils import Utils
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import uniform
import os


class Kolmogorov(Metric):

    def __init__(self, config, output_folder):
        self.name = "Kolmogorov"
        self.logger = Utils.get_logger("metric_kolmogorov")
        self.method_rks = {}
        self.config = config
        self.output_folder = output_folder
        pass

    def plot_cdf(self, data):
        x, y = Utils.ecdf(data)
        x = np.append(x, [1.0])
        y = np.append(y, [1.0])
        y = y - uniform.cdf(x)
        plt.plot(x, y)
        plt.xlabel('rank', fontsize=16)
        plt.ylabel('cdf(r)-r', fontsize=16)

    def add(self,
            pf,
            method_name,
            validation: GeneDiffValidation,
            result: GeneMethodResult):

        if validation.source not in self.method_rks:
            self.method_rks[validation.source] = {}

        method_rks = self.method_rks[validation.source]
        pf = pf.lower()
        self.logger.info(f"perturbation factor: {pf}")
        # we collect all results
        if method_name not in method_rks:
            method_rks[method_name] = []

        result_data = result.result
        valid_data = validation.data
        if pf not in valid_data:
            self.logger.error(f"{pf} not found in set {validation.source}")
            return

        genes = valid_data[pf]
        valid_genes = set([x.lower() for x in genes])

        rks = []
        number_of_genes = len(result_data)
        for index, gene_entry in enumerate(result_data):
            gene_name = gene_entry.gene_name.lower()
            if gene_name in valid_genes:
                rank = index / number_of_genes
                rks.append(rank)

        method_rks[method_name].extend(rks)

    def evaluate(self, group_name):
        fig = plt.figure()
        Utils.create_folder_if_not_exist(self.output_folder)
        save_path = os.path.join(self.output_folder, self.name)
        Utils.create_folder_if_not_exist(save_path)

        for validation_source in self.method_rks.keys():
            method_rks = self.method_rks[validation_source]
            methods = []
            for method_name, rks in method_rks.items():
                rks_array = np.sort(np.array(rks))
                self.plot_cdf(rks_array)
                methods.append(method_name)
            name = f"{group_name}_{validation_source}"
            path_file = f"{name}.png"
            path_method = os.path.join(save_path,
                                       path_file)
            plt.legend(methods)
            plt.savefig(path_method)
            plt.title(f"{name}")
            plt.clf()
            fig.clf()
        self.method_rks = {}
        self.logger.info("done")

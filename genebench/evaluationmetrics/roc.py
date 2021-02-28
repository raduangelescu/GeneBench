from genebench.evaluationmetrics.base import Metric
from genebench.datatypes import GeneDiffValidation, GeneMethodResult
from genebench.utils import Utils
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import uniform
from sklearn.metrics import roc_curve, roc_auc_score
import os


class ROC(Metric):

    def __init__(self, config, output_folder):
        self.name = "ROC"
        self.logger = Utils.get_logger("metric_roc")
        self.config = config
        self.output_folder = output_folder
        self.method_roc = {}
        pass

    def add(self,
            pf,
            method_name,
            validation: GeneDiffValidation,
            result: GeneMethodResult):

        if validation.source not in self.method_roc:
            self.method_roc[validation.source] = {}

        method_roc = self.method_roc[validation.source]
        pf = pf.lower()
        self.logger.info(f"perturbation factor: {pf}")
        # we collect all results
        if method_name not in method_roc:
            method_roc[method_name] = {'y': [], 'pred': []}

        result_data = result.result
        valid_data = validation.data
        if pf not in valid_data:
            self.logger.error(f"{pf} not found in set {validation.source}")
            return

        genes = valid_data[pf]
        valid_genes = set([x.lower() for x in genes])
        real_class = []
        pred = []
        for gene_entry in result_data:
            gene_name = gene_entry.gene_name.lower()
            if gene_name in valid_genes:
                real_class.append(1)
            else:
                real_class.append(0)
            pred.append(gene_entry.score)
        method_roc[method_name]['y'].extend(real_class)
        method_roc[method_name]['pred'].extend(pred)

    def evaluate(self, group_name):
        fig = plt.figure()
        Utils.create_folder_if_not_exist(self.output_folder)
        save_path = os.path.join(self.output_folder, self.name)
        Utils.create_folder_if_not_exist(save_path)

        for validation_source in self.method_roc.keys():
            method_roc = self.method_roc[validation_source]
            methods = []
            for method_name, roc in method_roc.items():
                pred = roc['pred']
                pred = np.nan_to_num(pred, True, 0.0, 1.0, 0.0)
                fpr, tpr, _ = roc_curve(roc['y'], pred)
                score = roc_auc_score(roc['y'], pred)
                tpr = tpr - uniform.cdf(fpr)
                plt.plot(fpr, tpr)
                plt.xlabel('False Positive Rate', fontsize=16)
                plt.ylabel('True Positive Rate', fontsize=16)
                methods.append(f"{method_name} AUC {score:.4f}")
            name = f"{group_name}_{validation_source}"
            path_file = f"{name}.png"
            path_method = os.path.join(save_path,
                                       path_file)
            plt.legend(methods)
            plt.savefig(path_method)
            plt.title(f"{name}")
            plt.clf()
            fig.clf()
        self.method_roc = {}
        self.logger.info("done")

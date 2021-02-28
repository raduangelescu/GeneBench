from genebench.evaluationmetrics.base import Metric
from genebench.datatypes import GeneDiffValidation, GeneMethodResult
from genebench.utils import Utils
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import uniform
from sklearn.metrics import f1_score
import os


class F1(Metric):

    def __init__(self, config, output_folder):
        self.name = "F1"
        self.logger = Utils.get_logger("metric_f1")
        self.config = config
        self.output_folder = output_folder
        self.method_f1 = {}
        pass

    def add(self,
            pf,
            method_name,
            validation: GeneDiffValidation,
            result: GeneMethodResult):

        if validation.source not in self.method_f1:
            self.method_f1[validation.source] = {}

        method_f1 = self.method_f1[validation.source]
        pf = pf.lower()
        self.logger.info(f"perturbation factor: {pf}")
        # we collect all results
        if method_name not in method_f1:
            method_f1[method_name] = {'y': [], 'pred': []}

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
            if gene_entry.score >= 0.5:
                pred.append(1)
            else:
                pred.append(0)
        method_f1[method_name]['y'].extend(real_class)
        method_f1[method_name]['pred'].extend(pred)

    def evaluate(self, group_name):
        Utils.create_folder_if_not_exist(self.output_folder)
        save_path = os.path.join(self.output_folder, self.name)
        Utils.create_folder_if_not_exist(save_path)

        for validation_source in self.method_f1.keys():
            method_f1 = self.method_f1[validation_source]
            methods = []
            for method_name, _f1 in method_f1.items():
                pred = _f1['pred']
                pred = np.nan_to_num(pred, True, 0.0, 1.0, 0.0)
                f1 = f1_score(_f1['y'], pred)
                methods.append(f"{method_name} F1 Score: {f1:.4f}")
            name = f"{group_name}_{validation_source}"
            path_file = f"{name}.txt"
            path_method = os.path.join(save_path,
                                       path_file)
            with open(path_method, mode='wt', encoding='utf-8') as out_scores:
                out_scores.write('\n'.join(methods))
        self.method_f1 = {}
        self.logger.info("done")

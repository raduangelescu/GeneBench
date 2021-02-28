import numpy as np
import logging
import json
import matplotlib.pyplot as plt
from scipy.stats import uniform
import os
import pandas as pd


class Utils:
    POS_INF = 10000000
    NEG_INF = 0

    @staticmethod
    def quantile_normalize(df_input):
        df = df_input.copy()
        dic = {}
        for col in df:
            dic.update({col: sorted(df[col])})
        sorted_df = pd.DataFrame(dic)
        rank = sorted_df.mean(axis=1).tolist()
        for col in df:
            t = np.searchsorted(np.sort(df[col]), df[col])
            df[col] = [rank[i] for i in t]
        return df

    @staticmethod
    def log_if_necessary(data):
        if np.max(data) - np.min(data) > 100:
            data = np.where(data == 0, 1, data)
            return np.log2(data)
        return data

    @staticmethod
    def apply_cutoff(ranked_genes, cutoff):
        """Applies a cutoff to both lists, assuming left-to-right
        least-to-greatest.
        """
        if cutoff is None:
            return ranked_genes
        return ranked_genes[:cutoff]

    @staticmethod
    def get_random_gene_names(num_genes):
        gene_names = [f"g_{x}" for x in range(0, num_genes)]
        return gene_names

    @staticmethod
    def get_random_tf_names(num_perturbation_factors):
        tf_names = [f"tf_{x}" for x in range(0, num_perturbation_factors)]
        return tf_names

    class SetEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, set):
                return list(obj)
            return json.JSONEncoder.default(self, obj)

    @staticmethod
    def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
        return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

    @staticmethod
    def delete_files_from_folder(pth):
        files = Utils.list_files_in_folder(pth)
        for file in files:
            os.unlink(file)

    @staticmethod
    def list_files_in_folder(pth):
        files = [os.path.join(pth, f) for f in os.listdir(pth) if os.path.isfile(os.path.join(pth, f))]
        return sorted(files)

    @staticmethod
    def list_filenames_in_folder(pth):
        files = [ f for f in os.listdir(pth) if os.path.isfile(os.path.join(pth, f))]
        return sorted(files)

    @staticmethod
    def list_folders_in_folder(pth):
        files = [f for f in os.listdir(pth) if os.path.isdir(os.path.join(pth, f))]
        return sorted(files)

    @staticmethod
    def create_folder_if_not_exist(folder_path):
        try:
            os.mkdir(folder_path)
        except OSError:
            return True
        else:
            return True

    @staticmethod
    def get_config(config_filename, section):
        with open(config_filename) as f:
            config_object = json.load(f)
            return config_object[section]

    @staticmethod
    def get_logger(name):
        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger(name)
        logs_folder = os.path.join('..', 'logs')
        Utils.create_folder_if_not_exist(logs_folder)
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        fh = logging.FileHandler(os.path.join(logs_folder, f"{name}.log"))
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        fh.setFormatter(formatter)
        ch.setFormatter(formatter)
        logger.addHandler(fh)
        logger.addHandler(ch)
        return logger

    @staticmethod
    def find_in_array(df, array):
        for type_label in array:
            if type_label in df:
                return type_label
        return 'unknown'

    def repair_nan_fast(nparray):
        mean = np.nanmean(nparray, axis=1)
        mean = np.nan_to_num(mean, copy=True, nan=0,
                             posinf=Utils.POS_INF,
                             neginf=Utils.NEG_INF)
        inds = np.where(np.isnan(nparray))
        nparray[inds] = np.take(mean, inds[0])
        return nparray

    @staticmethod
    def repair_nan(df):
        df_mean = df.mean(axis=1)
        df_mean_fill = df.T.fillna(df_mean).T
        df_mean_fill_zero = df_mean_fill.fillna(value=0.0)
        return df_mean_fill_zero

    @staticmethod
    def is_control(df, control_labels):
        no_find_exact = Utils.find_in_array(df, control_labels) != 'unknown'
        return no_find_exact

    @staticmethod
    def deduplicate_genes(genes_read):
        gene_set = set()
        unique_id = 99
        gene_list = []
        for gene in genes_read:
            if gene in gene_set:
                unique_id = unique_id + 1
                gene_set.add(f'{gene}_{unique_id}')
                gene_list.append(f'{gene}_{unique_id}')
            else:
                gene_set.add(gene)
                gene_list.append(gene)
        return gene_list

    @staticmethod
    def ecdf(data):
        x = np.sort(data)
        n = x.size
        y = np.arange(1, n+1) / n
        return(x, y)

    @staticmethod
    def plot_cdf_test(test):
        x, y = Utils.ecdf(test)
        x = np.append(x, [1.0])
        y = np.append(y, [1.0])
        y = y - uniform.cdf(x)
        plt.plot(x, y)
        plt.xlabel('rank', fontsize=16)
        plt.ylabel('cdf(r)-r', fontsize=16)

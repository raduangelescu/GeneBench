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
    def get_feature_vector_size(self):
        return 11

    @staticmethod
    def get_feature_vector(control_gene, perturbed_gene):
        # y_g = (y_g1, y_g2... y_gn) - column vector  whene n is the number of RNA samples
        # X = nXp design matrix representing experimental design
        # Beta_g = unknown coefficient vector that parametrizes the average expression levels in each experimental condition
        # y_g_i is assumed independant with wariance sigma_g^2
        c_g = control_gene
        p_g = perturbed_gene
        # number of same gene values (control + perturbation)
        n = c_g.shape[0] + p_g.shape[0]
        # number of situations
        p = 2

        X = np.append(np.zeros((c_g.shape[0], 1)),
                      np.ones((p_g.shape[0], 1)),
                      axis=0)
        yg = np.append(c_g, p_g)
        inverse_mtx = np.linalg.inv(np.matmul(np.transpose(X), X))
        inverse_transpose_mul = np.matmul(inverse_mtx, np.transpose(X))
        beta_g = np.matmul(inverse_transpose_mul, yg)
        dg = n - p
        mu = np.matmul(X, beta_g)
        residual = yg - mu
        sg2 = np.dot(residual, residual)/dg
        feature_vector = [beta_g[0],
                          sg2,
                          dg,
                          np.var(c_g),
                          np.var(p_g),
                          c_g[0],
                          c_g[1],
                          c_g[2],
                          p_g[0],
                          p_g[1],
                          p_g[2]]
        return np.array(feature_vector)

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

    @staticmethod
    def filter_data(logger, data):
        np_data = np.array(data)
        if np.isnan(np_data):
            logger.warning("Bad data, we need to fix NAN and Inf")
        np_data = np.nan_to_num(np_data,
                                nan=0.0,
                                posinf=99999.0,
                                neginf=-99999.0)
        np_data = log_if_necessary(np_data.T)

        if np.isnan(np_data).any():
            logger.error("Bad data, not log")
            return False

        pd_data = pd.DataFrame(np_data)
        pd_data_q = quantileNormalize(pd_data)
        if np.isnan(pd_data_q.to_numpy()).any():
            logger.error("Bad data, bad normalization")
            return False
        return pd_data_q

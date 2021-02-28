from genebench.utils import Utils
import numpy as np
import pickle
from genebench.diffmethods.base.diffmethod import DiffMethod


class MIDGET(DiffMethod):
    def __do_validation_source(self, data, source, geo_id_tf, tf_genes):
        for metadata in data:
            tf_name = metadata['tf'].lower()
            geo_id_tf[metadata['geoid']] = tf_name
            try:
                validation_data = self.storage.get_validation_data(source,
                                                                   tf_name)[0]
            except Exception as ex:
                msg = f'tf: {tf_name} src: {source} error: {str(ex)}'
                self.logger.warning(msg)
                continue
            if tf_name not in tf_genes:
                tf_genes[tf_name] = set()
            for valid_gene in validation_data['genes']:
                lower_name = valid_gene.lower()
                tf_genes[tf_name].add(lower_name)

    def filter_data(self, logger, data):
        np_data = np.array(data)
        if np.isnan(np_data):
            logger.warning("Bad data, we need to fix NAN and Inf")
        np_data = np.nan_to_num(np_data,
                                nan=0.0,
                                posinf=99999.0,
                                neginf=-99999.0)
        np_data = Utils.log_if_necessary(np_data.T)

        if np.isnan(np_data).any():
            logger.error("Bad data, not log")
            return False

        pd_data = pd.DataFrame(np_data)
        pd_data_q = Utils.quantile_normalize(pd_data)
        if np.isnan(pd_data_q.to_numpy()).any():
            logger.error("Bad data, bad normalization")
            return False
        return pd_data_q

    def __get_all_raw_geo(self):
        all_data = []
        geo_datas = self.storage.get_geo_data()
        for data in geo_datas:
            geo_id = data["name"]
            valid_c, control = Utils.filter_data(self.logger,
                                                 data["control"])
            valid_p, perturbed = Utils.filter_data(self.logger,
                                                   data["perturbed"])

            if valid_c is False or valid_p is False:
                self.logger.error(f"Bad data, skiping geo_id {geo_id}")
                continue

            gene_names = data["genes"]
            gene_names = [name.lower() for name in gene_names]
            all_data.append({
                'control': control,
                'perturbed': perturbed,
                'genes': gene_names,
                'geo_id': geo_id})

    def get_validation_data(self):
        tf_genes = {}
        geo_id_tf = {}
        db_valid = self.storage.get_validation_sources()
        validation_sources = [x for x in db_valid]
        data = self.storage.get_geo_tf_data()[0]['data']
        for validation_source in validation_sources:
            self.logger.info(f'Using validation source {validation_source}')
            self.__do_validation_source(data,
                                        validation_source,
                                        geo_id_tf,
                                        tf_genes)
        return {'tf_genes': tf_genes, 'geo_id_tf': geo_id_tf}

    def get_feature_vector(self, control_gene, perturbed_gene):
        # y_g = (y_g1, y_g2... y_gn) - column vector  whene n is
        #  the number of RNA samples
        # X = nXp design matrix representing experimental design
        # Beta_g = unknown coefficient vector that parametrizes
        # the average expression levels in each experimental condition
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

    def save_training_data(self, file_name):
        validation_data = self.get_validation_data()
        self.logger.info('training with data from storage')
        raw_data = self.__get_all_raw_geo()

        with open(file_name, 'wb') as f:
            pickle.dump([validation_data, raw_data], f)

    def load_training(self, file_name):
        self.logger.info(f'loading training data from {file_name}')
        with open(file_name, 'rb') as f:
            return pickle.load(f)

    def get_data(self, data):
        X = []
        y = []
        data_len = len(data)
        for row_idx in range(0, data_len):
            data_row = data[row_idx]
            X.append(data_row['x'])
            y.append(data_row['y'])
        return [np.array(X), np.array(y)]

    def do_run_feature_vectors(self, control, perturbed, genes):
        retX = []
        for row_idx in range(0, len(genes)):
            gene_control = control[row_idx]
            gene_perturbed = perturbed[row_idx]
            feature_vector = self.get_feature_vector(gene_control,
                                                     gene_perturbed)
            retX.append(feature_vector)
        retX = np.array(retX)
        return retX

    def save_all_feature_vectors(self, file_name_training, out_file_name):
        validation_data, raw_data = self.load_training(file_name_training)
        all_feature_vectors = []
        num_geos = len(raw_data)
        for id, data in enumerate(raw_data):
            experiment_name = data['geo_id']
            experiment_genes = data['genes']
            id_tf = validation_data['geo_id_tf']
            tf_genes = validation_data['tf_genes']
            control = data['control']
            perturbed = data['perturbed']

            progress_msg = f"geo {experiment_name} [{id}/{num_geos}]"
            self.logger.info(progress_msg)
            tf_for_experiment = id_tf[experiment_name]
            if tf_for_experiment not in tf_genes:
                continue

            valid_genes = tf_genes[tf_for_experiment]
            y = [int(gene in valid_genes) for gene in experiment_genes]
            y = np.array(y)
            for row_idx in range(0, len(experiment_genes)):
                gene_control = control.T[row_idx].to_numpy()
                gene_perturbed = perturbed.T[row_idx].to_numpy()
                feature_vector = self.get_feature_vector(gene_control,
                                                         gene_perturbed)
                fy = y[row_idx]
                all_feature_vectors.append({"x": feature_vector, "y": fy})

        self.logger.info(f"writing data to file {out_file_name}")
        with open(out_file_name, 'wb') as f:
            pickle.dump(all_feature_vectors, f)

    def buildModel(self):
        pass

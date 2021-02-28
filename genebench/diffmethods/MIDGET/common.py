from utils import Utils
import numpy as np
import pickle
import random
from diffmethods.base.diffmethod import DiffMethod


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

    def __get_all_raw_geo(self):
        all_data = []
        geo_datas = self.storage.get_geo_data()
        for data in geo_datas:
            geo_id = data["name"]
            valid_c, control = Utils.filter_data(self.logger, data["control"])
            valid_p, perturbed = Utils.filter_data(self.logger, data["perturbed"])

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

    def get_feature_vector_size(self):
        return Utils.get_feature_vector_size()

    def get_feature_vector(self, control_gene, perturbed_gene):
        return Utils.get_feature_vector(control_gene, perturbed_gene)

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

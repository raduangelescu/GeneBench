import time
import sys
import os
import json
import numpy as np
import GEOparse
import pandas as pd
from genebench.utils import Utils
from genebench.storage.storage import Storage
from genebench.datatypes import GeoData

class LabelingConfig:
    def __init__(self, control, type, gene_names,
                 no_column_title, no_column_control, no_column_accession):
        self.control = control
        self.type = type
        self.gene_names = gene_names
        self.no_column_control = no_column_control
        self.no_column_title = no_column_title
        self.no_column_accession = no_column_accession


class InputConfig:
    def __init__(self,
                 pf_field,
                 name,
                 file):
        self.name = name
        self.file = file
        self.data = None
        self.pf_field = pf_field

    def load(self):
        with open(self.file, "r") as read_stream:
            self.data = json.load(read_stream)['data']


class GEOImporterConfig:
    def __init__(self,
                 data_path,
                 input_data,
                 labeling,
                 geo_cache_folder_name,
                 download_retry_count,
                 download_wait_seconds_before_retry):
        self.data_path = data_path
        self.input_data = []
        for input in input_data:
            self.input_data.append(InputConfig(**input))
        self.labeling = LabelingConfig(**labeling)
        self.cache_folder = geo_cache_folder_name
        self.download_retry_count = download_retry_count
        self.retry_wait = download_wait_seconds_before_retry


class GEOImporter:
    def __init__(self, config_filename):
        GEOparse.set_verbosity("ERROR")
        self.config_filename = config_filename
        data_section = Utils.get_config(config_filename, 'GEOImporter')
        self.config = GEOImporterConfig(**data_section)
        self.logger = Utils.get_logger('GEOImporter')
        self.storage = Storage(config_filename)
        self.labels = self.config.labeling
        self.inputs = self.config.input_data
        self.control_labels = self.labels.control
        self.type_labels = self.labels.type
        self.gene_names = self.labels.gene_names
        self.path = self.config.data_path
        self.experiment_collumns = {}

    def __get_genes(self, gse):
        gene_label = Utils.find_in_array(self.gene_names, gse.table.columns)
        genes_read = gse.table[gene_label].tolist()
        return Utils.deduplicate_genes(genes_read)

    def __split_control_perturbed(self, gse, column_entry, type_idx):
        control_series = []
        perturbed_series = []
        for series_name, description in gse.columns.iterrows():
            if Utils.is_control(description[type_idx], self.control_labels):
                column_entry['control'] = description[type_idx]
                control_series.append(series_name)
            else:
                column_entry['perturbed'] = description[type_idx]
                perturbed_series.append(series_name)

        return [control_series, perturbed_series]

    def __fix_bad_data(self, data: GeoData):
        max_replicates = 6
        # limit data to most important because we don't have
        # enough memory to run the R methods with full data
        if len(data.control_array) < len(data.genes):
            # some experiments need to be transposed
            control = np.array(data.control_array).T.tolist()
            perturbed = np.array(data.perturbed_array).T.tolist()
        else:
            control = data.control_array
            perturbed = data.perturbed_array
        # for control pick the first replicates
        control = [x[:max_replicates] for x in control]
        # for perturbed pick the last replicates
        perturbed = [x[-max_replicates:] for x in perturbed]
        # the above pick was done to favorize timeseries experiments
        control = Utils.log_if_necessary(np.array(control))
        perturbed = Utils.log_if_necessary(np.array(perturbed))

        control = Utils.quantile_normalize(pd.DataFrame(control))
        perturbed = Utils.quantile_normalize(pd.DataFrame(perturbed))

        data.control_array = control.to_numpy().tolist()
        data.perturbed_array = perturbed.to_numpy().tolist()
        return data

    def __do_data_item(self, gse, geo_id, source_name, pf):
        if not hasattr(gse, 'columns'):
            geo_data = self.__do_no_colums_item(gse,
                                                geo_id,
                                                source_name,
                                                pf)
        else:
            geo_data = self.__do_columns_item(gse,
                                              geo_id,
                                              source_name,
                                              pf)

        return geo_data

    def __do_no_colums_item(self, gse, geo_id, source_name, pf):
        control = self.labels.no_column_control
        phenotype_data = gse.phenotype_data
        columns = phenotype_data.columns
        info_experiment_idx = columns.get_loc(self.labels.no_column_title)
        gsm_ids_idx = columns.get_loc(self.labels.no_column_accession)
        gsm_type = list(phenotype_data.values[:, info_experiment_idx])
        gsm_ids = list(phenotype_data.values[:, gsm_ids_idx])
        control_gsms = []
        perturbation_gsms = []
        raw_control_data = []
        raw_perturbed_data = []
        for idx in range(0, len(gsm_type)):
            gsm_id = gsm_ids[idx]
            table = gse.gsms[gsm_id].table
            value_idx = table.columns.get_loc('VALUE')
            values = gse.gsms[gsm_id].table.values[:, value_idx].tolist()
            if Utils.find_in_array(gsm_type[idx], control) != 'unknown':
                control_gsms.append(gsm_id)
                raw_control_data.append(values)
            else:
                perturbation_gsms.append(gsm_id)
                raw_perturbed_data.append(values)

        if not control_gsms:
            self.logger('[no col]no control for {geo_id}')
            return None

        genes = gse.gsms[control_gsms[0]].table.values[:, 0]
        np_control_raw = np.array(raw_control_data)
        np_perturbed_raw = np.array(raw_perturbed_data)

        control = Utils.repair_nan_fast(np_control_raw)
        perturbed = Utils.repair_nan_fast(np_perturbed_raw)

        self.logger.info(f'finished {geo_id}')
        geo_data = GeoData({"name": geo_id,
                            "genes": genes.tolist(),
                            "source": source_name,
                            "perturbed_series_names": perturbation_gsms,
                            "control_series_names": control_gsms,
                            "extra_info": gse.metadata,
                            "perturbed_array": perturbed.tolist(),
                            "control_array": control.tolist(),
                            "pf": pf})
        return geo_data

    def __do_columns_item(self, gse, geo_id, source_name, pf):
        iter_labels = list(self.type_labels)
        type_labels = self.type_labels
        column_entry = {}
        column_entry['all'] = iter_labels

        while(iter_labels):
            type_label = Utils.find_in_array(gse.columns, iter_labels)
            iter_labels.pop(0)

            if type_label == 'unknown':
                error_msg = f'no label geoid {geo_id} labels:{type_labels}'
                self.logger.error(error_msg)
                continue

            type_idx = gse.columns.columns.get_loc(type_label)
            gene_label = self.__get_genes(gse)

            if gene_label == 'unknown':
                self.logger.error(
                    f'no gene label for geoid {geo_id} labels:{type_labels}')
                continue

            control_series, perturbed_series = self.__split_control_perturbed(
                gse, column_entry, type_idx)

            if not control_series:
                continue

            if not perturbed_series:
                continue

            np_control_raw = gse.table[control_series].to_numpy()
            np_perturbed_raw = gse.table[perturbed_series].to_numpy()

            control = Utils.repair_nan_fast(np_control_raw)
            perturbed = Utils.repair_nan_fast(np_perturbed_raw)

            self.logger.info(f'finished {geo_id}')
            geo_data = GeoData({"name": geo_id,
                                "genes": gene_label,
                                "source": source_name,
                                "perturbed_series_names": perturbed_series,
                                "control_series_names": control_series,
                                "extra_info": gse.metadata,
                                "perturbed_array": perturbed.tolist(),
                                "control_array": control.tolist(),
                                "pf": pf})

            self.experiment_collumns[geo_id] = column_entry

            return geo_data

        error_msg = (
            f"could not split {geo_id} in 2 classes"
            f"high cols: {gse.columns}"
        )
        self.logger.error(error_msg)
        return None

    def __download_retry(self, geo_id, cache_path):
        gse = None
        retry_count = 0
        while True:
            self.logger.info(f"Downloading {geo_id}")
            self.logger.info(f"--Retry count {retry_count}")
            try:
                gse = GEOparse.get_GEO(geo=geo_id, destdir=cache_path)
            except IOError as err:
                self.logger.warning(f"Error downloading geo data {err}")
                if retry_count > self.config.download_retry_count:
                    gse = None
                    self.logger.error(f"Could not download {geo_id}")
                    break
                self.logger.warning(
                    f"Waiting for {self.config.retry_wait} seconds")
                time.sleep(self.config.retry_wait)
                continue
            break
        return gse

    def __do_input(self, input):
        cache_folder = self.config.cache_folder
        cache_path = os.path.join(self.path, cache_folder, input.name)
        created_folder = Utils.create_folder_if_not_exist(cache_path)
        if created_folder:
            self.logger.info(f"created directory {cache_path}")
        log_data = {}
        for data_item in input.data:
            geo_id = data_item['geoid']
            pf = data_item[input.pf_field]
            info_msg = f'Getting GEO: {geo_id} in cache folder {cache_path}'
            self.logger.info(info_msg)
            gse = self.__download_retry(geo_id, cache_path)
            if gse is None:
                sys.exit(f"Failed to download data for {geo_id}")
            geo_data = self.__do_data_item(gse,
                                           geo_id,
                                           input.name,
                                           pf)
            if geo_data:
                geo_data = self.__fix_bad_data(geo_data)
                self.storage.insert_geo(geo_data)
        self.logger.info('Writing collumns to json file')
        log_str = json.dumps(log_data, sort_keys=True, indent=4)
        self.logger.info(f"{log_str}")
        self.logger.info(f'Finished importing GEO data for {input.name}')

    def importGEOData(self):
        labels = self.config.labeling
        inputs = self.config.input_data
        control_labels = self.labels.control
        type_labels = labels.type
        gene_names = labels.gene_names

        log_message = (
            f"reading geos"
            f"--control labels:{control_labels}"
            f"--type labels: {type_labels}"
            f"--gene names: {gene_names}"
        )
        self.logger.info(log_message)
        for input in inputs:
            self.logger.info(f"Loading data from file {input.file}")
            input.load()
            self.__do_input(input)

import csv
import os
import re
from genebench.storage.storage import Storage
from genebench.utils import Utils
from genebench.datatypes import GeneDiffValidation


class ValidationDataImporterConfig:
    def __init__(self, base_path, sources):
        self.base_path = base_path
        self.sources = sources


class ValidationDataImporter:

    def __init__(self, config):
        self.collect = {}
        self.storage = Storage(config)
        validation_section = Utils.get_config(config, 'ValidationDataImporter')
        self.config = ValidationDataImporterConfig(**validation_section)
        self.logger = Utils.get_logger('ValidationDataImporter')

    def importValidationData(self):
        self.logger.info("importing validation data")
        sources = self.config.sources

        self.logger.info("importing drug gene interaction")
        drug_gene_interaction = sources['DRUG GENE INTERACTION']
        drug_file = drug_gene_interaction['data_file']
        dgidb_data = self.drug_gene_interaction_db(self.config.base_path,
                                                   drug_file)

        self.logger.info("importing encode")
        encode_source = sources['ENCODE']
        encode_data = self.__from_attribute_matrix('encode',
                                                   encode_source['dict_file'],
                                                   encode_source['data_file'],
                                                   self.config.base_path)

        self.logger.info("importing chea")
        chea_source = sources['CHEA']
        chea_data = self.__from_attribute_matrix('chea',
                                                 chea_source['dict_file'],
                                                 chea_source['data_file'],
                                                 self.config.base_path)

        self.logger.info("importing pp interaction")
        ppi_int = sources['PP INTERACTION']
        ppi_int_data = self.ppi_interaction(self.config.base_path,
                                            ppi_int['data_file'])

        self.logger.info("importing ppi study")
        ppi_tfs = sources['PPI STUDY']
        ppi_tfs_data = self.ppi_tf(self.config.base_path,
                                   ppi_tfs['data_file'])

        store_to_db = {
            'tf_encode': encode_data,
            'tf_chea':  chea_data,
            'pp_interaction': ppi_int_data,
            'ppis_study': ppi_tfs_data
        }
        store_to_db['all'] = self.merge_all(store_to_db)
        store_to_db['drug_gene_interaction'] = dgidb_data
        self.logger.info('storing to validation data..')
        for key, data in store_to_db.items():
            data = GeneDiffValidation({'source': key, 'data': data})
            self.storage.insert_validation(data)
        self.logger.info('done..')

    def __check_add_collection(self, collection, p1, p2):
        def check_add_one(p1, p2, collection):
            if p1 in collection:
                if p2 not in collection[p1]:
                    collection[p1].add(p2)
                    return True
            return False

        if p1 not in collection and p2 not in collection:
            collection[p1] = set([p2])
            return True

        if check_add_one(p1, p2, collection):
            return True

        if check_add_one(p2, p1, collection):
            return True

        return False

    def __load_dict_attributes(self, dict_path):
        attributes = {}
        with open(dict_path) as tsvfile:
            reader = self.__get_tsv_reader(tsvfile)
            for row in reader:
                attributes[row['GeneID']] = row['GeneSym']
        return attributes

    def __get_tsv_reader(self, tsvfile):
        return csv.DictReader(
            filter(lambda row: row[0] != '#', tsvfile), dialect='excel-tab')

    def __from_attribute_matrix(self, name, dict_file_path,
                                data_file_path, base_path):

        dict_path = os.path.join(base_path, dict_file_path)
        data_path = os.path.join(base_path, data_file_path)
        collect = {}
        new_links = 0
        new_tfs = 0
        attributes = {}
        self.logger.info(f"Doing db {name}")
        attributes = self.__load_dict_attributes(dict_path)

        with open(data_path) as tsvfile:
            reader = self.__get_tsv_reader(tsvfile)
            for row in reader:
                popGene = row['GeneSym'].lower()
                for key, value in row.items():
                    if key not in attributes:
                        continue
                    value_float = float(value)
                    if Utils.isclose(value_float, 0.0):
                        continue
                    tf_name = attributes[key].lower()
                    if tf_name in collect:
                        if popGene not in collect[tf_name]:
                            collect[tf_name].add(popGene)
                            new_links = new_links + 1
                    else:
                        collect[tf_name] = set([popGene])
                        new_tfs = new_tfs + 1
        message = (f"TF db:{name }",
                   f"new tfs: {new_tfs}",
                   f"new links: {new_links}")
        self.logger.info(message)
        return collect

    def ppi_interaction(self, base_path, data_file) -> GeneDiffValidation:
        collect = {}
        new_entries = 0
        pina_filename = os.path.join(base_path, data_file)
        with open(pina_filename) as tsvfile:
            reader = self.__get_tsv_reader(tsvfile)
            for row in reader:
                try:
                    proteinALong = row['Alt. ID(s) interactor A']
                    proteinBLong = row['Alt. ID(s) interactor B']
                    proteinAShort = re.search('uniprotkb:(.*)\(gene name\)',
                                              proteinALong).group(1).lower()
                    proteinBShort = re.search('uniprotkb:(.*)\(gene name\)',
                                              proteinBLong).group(1).lower()

                    if self.__check_add_collection(collect, proteinAShort,
                                                   proteinBShort):
                        new_entries = new_entries + 1
                except Exception:
                    print(str(row))
        self.logger.info(f'PP interaction network {new_entries} new links')
        return collect

    def ppi_tf(self, base_path, data_file) -> GeneDiffValidation:
        collect = {}
        new_entries = 0
        pina_filename = os.path.join(base_path, data_file)
        with open(pina_filename) as tsvfile:
            reader = csv.reader(tsvfile, delimiter="\t")
            for row in reader:
                for index, gene_name in enumerate(row):
                    if index == 0:
                        continue
                    if len(gene_name) == 0:
                        continue

                    if self.__check_add_collection(collect,
                                                   row[0].lower(),
                                                   gene_name.lower()):
                        new_entries = new_entries + 1
        self.logger.info(f'PP TF network {new_entries} new links')
        return collect

    def drug_gene_interaction_db(self, base_path, data_file):
        collect = {}
        DRUG_NAME_IDX = 7
        GENE_NAME_IDX = 0
        DRUG_CLAIM_NAME_IDX = 6
        filename = os.path.join(base_path, data_file)
        with open(filename) as tsvfile:
            tsvreader = csv.reader(tsvfile, delimiter="\t")
            for idx, line in enumerate(tsvreader):
                if idx == 0:
                    continue
                drug = line[DRUG_NAME_IDX].lower()
                gene = line[GENE_NAME_IDX].lower()
                if drug.isspace() or drug == '':
                    drug = line[DRUG_CLAIM_NAME_IDX].lower()
                if drug not in collect:
                    collect[drug] = []
                collect[drug].append(gene)
        return collect

    def merge_all(self, all_tf_dicts):
        collect = {}
        for source, source_data in all_tf_dicts.items():
            for tf_name, genes in source_data.items():
                if tf_name not in collect:
                    collect[tf_name] = genes
                else:
                    collect[tf_name] = collect[tf_name].union(genes)
        return collect


def main():
    importer = ValidationDataImporter('config.json')
    importer.importValidationData()


if __name__ == "__main__":
    main()

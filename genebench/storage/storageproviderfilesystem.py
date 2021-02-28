import os
import json
import glob
from typing import List
from genebench.utils import Utils
from genebench.datatypes import GeneDiffValidation, GeoData, GeneMethodResult


class FileSystemConfig:
    def __init__(self, dict):
        self.base_path = dict['base_path']
        self.validation_folder = dict['validation_folder']
        self.geo_folder = dict['geo_folder']
        self.results_folder = dict['results_folder']


class StorageProviderFileSystem:

    def __init__(self, config: dict):
        self.config = FileSystemConfig(config)
        base_path = self.config.base_path
        Utils.create_folder_if_not_exist(base_path)
        self.geo_path = os.path.join(base_path,
                                     self.config.geo_folder)
        Utils.create_folder_if_not_exist(self.geo_path)
        self.validation_path = os.path.join(base_path,
                                            self.config.validation_folder)
        Utils.create_folder_if_not_exist(self.validation_path)
        self.results_path = os.path.join(base_path,
                                         self.config.results_folder)
        Utils.create_folder_if_not_exist(self.results_path)

    def insert_validation(self, data: GeneDiffValidation):
        Utils.create_folder_if_not_exist(self.validation_path)
        out_path = os.path.join(self.validation_path, f"{data.source}.json")
        with open(out_path, "w") as out:
            json.dump({"source": data.source,
                       "data": data.data},
                      out,
                      cls=Utils.SetEncoder)

    def delete_validation(self, source):
        out_path = os.path.join(self.validation_path, f"{source}.json")
        Utils.delete_files_from_folder(out_path)

    def get_validation_data(self, source, tf) -> GeneDiffValidation:
        out_path = os.path.join(self.validation_path, f"{source}.json")
        if not os.path.isfile(out_path):
            return None
        data = self.__load_data_from_file(out_path, GeneDiffValidation)
        return data

    def get_validation_sources(self):
        files = Utils.list_filenames_in_folder(self.validation_path)
        return [x[:-5] for x in files]

    def insert_geo(self, data: GeoData):
        Utils.create_folder_if_not_exist(self.geo_path)
        folder = os.path.join(self.geo_path, data.source)
        Utils.create_folder_if_not_exist(folder)
        out_path = os.path.join(folder, f"{data.name}.json")
        with open(out_path, "w") as out:
            json.dump(data.get_as_dict(), out)

    def delete_geo(self):
        Utils.delete_files_from_folder(self.geo_path)

    def __load_data_from_file(self, path, class_type):
        if os.path.isfile(path):
            with open(path, "r") as in_file:
                json_l = json.load(in_file)
                return class_type(json_l)
        return None

    def __load_data_from_files(self, files, class_type):
        data = []
        for file_path in files:
            abs_file_path = os.path.join(file_path)
            content = self.__load_data_from_file(abs_file_path, class_type)
            if content:
                data.append(content)
        return data

    def __load_data_from_folder(self, path, class_type):
        files = Utils.list_files_in_folder(path)
        return self.__load_data_from_files(files, class_type)

    def get_geo(self, filter) -> List[GeoData]:
        if 'name' in filter:
            file_path = os.path.join(self.geo_path,
                                     filter['source'],
                                     f"{filter['name']}.json")
            return [self.__load_data_from_file(file_path, GeoData)]

        if 'source' in filter:
            folder = os.path.join(self.geo_path,
                                  filter['source'])
        else:
            folder = self.geo_path
            folders = Utils.list_folders_in_folder(folder)
            ret_data = []
            for fld in folders:
                fld_path = os.path.join(folder, fld)
                dt = self.__load_data_from_folder(fld_path, GeoData)
                ret_data.extend(dt)
            return ret_data

        return self.__load_data_from_folder(folder, GeoData)

    def insert_method_results(self,
                              result: GeneMethodResult,
                              method_name: str,
                              experiment_name: str):
        Utils.create_folder_if_not_exist(self.results_path)
        method_folder = os.path.join(self.results_path, method_name)
        Utils.create_folder_if_not_exist(method_folder)
        file_output = os.path.join(method_folder, f"{experiment_name}.json")

        with open(file_output, "w") as out:
            output = result.to_dict()
            json.dump(output, out)

    def get_method_results(self,
                           filter):
        path = self.results_path
        has_method_name = False
        has_experiment_name = False

        if 'method_name' in filter:
            path = os.path.join(self.results_path,
                                filter['method_name'])
            has_method_name = True

        if 'experiment_name' in filter:
            has_experiment_name = True

        if has_method_name is True and has_experiment_name is True:
            # single file in path
            exp = os.path.join(path,
                               f"{filter['experiment_name']}.json")
            return [self.__load_data_from_file(exp, GeneMethodResult)]

        elif has_method_name is True and has_experiment_name is False:
            # find experiment name in all folders
            exp = os.path.join(self.results_path,
                               f"*{filter['experiment_name']}.json")
            files = glob.glob(exp, recursive=True)
            return self.__load_data_from_files(files, GeneMethodResult)

        elif has_method_name is False and has_experiment_name is False:
            # get all data
            return self.__load_data_from_folder(path, GeneMethodResult)

        return []

    def has_method_results(self, method_name: str,
                           experiment_name: str):
        path = os.path.join(self.results_path,
                            method_name,
                            f"{experiment_name}.json")
        return os.path.isfile(path)

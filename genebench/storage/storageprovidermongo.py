from pymongo import MongoClient
import gridfs
import copy
import json
from genebench.storage.storageprovider import StorageProvider
from genebench.datatypes import GeoData, GeneMethodResult, GeneDiffValidation
from typing import List
from genebench.utils import Utils


class MongoConfig:
    def __init__(self,
                 user,
                 password,
                 host,
                 port,
                 database_name,
                 validation_collection_name,
                 geo_data_collection_name,
                 results_collection_name,
                 anon):
        self.user = user
        self.password = password
        self.host = host
        self.port = port
        self.database_name = database_name
        self.validation_collection_name = validation_collection_name
        self.geo_data_collection_name = geo_data_collection_name
        self.results_collection_name = results_collection_name
        self.anon = anon

class StorageProviderMongo(StorageProvider):
    def __init__(self, config: dict):
        self.config = MongoConfig(**config)
        if self.config.anon:
            self.client = MongoClient(host=self.config.host,
                                      port=self.config.port)
        else:
            self.client = MongoClient(host=self.config.host,
                                      port=self.config.port,
                                      username=self.config.user,
                                      password=self.config.password)

        self.database = self.client[self.config.database_name]
        self.filesystem = gridfs.GridFS(self.database)

    def __delete(self, section_name):
        self.database[section_name].delete_many({})

    def __insert(self, section_name, data):
        self.database[section_name].insert_one(data)

    def __insert_big_data(self, section_name, meta_data, big_data):
        json_string = json.dumps(big_data, cls=Utils.SetEncoder)
        file_id = self.filesystem.put(json_string.encode('utf-8'))
        save_meta_data = copy.copy(meta_data)
        save_meta_data['file_id'] = file_id
        self.__insert(section_name, save_meta_data)

    def __get_big_data(self, section_name, filter):
        all_data = list(self.__get_data(section_name, filter))
        for data in all_data:
            data_link = data['file_id']
            raw_binary = self.filesystem.get(data_link)
            utf8_data = raw_binary.read().decode("utf-8")
            data['file'] = json.loads(utf8_data)
        return all_data

    def __get_data(self, section_name, filter):
        return self.database[section_name].find(filter)

    def __get_data_distinct(self, section_name, filter, field):
        return self.database[section_name].find(filter).distinct(field)

    def insert_validation(self, data: GeneDiffValidation):
        section = self.config.validation_collection_name
        self.__insert_big_data(section, {"source": data.source},
                               data.data)

    def delete_validation(self):
        section = self.config.validation_collection_name
        self.__delete(section)

    def get_validation_data(self, source, tf) -> GeneDiffValidation:
        section = self.config.validation_collection_name
        filter = {'tf': tf, 'source': source}
        data = self.__get_data(section, filter)
        ret_data = GeneDiffValidation()
        ret_data.data = data['data']
        ret_data.source = data['source']
        return ret_data

    def get_validation_sources(self):
        section = self.config.validation_collection_name
        return self.__get_data_distinct(self, section, {}, "source")

    def delete_geo(self):
        section = self.config.geo_data_collection_name
        self.__delete(section)

    def insert_geo(self, data: GeoData):
        section = self.config.geo_data_collection_name
        self.__insert_big_data(section,
                               data.get_meta_data(),
                               data.get_big_data())

    def get_geo(self, filter) -> List[GeoData]:
        section = self.config.geo_data_collection_name
        return self.__get_big_data(self, section, filter)

    def insert_method_results(self,
                              result: GeneMethodResult,
                              method_name: str,
                              experiment_name: str):
        section = self.config.results_collection_name
        meta_data = {
                    'method_name': method_name,
                    'experiment_name': experiment_name}
        self.__insert_big_data(section,
                               meta_data,
                               result.to_dict())

    def get_method_results(self,
                           filter):
        section = self.config.results_collection_name
        data = self.__get_big_data(self, section, filter)['file']
        return GeneMethodResult(data)

    def has_method_results(self, method_name: str,
                           experiment_name: str):
        section = self.config.results_collection_name
        filter = {'method_name': method_name,
                  'experiment_name': experiment_name}
        results = self.database[section].find(filter)
        count = results.count()
        if count > 0:
            return True
        return False

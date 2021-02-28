from genebench.utils import Utils
from genebench.storage.storageprovidermongo import StorageProviderMongo
from genebench.storage.storageproviderfilesystem import StorageProviderFileSystem
from genebench.datatypes import GeneDiffValidation, GeoData, GeneMethodResult


class StorageConfig:
    def __init__(self,  dict):
        self.load_order = dict['load_order']
        self.providers = dict['providers']


class Storage:

    @staticmethod
    def create_provider(name, config):
        if name == 'mongo':
            return StorageProviderMongo(config)
        if name == 'filesystem':
            return StorageProviderFileSystem(config)
        return None

    def __init__(self, config_filename):
        self.logger = Utils.get_logger('Storage')
        config_section = Utils.get_config(config_filename, "Storage")
        self.config = StorageConfig(config_section)
        self.providers = {}
        providers = self.config.providers
        for provider_name, provider_json_config in providers.items():
            self.providers[provider_name] = Storage.create_provider(
                provider_name, provider_json_config)
        self.logger.info(f"Started storage with config {config_section}")

    def insert_validation(self, data: GeneDiffValidation):
        for provider_name in self.config.load_order:
            self.providers[provider_name].insert_validation(data)

    def delete_validation(self):
        for provider_name in self.config.load_order:
            self.providers[provider_name].delete_validation()

    def get_validation_data(self, source, tf) -> GeneDiffValidation:
        for provider_name in self.config.load_order:
            provider = self.providers[provider_name]
            return provider.get_validation_data(source, tf)

    def get_validation_sources(self):
        for provider_name in self.config.load_order:
            return self.providers[provider_name].get_validation_sources()

    def delete_geo(self):
        for provider_name in self.config.load_order:
            self.providers[provider_name].delete_geo()

    def insert_geo(self, data: GeoData):
        for provider_name in self.config.load_order:
            self.providers[provider_name].insert_geo(data)

    def get_geo(self, filter) -> GeoData:
        for provider_name in self.config.load_order:
            return self.providers[provider_name].get_geo(filter)

    def has_method_results(self,
                           method_name: str,
                           experiment_name: str):
        for provider_name in self.config.load_order:
            provider = self.providers[provider_name]
            has = provider.has_method_results(method_name,
                                              experiment_name)
            if has:
                return True
        return False

    def insert_method_results(self,
                              result: GeneMethodResult,
                              method_name: str,
                              experiment_name: str):

        for provider_name in self.config.load_order:
            provider = self.providers[provider_name]
            provider.insert_method_results(result,
                                           method_name,
                                           experiment_name)

    def get_method_results(self, filter):
        for provider_name in self.config.load_order:
            return self.providers[provider_name].get_method_results(filter)

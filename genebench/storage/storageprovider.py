from genebench.datatypes import GeneMethodResult, GeoData, GeneDiffValidation


class StorageProvider:
    def insert_validation(self, data: GeneDiffValidation):
        pass

    def delete_validation(self):
        pass

    def get_validation_data(self, source, tf) -> GeneDiffValidation:
        pass

    def get_validation_sources(self) -> list:
        pass

    def delete_geo(self):
        pass

    def insert_geo(self, data: GeoData):
        pass

    def get_geo(self, filter) -> GeoData:
        pass

    def insert_method_results(self,
                              result: GeneMethodResult,
                              method_name: str,
                              experiment_name: str):
        pass

    def get_method_results(self,
                           method_name: str,
                           experiment_name: str):
        pass

    def has_method_results(method_name: str,
                           experiment_name: str):
        pass

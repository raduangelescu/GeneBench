
from genebench.datatypes import GeneData, GeneDiffValidation


class SilicoDataGeneratorConfig:
    def __init__(self,
                 name,
                 source,
                 module_name,
                 class_name,
                 params):
        self.source = source
        self.module_name = module_name
        self.class_name = class_name
        self.params = params
        self.name = name


class SilicoDataGenerator:
    def __init__(self, config):
        self.config = SilicoDataGeneratorConfig(**config)

    def generate_experiments(self,
                             validation_data,
                             num_experiments,
                             num_genes) -> GeneData:
        return GeneData()

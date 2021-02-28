from genebench.datatypes import GeneData
from genebench.datatypes import GeoData
from genebench.datatypes import GeneDiffValidation
from genebench.utils import Utils
import importlib
import random


class SilicoGeneratorEntry:
    def __init__(self,
                 name,
                 source,
                 module_name,
                 class_name,
                 params):
        self.module_name = module_name
        self.class_name = class_name
        self.params = params
        self.name = name
        self.source = source


class SilicoGeneratorsManagerConfig:

    def __init__(self, dict):
        self.generators = []
        for name, config_item in dict.items():
            self.generators.append(SilicoGeneratorEntry(name=name,
                                                        **config_item))


class SilicoGeneratorsManager:

    def __init__(self, config_filename):
        self.logger = Utils.get_logger("SilicoGeneratorsManager")
        self.config = Utils.get_config(config_filename, "SilicoGenerators")
        self.config = SilicoGeneratorsManagerConfig(self.config)
        self.generator_instances = {}
        self.module_names = []

    def register_diff_method(self, config):
        self.generators.append(SilicoGeneratorEntry(**config))

    def setup(self):
        if self.module_names:
            for module_name in self.module_names:
                importlib.import_module(module_name)

        for gen_entry in self.config.generators:
            module = importlib.import_module(gen_entry.module_name)
            class_desc = getattr(module, gen_entry.class_name)
            class_instance = class_desc(gen_entry.__dict__)
            self.generator_instances[gen_entry.name] = class_instance

    def generate_validation_data(self,
                                 num_genes,
                                 num_pfs) -> GeneDiffValidation:
        gene_names = Utils.get_random_gene_names(num_genes)
        pf_names = Utils.get_random_tf_names(num_pfs)
        data = {}
        for pf in pf_names:
            num = random.randint(2, len(gene_names))
            choice_list = random.choices(gene_names, k=num)
            choice_set = set(choice_list)
            choice_list = list(choice_set)
            data[pf] = choice_list

        validation_data = GeneDiffValidation()
        validation_data.source = "silico"
        validation_data.data = data
        return validation_data

    def generate_experiments(self,
                             generator_name,
                             validation_data: GeneDiffValidation,
                             num_experiments,
                             num_genes):
        if generator_name not in self.generator_instances:
            self.logger.info(f"Generator {generator_name} not found!")
            return GeneData()

        instance = self.generator_instances[generator_name]
        return instance.generate_experiments(validation_data,
                                             num_experiments,
                                             num_genes)

    def generate_all_experiments(self,
                                 validation_data: GeneDiffValidation,
                                 num_experiments,
                                 num_genes):
        all_data_geo = []
        for generator_name in self.generator_instances.keys():
            self.logger.info(f"Generating data with {generator_name}")
            gene_data = self.generate_experiments(generator_name,
                                                  validation_data,
                                                  num_experiments,
                                                  num_genes)
            all_data_geo.append(gene_data)
        return all_data_geo

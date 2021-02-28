from genebench.datatypes import GeneDiffInput
from genebench.datatypes import GeneMethodResult
from genebench.utils import Utils
import importlib


class DiffMethodEntry:
    def __init__(self,
                 module_name,
                 class_name,
                 config,
                 feature_file=''):
        self.module_name = module_name
        self.class_name = class_name
        self.config = config
        self.feature_file = feature_file


class DiffMethodsManagerConfig:

    def __init__(self, dict):
        self.methods = {}
        for key, config_item in dict.items():
            self.methods[key] = DiffMethodEntry(**config_item)


class DiffMethodsManager:

    def __init__(self, config_filename):
        self.logger = Utils.get_logger("DiffMethodManager")
        self.config = Utils.get_config(config_filename, "Methods")
        self.config = DiffMethodsManagerConfig(self.config)
        self.method_instances = {}
        self.module_names = []

    def register_diff_method(self, config):
        self.methods.append(DiffMethodEntry(**config))

    def setup(self):
        if self.module_names:
            for module_name in self.module_names:
                importlib.import_module(module_name)

        for name, method_entry in self.config.methods.items():
            module = importlib.import_module(method_entry.module_name)
            class_desc = getattr(module, method_entry.class_name)
            class_instance = class_desc()
            self.method_instances[name] = class_instance
            class_instance.setup(method_entry.config)

    def train(self, method_name: str):
        if method_name not in self.method_instances:
            self.logger(f"Method {method_name} not found!")
            return False

        instance = self.method_instances[method_name]
        instance.train()
        return True

    def run(self, gene_input: GeneDiffInput, method_name: str):
        if method_name not in self.method_instances:
            self.logger(f"Method {method_name} not found!")
            return GeneMethodResult()

        instance = self.method_instances[method_name]
        return instance.run(gene_input)

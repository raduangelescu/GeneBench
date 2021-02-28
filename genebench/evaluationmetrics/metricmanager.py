from genebench.datatypes import GeneData
from genebench.datatypes import GeoData
from genebench.datatypes import GeneDiffValidation
from genebench.utils import Utils
import importlib


class MetricEntry:
    def __init__(self,
                 name,
                 module_name,
                 class_name,
                 params):
        self.module_name = module_name
        self.class_name = class_name
        self.params = params
        self.name = name


class MetricManagerConfig:

    def __init__(self, dict):
        self.metrics = []
        self.output_folder = dict['output_folder']
        for _, config_item in dict['metrics'].items():
            self.metrics.append(MetricEntry(**config_item))


class MetricManager:

    def __init__(self, config_filename):
        self.logger = Utils.get_logger("MetricManager")
        self.config = Utils.get_config(config_filename, "AccuracyMetrics")
        self.config = MetricManagerConfig(self.config)
        self.metric_instances = {}
        self.module_names = []

    def register_metric_method(self, config):
        self.config.metrics.append(MetricEntry(**config))

    def setup(self):
        if self.module_names:
            for module_name in self.module_names:
                importlib.import_module(module_name)

        for metric_entry in self.config.metrics:
            module = importlib.import_module(metric_entry.module_name)
            class_desc = getattr(module, metric_entry.class_name)
            class_instance = class_desc(metric_entry.__dict__,
                                        self.config.output_folder)
            self.metric_instances[metric_entry.name] = class_instance

    def add(self,
            pf,
            method_name,
            valid,
            res):
        for _, instance in self.metric_instances.items():
            instance.add(pf,
                         method_name,
                         valid,
                         res)

    def evaluate(self, group_name):
        for _, instance in self.metric_instances.items():
            instance.evaluate(group_name)

from genebench.datatypes import GeneDiffInput, GeneData
from genebench.utils import Utils
from genebench.diffmethods.diffmethodsmanager import DiffMethodsManager
from genebench.storage.storage import Storage
from genebench.evaluationmetrics.metricmanager import MetricManager


class BenchmarkMethodGroupConfig:
    def __init__(self, name, methods):
        self.name = name
        self.methods = methods


class BenchmarkRunConfig:
    def __init__(self,
                 config):
        self.name = config["name"]
        self.data_sources = config["data_sources"]
        self.is_valid = True
        self.method_group_ids = config["method_group_ids"]
        self.validation_sets = config["validation_sets"]


class BenchamarkDiffMethodsConfig:
    def __init__(self, logger, method_groups, runs):
        self.method_groups = {}
        for name, method_group in method_groups.items():
            self.method_groups[name] = BenchmarkMethodGroupConfig(
                **method_group)

        self.runs = []
        for run in runs:
            config = BenchmarkRunConfig(run)
            if config.is_valid:
                self.runs.append(config)
            else:
                logger.error(f"run {config.name} has a bad type {config.type}")


class BenchmarkDiffMethods:

    def __init__(self, config_filename):
        self.logger = Utils.get_logger("Benchmark")
        self.config = Utils.get_config(config_filename, "Benchmark")
        self.config = BenchamarkDiffMethodsConfig(logger=self.logger,
                                                  **self.config)
        self.method_manager = DiffMethodsManager(config_filename)
        self.method_manager.setup()
        self.storage = Storage(config_filename)
        self.metric_manager = MetricManager(config_filename)
        self.metric_manager.setup()

    def run_method(method_name: str, input: GeneDiffInput):
        pass

    def get_execution_map(self):
        # collect all methods so we don't run them twice on same inputs
        execution_map = {}

        for run in self.config.runs:
            for _, data in self.config.method_groups.items():
                for method_name in data.methods:
                    if method_name not in execution_map:
                        entry = set()
                        execution_map[method_name] = entry
                    else:
                        entry = execution_map[method_name]
                    new_sources = set(run.data_sources)
                    execution_map[method_name] = entry.union(new_sources)
        return execution_map

    def generate_method_results(self):
        execution_map = self.get_execution_map()
        data_source_cache = {}
        for method_name, execution_data in execution_map.items():
            for source in execution_data:
                self.logger.info(f"get data from: {source}")
                if source in data_source_cache:
                    experiments = data_source_cache[source]
                geos = self.storage.get_geo({'source': source})
                experiments = []
                for exp in geos:
                    gen_data = GeneData()
                    gen_data.gene_input = GeneDiffInput.from_geo_data(exp)
                    gen_data.name = exp.name
                    experiments.append(gen_data)
                data_source_cache[source] = experiments
                num_experiments = len(experiments)

                for id, gene_data in enumerate(experiments):
                    if self.storage.has_method_results(method_name,
                                                       gene_data.name):
                        exp_name = gene_data.name
                        self.logger.info(f"already computed [{exp_name}]")
                        continue
                    self.logger.info(f"running {id+1}/{num_experiments}")
                    self.logger.info(f"{method_name}[{gene_data.name}]")

                    res = self.method_manager.run(gene_data.gene_input,
                                                  method_name)
                    self.storage.insert_method_results(res,
                                                       method_name,
                                                       gene_data.name)

    def generate_comparison_single(self, method_name, geodata, run, cache):
        filter = {'method_name': method_name,
                  'experiment_name': geodata.name}
        res = self.storage.get_method_results(filter)[0]
        for validation_set in run.validation_sets:
            self.logger.info(f"Getting validation data for {validation_set}")
            if validation_set not in cache:
                valid = self.storage.get_validation_data(validation_set,
                                                         geodata.pf)
                cache[validation_set] = valid

            valid = cache[validation_set]

            self.logger.info("Adding data to metrics")
            self.metric_manager.add(geodata.pf,
                                    method_name,
                                    valid,
                                    res)

    def generate_comparisons(self):
        validation_cache = {}
        self.logger.info("collecting metrics")
        for run in self.config.runs:
            all_geos = []
            for data_source in run.data_sources:
                all_geos.extend(self.storage.get_geo({
                    "source": data_source}))
            for key, data in self.config.method_groups.items():
                for method_name in data.methods:
                    for geodata in all_geos:
                        self.generate_comparison_single(method_name,
                                                        geodata,
                                                        run,
                                                        validation_cache)
                self.logger.info(f"evaluating metrics for {key}")
                self.metric_manager.evaluate(key)
        self.logger.info("finished generating metrics")

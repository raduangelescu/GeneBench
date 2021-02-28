from genebench.silico.silicogeneratorsmanager import SilicoGeneratorsManager
from genebench.storage.storage import Storage
from genebench.utils import Utils


class GenerateSilicoConfig():
    def __init__(self,
                 num_genes,
                 num_experiments,
                 num_pfs):
        self.num_genes = num_genes
        self.num_experiments = num_experiments
        self.num_pfs = num_pfs


class GenerateSilicoData():
    def __init__(self, config_filename):
        self.logger = Utils.get_logger('generate_and_store_silico')
        config_raw = Utils.get_config(config_filename, "SilicoData")
        self.config = GenerateSilicoConfig(**config_raw)
        self.store = Storage(config_filename)
        self.generator_manager = SilicoGeneratorsManager(config_filename)
        self.generator_manager.setup()

    def run(self):
        self.logger.info("generating silico data")
        num_genes = self.config.num_genes
        num_experiments = self.config.num_experiments
        num_pfs = self.config.num_pfs
        self.logger.info(f"    --number of genes:{num_genes}")
        self.logger.info(f"    --number of experiments:{num_experiments}")
        self.logger.info(f"    --number of pfs:{num_pfs}")

        self.logger.info("generating validation data")
        validation = self.generator_manager.generate_validation_data(
            num_genes=num_genes,
            num_pfs=num_pfs)
        self.logger.info("storing validation data")
        self.store.insert_validation(validation)
        self.logger.info("generating all fake experiments")
        providers = self.generator_manager.generate_all_experiments(
                                                             validation,
                                                             num_experiments,
                                                             num_genes)
        self.logger.info("storing experiments data")
        for id, provider in enumerate(providers):
            self.logger.info(f"silico provider: {id+1}/{len(providers)}")
            for id, experiment in enumerate(provider):
                self.logger.info(f"silico experiment {id+1}/{len(provider)}")
                self.store.insert_geo(experiment)
        self.logger.info("done")


def main():
    generate = GenerateSilicoData("config.json")
    generate.run()


if __name__ == "__main__":
    main()

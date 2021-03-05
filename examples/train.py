from genebench.diffmethods.MIDGET.neural import MIDGETNeural
from genebench.diffmethods.MIDGET.xgb import MIDGETXgBoost
from genebench.storage.storage import Storage
from genebench.utils import Utils
from genebench.diffmethods.diffmethodsmanager import DiffMethodsManager
from genebench.benchmark import BenchmarkDiffMethods

import numpy as np
import pickle


def generate_all_feature_vectors(out_file_name):
    storage = Storage("config.json")
    logger = Utils.get_logger("generate_feature_vectors")
    logger.info("getting all validation sources")
    valid_sources = list(storage.get_validation_sources())
    logger.info(f"we have {len(valid_sources)} validation sources")
    logger.info("filtering validation sources")

    valid_sources = list(filter(lambda x: x != 'all',
                                valid_sources))
    logger.info(f"we now have {len(valid_sources)} validation sources")
    validation_cache = {}
    logger.info("getting all geos")
    geos = storage.get_geo({})
    logger.info(f"we have total {len(geos)}")
    logger.info("filtering geos")
    #geos = list(filter(lambda x: 'silico' not in x.source, geos))
    logger.info(f"after filtering geos we have {len(geos)}")
    logger.info("populating validation cache")
    for idx, valid_source in enumerate(valid_sources):
        logger.info(f"adding to cache {idx}/{len(valid_sources)}")
        # collect all pfs
        pfs = set()
        for geo in geos:
            pf = geo.pf.lower()
            pfs.add(pf)
        
        for pf in pfs:
            valid = storage.get_validation_data(valid_source, pf)
            valid_data = valid.data
            if pf in valid_data:
                if pf not in validation_cache:
                    validation_cache[pf] = set()
                genes = valid_data[pf]
                genes = [x.lower() for x in genes]
                validation_cache[pf].update(genes)

    logger.info("generating training data from geos")
    all_data = []
    midget_instance = MIDGETNeural()
    running_sum = 0
    for idx, geo in enumerate(geos):
        pf = geo.pf.lower()
        if pf not in validation_cache:
            logger.error(f"could not find {pf} in validation set, skipping")
            continue
        logger.info(f"parsing geo {idx}/{len(geos)}")
        control = np.array(geo.control_array)
        perturbed = np.array(geo.perturbed_array)
        valid_genes = validation_cache[pf]
        geo_genes = []
        for gene in geo.genes:
            if isinstance(gene, str):
                geo_genes.append(gene.lower())
            else:
                geo_genes.append("__")
        y = [int(gene in valid_genes) for gene in geo_genes]
        y = np.array(y)
        running_sum += np.sum(y)
        num_genes = len(geo.genes)
        for row_idx in range(0, num_genes):
            feature_vector = midget_instance.get_feature_vector(control[row_idx],
                                                      perturbed[row_idx])
            fy = y[row_idx]
            all_data.append({"x": feature_vector, "y": fy})
        logger.info(f"current positives {running_sum}")
    logger.info(f"writing file with {running_sum} positives")
    with open(out_file_name, 'wb') as f:
        pickle.dump(all_data, f)
    logger.info("done")


def train():
    methods_manager = DiffMethodsManager("training_config.json")
    methods_manager.setup()
    methods_manager.train_all()

def train_and_bench():
    methods_manager = DiffMethodsManager("training_config.json")
    methods_manager.setup()
    methods_manager.train_all()

    benchmark = BenchmarkDiffMethods("config.json")
    # run all methods and store results
    benchmark.generate_method_results()
    # generate the actual metrics and plots for comparison
    benchmark.generate_comparisons()

def main():
    generate_all_feature_vectors("data/feature_vectors.pkl")
    train_and_bench()


if __name__ == "__main__":
    main()

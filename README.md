![Image logo](https://raw.githubusercontent.com/raduangelescu/GeneBench/main/logo.svg)

## Overview

Gene Bench is a benchmark-ing framework used in analyzing methods that detect differentially expressed genes from biological samples.
Besides being a benchmarking framework it also contains some commonly used Differential genes detection algorithms:
  
  - [LIMMA](https://bioconductor.org/packages/release/bioc/html/limma.html) (R implementation)
  - [Characteristic Direction Py](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-79) (a python characteristic direction that is good for big data, implementation using Nipals PCA)
  - [Characteristic Direction](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-79) (R implementation **crashes on large data** with out of memory because of normal PCA )
  - [T-Test](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html)
  - [SAM](https://cran.r-project.org/web/packages/samr/index.html)

Our own methods based on machine learning:
  - **MIDGET-xgb** (2 variations)
  - **MIDGET-neural** (6 variations)
and a baseline random algorithm.

The framework also provides a base for calling methods from R (because most medical laboratories use R for statistical analysis). 

We provide already implemented evaluation metrics:
 - **Kolmogorov** based visual evaluation
 - **F1** score
 - **ROC/AUC** curve analysis

Storage providers:
 - [mongodb](https://www.mongodb.com/) usefull for sharing amongst resourcers and fully configurable 
 - filesystem (via json files and folder structure)

And a silico data generator that is based on an easy to understand linear method.

You can also provide your own implementations for
 - methods
 - evaluation metrics
 - silico data generators
 - storage providers

You may configure any combination of evaluation metrics, storage providers, silico data generators and methods for your analysis and benchmark.
**Note 0** Please cites this framework in your papers if you plan to use it.
**Note 1** The package was tested on Linux Ubuntu 20.0. It should work in Windows too but you need to install dependencies with the normal installers. 
**Note 2** If you want good performance in evaluating MIDGET Neural methods we advise you use a GPU with CUDA suport and [do the necessary steps to activate it with tensor flow](https://www.tensorflow.org/install/gpu)

## Installation
Install using pip 
```bashs
pip install genebench
```

or install latest from github

```bash
git clone https://github.com/raduangelescu/GeneBench.git
python setup.py install

```

## Usage
The most important thing to note is that you need to have a correct config json file that describes what you are trying to do. The next sections will explain the code and config needed to access the framework features.
For easy usage, you may find the below examples in [this folder](https://github.com/raduangelescu/GeneBench/tree/main/examples)

### Quick start
To get the default working you only need to install the package.
- [Install MongoDb](https://docs.mongodb.com/manual/installation/)
```bash
wget -qO - https://www.mongodb.org/static/pgp/server-4.4.asc | sudo apt-key add -
sudo apt-get install gnupg
wget -qO - https://www.mongodb.org/static/pgp/server-4.4.asc | sudo apt-key add -
echo "deb [ arch=amd64,arm64 ] https://repo.mongodb.org/apt/ubuntu focal/mongodb-org/4.4 multiverse" | sudo tee /etc/apt/sources.list.d/mongodb-org-4.4.list
sudo apt-get update
sudo apt-get install -y mongodb-org
```
- Install R for the R-based methods (SAM, Limma)
```bash
sudo apt install r-base
```

-  Open R and Install packages
	- [Install Bioconductor](https://www.bioconductor.org/install/)

	```r
 		install.packages("BiocManager")
	```

	- [Install LIMMA in R](https://bioconductor.org/packages/release/bioc/html/limma.html)
	```r
   		BiocManager::install("limma")
	```
	- [Install SAM in R and impute](http://www.sthda.com/english/articles/2-r/1-install-samr-package/)
	```r
		install.packages('samr')
		source("http://bioconductor.org/biocLite.R")
		biocLite("impute")
	```
- Copy [this config.json](https://raw.githubusercontent.com/raduangelescu/GeneBench/main/examples/config.json) into your project folder
- make a folder named "data"  in your project folder
- Run startup function
```python
from genebench.startup import setup_default_data
setup_default_data("config.json") 
```
### Storage config
Any benchmark run via this framework needs storage. You may use any of the storage providers implemented in the framework, or both. 
The config file needs the **"Storage"** section in order to know how you will set and read the data. The below config is an example of using both providers.
```json
"Storage":{
	"providers":{
		"mongo": {
			"host":"localhost",
			"port":27017,
			"anon":true,
			"user":"",
			"password":"",
			"database_name":"gene_expression_bench_test",
			"validation_collection_name":"validation",
			"geo_data_collection_name":"geo",
			"results_collection_name":"results"
		},
		"filesystem": {
			"base_path":"data/db_test",
			"validation_folder":"validation",
			"geo_folder": "geo",
			"results_folder": "results"
		}
	},
	"load_order":[
		"filesystem",
		"mongo"
	]
}
```
If you only need to use one provider just remove the other from the above config form "providers" and "load_order". In most cases, when you don't have any mongodb installed you will probably only want to use the filesystem provider which should work locally.

### Generating silico data
The silico data generator tries to mimick a biological experiment: it will generate fake genes and fake perturbation factors and then create links between them. Based on the above biological network it will use silico providers to generate the fake data.

To generate the silico data you need to have the **"SilicoData"** and **"SilicoGenerators"** config sections in your config.json file
```json
"SilicoData":{
	"num_genes":10000,
	"num_experiments":100,
	"num_pfs":100
},
"SilicoGenerators":{
	"silico.linear":{
		"source": "silico_linear",
		"module_name": "genebench.silico.generators.linear",
		"class_name":"LinearDataGenerator",
		"params": {
					"diff_factor":3.0,
					"noise_factor":0.5,
					"num_replicates":3
		}
	}
}
```

After providing the correct config file you only need to run the below code:
```python
# examples/silicotest.py 
from genebench.generatesilico import GenerateSilicoData

generate = GenerateSilicoData("config.json")
generate.run()
```
The generator run will populate the data in the storage provider. The data itself will be in a generic format, the only way to know it was generated by a silico provider is the source name from the config. **Be sure to add an unique source name for all silico generators so you can identify experiments later**


### Importing in vivo GEO data

For this feature you need the **GEOImporter** section in the config file:
Modifying the "labeling" section is not advised, and only required if you want to add new GEO data and filter it automatically.
To add new GEO input data you may create files similarly to **GEOdatasets_tf_perturbations.json** and add them into the input_data field, more details will be provided in the full documentation.
```json
"GEOImporter":{
	"data_path":"data",
	"download_retry_count": 10,
	"download_wait_seconds_before_retry": 5,
	"geo_cache_folder_name":"cache",
	"input_data":[
		{
			"pf_field": "tf",
			"name":"tf_perturbations",
			"file":"data/default_experiments/GEOdatasets_tf_perturbations.json"
		},
		{
			"pf_field": "drug",
			"name":"drug_perturbations",
			"file":"data/default_experiments/GEOdatasets_drugs_perturbations.json"
		}
	],
	"labeling":
	{
		"control":["control","Control", "wild type", "wildtype","Wild Type","Wildtype","wild-type","baseline", "control virus","uninduced", "untreated", "0 pM", "water", "Time of treatment 0 weeks", "none","pretreatment", "vehicle", "control virus"],
		"type":["genotype/variation","disease state", "protocol","growth protocol","dose", "agent", "description", "time"],
		"gene_names":"IDENTIFIER",
		"no_column_control":["MCF7/BUS, 0 pM","baseline", "minus", "Resistant", "untreated", "Untreated", "Control", "control", "nonresponder", "resistant", "vehicle", "saline","0 pM E2","DMS", "AO ", " AG", "_Veh_","N1+2","N3+4","N5+6"],
		"no_column_title": "title",
		"no_column_accession": "geo_accession"
	}
}
```
It is usually the case that you will not modify this config in a normal run so you can paste it from above, removing comments.
To generate and store the experiment data you just run the below code, using whatever file-name your config is to replace 'config.json'.
```python
# examples/importgeodata.py
from genebench.geoimporter import GEOImporter
importer = GEOImporter('config.json')
importer.importGEOData()
```

### Importing validation data
To be able to test if the methods you create/benchmark are better on the set of data, we need to have some validation data. To do this we provided multiple validation sources and a way to import them. 
The config section needed for this feature is **"ValidationDataImporter"**

```json
"ValidationDataImporter":{
	"base_path":"data/validation_sources",
	"sources":{
		"ENCODE":{
			"data_file":"encode_gene_attribute_matrix.txt",
			"dict_file":"encode_attribute_list_entries.txt"
		},
		"CHEA":{
			"data_file":"chea_gene_attribute_matrix.txt",
			"dict_file":"chea_attribute_list_entries.txt"
		},
		"PP INTERACTION":{
			"data_file":"PPI_interaction.tsv"
		},
		"PPI STUDY":{
			"data_file":"PPI_transcription_factors.txt"
		},
		"DRUG GENE INTERACTION":{
			"data_file":"drug_gene_interaction_db.tsv"
		}

	}
}
```

And the coded to import the data:

```python
# examples/importvalidationdata.py
from genebench.validationdataimporter import ValidationDataImporter
importer = ValidationDataImporter('config.json')
importer.importValidationData()
```

### Running benchmarks
After you have all your validation and experiment data in your storage providers (general format) you will now configure the benchmarks you want to run via the **"AccuracyMetrics"**, **"Methods"** and **"Benchmark"** config sections:
- In **AccuracyMetrics** we configure the metrics we want to generate for all benchmarks. Below we provide an example for using the default provided ones
```json
"AccuracyMetrics":{
	"output_folder":"data/metrics" ,
	"metrics":{
		"Kolmogorov": {
			"name": "Kolmogorov",
			"module_name": "genebench.evaluationmetrics.kolmogorov",
			"class_name":"Kolmogorov",
			"params":{
			}
		},
		"ROC": {
			"name": "ROC",
			"module_name": "genebench.evaluationmetrics.roc",
			"class_name":"ROC",
			"params":{}
		},
		"F1": {
			"name": "F1",
			"module_name": "genebench.evaluationmetrics.f1",
			"class_name":"F1",
			"params":{}
		}
	}
}
```
-  In **Methods** we add all the methods we want to use for benchmarking: Each entry has: a key by which the program will address the code, module_name: the module in which the method relies, class_name: the actual method class name and a custom config field so you can use parameters in your own method implementation. Below we provide a config that accounts for all our currently supported methods:
```json
"Methods":{
	"TTest":{
		"module_name":"genebench.diffmethods.ttest",
		"class_name":"TTest",
		"config":{}
	},
	"LIMMA":{
		"module_name":"genebench.diffmethods.limma",
		"class_name":"LIMMA",
		"config":{}
	},
	"SAM":{
		"module_name":"genebench.diffmethods.sam",
		"class_name":"SAM",
		"config":{}
	},
	"Random":{
		"module_name":"genebench.diffmethods.random",
		"class_name":"Random",
		"config": {}
	},
	"Characteristic direction py":{
		"module_name":"genebench.diffmethods.chdirpy",
		"class_name":"ChDirPy",
		"config":{}
	},
	"MIDGET Neural[n1]":{
		"module_name":"genebench.diffmethods.MIDGET.neural",
		"class_name":"MIDGETNeural",
		"config": {
			"data_split": 0.25,
			"method_name":"",
			"batch_size": 40000,
			"test_size": 40000,
			"number_of_epochs": 1000,
			"save_point_num_batches": 10,
			"model_name": "n1",
			"output_folder":"data/MIDGETNeural" ,
			"feature_file":"features.pkl"
		}
	},
	"MIDGET Neural[n2]":{
		"module_name":"genebench.diffmethods.MIDGET.neural",
		"class_name":"MIDGETNeural",
		"config": {
			"data_split": 0.25,
			"method_name":"",
			"batch_size": 40000,
			"test_size": 40000,
			"number_of_epochs": 1000,
			"save_point_num_batches": 10,
			"model_name": "n2",
			"output_folder":"data/MIDGETNeural",
			"feature_file":"features.pkl"
		}
	},
	"MIDGET Neural[n3]":{
		"module_name":"genebench.diffmethods.MIDGET.neural",
		"class_name":"MIDGETNeural",
		"config": {
			"data_split": 0.25,
			"method_name":"",
			"batch_size": 40000,
			"test_size": 40000,
			"number_of_epochs": 1000,
			"save_point_num_batches": 10,
			"model_name": "n3",
			"output_folder":"data/MIDGETNeural",
			"feature_file":"features.pkl"
		}
	},
	"MIDGET Neural[n4]":{
		"module_name":"genebench.diffmethods.MIDGET.neural",
		"class_name":"MIDGETNeural",
		"config": {
			"data_split": 0.25,
			"method_name":"",
			"batch_size": 40000,
			"test_size": 40000,
			"number_of_epochs": 1000,
			"save_point_num_batches": 10,
			"model_name": "n4",
			"output_folder":"data/MIDGETNeural",
			"feature_file":"features.pkl"
		}
	},
	"MIDGET Neural[n5]":{
		"module_name":"genebench.diffmethods.MIDGET.neural",
		"class_name":"MIDGETNeural",
		"config": {
			"data_split": 0.25,
			"method_name":"",
			"batch_size": 40000,
			"test_size": 40000,
			"number_of_epochs": 1000,
			"save_point_num_batches": 10,
			"model_name": "n5",
			"output_folder":"data/MIDGETNeural",
			"feature_file":"features.pkl"
		}
	},
	"MIDGET Neural[n6]":{
		"module_name":"genebench.diffmethods.MIDGET.neural",
		"class_name":"MIDGETNeural",
		"config": {
			"data_split": 0.25,
			"method_name":"",
			"batch_size": 40000,
			"test_size": 40000,
			"number_of_epochs": 1000,
			"save_point_num_batches": 10,
			"model_name": "n6",
			"output_folder":"data/MIDGETNeural",
			"feature_file":"features.pkl"
		}
	},
	"MIDGET XGB[xgb1]":{
		"module_name":"genebench.diffmethods.MIDGET.xgb",
		"class_name":"MIDGETXgBoost",
		"config": {
			"feature_file":"features.pkl",
			"data_split": 0.25,
			"num_round": 1000,
			"model_name": "xgb1",
			"model_folder": "data/MIDGETXGB",
			"param": {
				"verbosity": 2,
				"nthread": 8,
				"max_depth":6, 
				"eta":0.2, 
				"objective":"binary:logistic"
			}
		}
	},
	"MIDGET XGB[xgb2]":{
		"module_name":"genebench.diffmethods.MIDGET.xgb",
		"class_name":"MIDGETXgBoost",
		"config": {
			"feature_file":"features.pkl",
			"data_split": 0.25,
			"num_round": 1000,
			"model_name": "xgb2",
			"model_folder": "data/MIDGETXGB",
			"param":{
				"verbosity": 2,
				"nthread": 8,
				"max_depth":12, 
				"lambda": 0.7,
				"eta":0.2, 
				"objective":"binary:logistic"
			}
		}
	}
},
```
- in **Benchmark** we configure the method groups for which we want to generate the benchmarks and evaluations by using the **method_groups** field. In the **runs** field we specify which validation set and data set to use with on the provided method groups. Below you may find an example config which runs tests on 4 method groups separating **silico data** from **transcription factor** data and **drug-gene data**.
```json
"Benchmark":{
	"method_groups":
		{
			"classic":{
				"name":"Transcription factor benchmark [classic]",
				"methods":[
					"TTest",
					"LIMMA",
					"SAM",
					"Random",
					"Characteristic direction py"
				]
			},
			"original_neural":{
				"name":"Transcription factor benchmark [original neural]",
				"methods":[
					"MIDGET Neural[n1]",
					"MIDGET Neural[n2]",
					"MIDGET Neural[n3]",
					"MIDGET Neural[n4]",
					"MIDGET Neural[n5]",
					"MIDGET Neural[n6]"
				]
			},
			"original_xgb":{
				"name":"Transcription factor benchmark [original xgb]",
				"methods":[
					"MIDGET XGB[xgb1]",
					"MIDGET XGB[xgb2]"
				]
			},
			"best_methods":{
				"name":"Transcription factor benchmark [original xgb]",
				"methods":[
					"Characteristic direction py",
					"MIDGET XGB[xgb2]",
					"MIDGET Neural[n2]",
					"MIDGET Neural[n3]",
					"LIMMA",
					"TTest"
				]
			}
	},
	"runs": [
		{
			"name":"Silico",
			"data_sources":["silico_linear"],
			"validation_sets":["silico"],
			"method_group_ids":[ "classic","original_neural","original_xgb", "best_methods"]
		},
		{
			"name":"Transcription Factor",
			"data_sources":["tf_perturbations"],
			"validation_sets":[
				"tf_encode",
				"tf_chea",
				"pp_interaction",
				"ppis_study",
				"all"],
			"method_group_ids":[ "classic","original_neural","original_xgb","best_methods"]
		},
		{
			"name":"Drugs Gene Interaction",
			"data_sources":["drug_perturbations"],
			"validation_sets":["drug_gene_interaction"],
			"method_group_ids":[ "classic","original_neural","original_xgb", "best_methods"]
		}
		
	]
}
```

### Providing custom implementations
- **Custom Method**
To use your own method with this framework you need to inherit from DiffMethod located in the diffmethods.base module

```python
# example for a random method implementation
import numpy as np
from genebench.datatypes import GeneDiffInput, GeneMethodResult
from genebench.diffmethods.base.diffmethod import DiffMethod


class CustomRandom(DiffMethod):
	def setup(self, config):
        # this code only runs once, used for configuration
		pass

    def train(self, input: GeneDiffInput):
		# this code is used to train your method
        pass

    def run(self, input: GeneDiffInput) -> GeneMethodResult:
		# this is the actual method
        rankings = np.random.rand(len(input.genes))
        genes, scores = self.sort(input.genes, rankings)
        return GeneMethodResult.from_separate_lists(genes, scores)
  ```
Then either add method in config specifying the correct module and class name.

- **Custom Metric**
To use your own method with this framework you need to inherit from Metric located in the evaluationmetrics.base module

```python
# example of custom metric implementation
from genebench.evaluationmetrics.base import Metric
from genebench.datatypes import GeneDiffValidation, GeneMethodResult
from genebench.utils import Utils
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import uniform
from sklearn.metrics import f1_score
import os


class CustomMetric(Metric):

    def __init__(self, config, output_folder):
        self.name = "F1"
        self.logger = Utils.get_logger("custom")
        self.config = config
        self.output_folder = output_folder
        self.custom_method = {}
        pass
	
	def __calculate_score(y_real, y_pred):
		return y_real - y_pred

    def add(self,
            pf,
            method_name,
            validation: GeneDiffValidation,
            result: GeneMethodResult):
		# this method is called for each experiment
        if validation.source not in self.custom_method:
            self.custom_method[validation.source] = {}

        custom_method = self.custom_method[validation.source]
        pf = pf.lower()
        self.logger.info(f"perturbation factor: {pf}")
        # we collect all results
        if method_name not in custom_method:
            custom_method[method_name] = {'y': [], 'pred': []}

        result_data = result.result
        valid_data = validation.data
        if pf not in valid_data:
            self.logger.error(f"{pf} not found in set {validation.source}")
            return

        genes = valid_data[pf]
        valid_genes = set([x.lower() for x in genes])
        real_class = []
        pred = []
        for gene_entry in result_data:
            gene_name = gene_entry.gene_name.lower()
            if gene_name in valid_genes:
                real_class.append(1)
            else:
                real_class.append(0)
            if gene_entry.score >= 0.5:
                pred.append(1)
            else:
                pred.append(0)
        custom_method[method_name]['y'].extend(real_class)
        custom_method[method_name]['pred'].extend(pred)

    def evaluate(self, group_name):
		# this method is called and aggregates all experiment results
		# to generate final score per data set
        Utils.create_folder_if_not_exist(self.output_folder)
        save_path = os.path.join(self.output_folder, self.name)
        Utils.create_folder_if_not_exist(save_path)

        for validation_source in self.method_f1.keys():
            custom_method = self.custom_method[validation_source]
            methods = []
            for method_name, _val in custom_method.items():
                pred = _val['pred']
                pred = np.nan_to_num(pred, True, 0.0, 1.0, 0.0)
                val = self.__calculate_score(_val['y'], pred)
                methods.append(f"{method_name} Custom Score: {f1:.4f}")
            name = f"{group_name}_{validation_source}"
            path_file = f"{name}.txt"
            path_method = os.path.join(save_path,
                                       path_file)
            with open(path_method, mode='wt', encoding='utf-8') as out_scores:
                out_scores.write('\n'.join(methods))
        self.custom_method = {}
        self.logger.info("done")
```

Then add metric in config specifying the correct module and class name.
Feel free to check the codebase in the evaluationmetrics folder to see example implementations

- **Custom silico generator**
To use your own silico generator with this framework you need to inherit from SilicoDataGenerator located in the silico.generators.base module:

```python
# example of custom silico data generator
from numpy.lib.function_base import diff
from genebench.silico.generators.base import SilicoDataGenerator
from genebench.datatypes import GeneData
from genebench.datatypes import GeneDiffValidation
from genebench.datatypes import GeoData
from genebench.utils import Utils
import numpy as np
import random

# this structure will contain the config from the "config" field
# of the generator entry 
class CustomDataGeneratorParam():
    def __init__(self,
                 diff_factor,
                 noise_factor,
                 num_replicates):
        self.diff_factor = diff_factor
        self.noise_factor = noise_factor
        self.num_replicates = num_replicates
        pass

# this example is the actual LinearDataGenerator
class CustomDataGenerator(SilicoDataGenerator):
    def __init__(self, config):
        super().__init__(config)
        self.param = LinearDataGeneratorParam(**self.config.params)
        self.logger = Utils.get_logger('LinearDataGenerator')
        pass

    def generate_single(self, validation_data, id, num_genes) -> GeoData:
        all_tfs = list(validation_data.data.keys())
        picked_tf = random.choice(all_tfs)
        perturbed_genes = set(validation_data.data[picked_tf])
        genes = Utils.get_random_gene_names(num_genes)
        mask = []
        for gene in genes:
            if gene in perturbed_genes:
                mask.append(1)
            else:
                mask.append(0)
        mask = np.array(mask)
        df_factor = self.param.diff_factor
        validation = []
        for index, mask_value in enumerate(mask.tolist()):
            if mask_value == 1:
                validation.append(genes[index])
        num_replicates = self.param.num_replicates
        mask = np.array([mask])
        mask = np.repeat(mask, num_replicates, axis=0).T
        control = np.random.rand(num_genes, num_replicates)
        effect = np.random.rand(num_genes, num_replicates) * df_factor
        perturbation = control + np.multiply(mask, effect)

        gene_data = GeoData({
            "name": f"CUSTOM_SIL_{id}",
            "perturbed_series_names": ['fakeseries'],
            "control_series_names": ['fakeseries'],
            "extra_info": {"none": "none"},
            "perturbed_array": perturbation.tolist(),
            "control_array": control.tolist(),
            "source": self.config.source,
            "genes": genes,
            "pf": picked_tf
        })

        return gene_data

    def generate_experiments(self,
                             validation_data: GeneDiffValidation,
                             num_experiments,
                             num_genes):
        experiments = []
        for id in range(0, num_experiments):
            experiments.append(self.generate_single(validation_data,
                                                    id,
                                                    num_genes))
        return experiments
```
- **Training MIDGET and benching**
  To do this you may check the train.py file in the examples folder. Note that training and benching may take a whole day with the default configurations.
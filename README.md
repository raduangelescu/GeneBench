![Image logo](https://raw.githubusercontent.com/raduangelescu/GeneBench/main/logo.svg)

##Overview

Gene Bench is a benchmark-ing framework used in analyzing methods that detect differentially expressed genes from biological samples.
Besides being a benchmarking framework it also contains some commonly used Differential genes detection algorithms:
  
  - [LIMMA](https://bioconductor.org/packages/release/bioc/html/limma.html) (R implementation)
  - [Characteristic Direction Py](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-79) (a python characteristic direction that is good for big data, implementation using Nipals PCA)
  - [Characteristic Direction](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-79) (R implementation **crashes on large data** with out of memory because of normal PCA )
  - [T-Test](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html)
  - [SAM](https://cran.r-project.org/web/packages/samr/index.html)
Our own methods based on machine learning:
  - MIDGET-xgb (2 variations)
  - MIDGET-neural (6 variations)
and a baseline random algorithm.

The framework also provides a base for calling methods from R (because most medical laboratories use R for statistical analysis). 

We provide already implemented evaluation metrics:
 - Kolmogorov based visual evaluation
 - F1 score
 - ROC/AUC curve analysis

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
## Installation
Install using pip 
```
pip install genebench
```

or install latest from github

```
git clone https://github.com/raduangelescu/GeneBench.git
python setup.py install

```

## Usage
The most important thing to note is that you need to have a correct config json file that describes what you are trying to do. The next sections will explain the code and config needed to access the framework features.

### Storage config
Any benchmark run via this framework needs storage. You may use any of the storage providers implemented in the framework, or both. 
The config file needs the **"Storage"** section in order to know how you will set and read the data. The below config is an example of using both providers.
```
	"Storage":{
		"providers":{
			"mongo": {
				"host":"insert the ip of your mongodb server here",
				"port": 27017,
				"anon":true, -> set this to true if you have anonymous login
				"user":"username for mongodb",
				"password":"password for mongodb",
				"database_name":"your database name",
				"validation_collection_name":"the validation collection name",
				"geo_data_collection_name":"geo data collection name",
				"results_collection_name":"results collection name"
			},
			"filesystem": {
				"base_path":"base folder for storing all info that must exist in filesystem", 
				"validation_folder":"the validation folder name",
				"geo_folder": "geo folder name",
				"results_folder": "results folder name"
			}
		},
		"load_order":[ -> order in which system searches for data
			"filesystem",
			"mongo"
		]
	},
```
If you only need to use one provider just remove the other from the above config form "providers" and "load_order". In most cases, when you don't have any mongodb installed you will probably only want to use the filesystem provider which should work locally.

### Generating silico data
The silico data generator tries to mimick a biological experiment: it will generate fake genes and fake perturbation factors and then create links between them. Based on the above biological network it will use silico providers to generate the fake data.

To generate the silico data you need to have the **"SilicoData"** and **"SilicoGenerators"** config sections in your config.json file
```
	"SilicoData":{
		"num_genes":10000, -> number of genes in experiments
		"num_experiments":100, -> number of experiments
		"num_pfs":100 -> number of perturbation factors
	},
	"SilicoGenerators":{
		"silico.linear":{ -> add an entry for every silico data generator
			"source": "silico_linear", -> name of the generator source for later identification
			"module_name": "genebench.silico.generators.linear", -> module name from which to import the generator
			"class_name":"LinearDataGenerator", -> the generator class name
			"params": { -> custom generator parameters
						"diff_factor":3.0,
						"noise_factor":0.5,
						"num_replicates":3
			}
		}
	},
```

After providing the correct config file you only need to run the below code:
```
from genebench.generatesilico import GenerateSilicoData

generate = GenerateSilicoData("config.json")
generate.run()
```
The generator run will populate the data in the storage provider. The data itself will be in a generic format, the only way to know it was generated by a silico provider is the source name from the config. **Be sure to add an unique source name for all silico generators so you can identify experiments later**


### Importing in vivo GEO data

For this feature you need the **GEOImporter** section in the config file:
```
	"GEOImporter":{
		"data_path":"../data", -> base folder to store temp data
		"download_retry_count": 10, -> retry download 10 times before failing (sometimes GEO will fail because of overload)
		"download_wait_seconds_before_retry": 5, -> if GEO download fails, wait 5 seconds
		"geo_cache_folder_name":"cache", -> the folder to store the cache 
		"input_data":[
		 	{
				"pf_field": "tf", -> field than marks the perturbation factor
				"name":"tf_perturbations", -> source name in storage
				"file":"GEOdatasets_tf_perturbations.json" -> input file with hand picked GEO experiments
	 		},
	 		{
				"pf_field": "drug",
		 		"name":"drug_perturbations",
		 		"file":"GEOdatasets_drugs_perturbations.json"
	 		}
		],
		"labeling": -> this is used to automatically filter GEO data in control/perturbation categories
		{
			"control":["control","Control", "wild type", "wildtype","Wild Type","Wildtype","wild-type","baseline", "control virus","uninduced", "untreated", "0 pM", "water", "Time of treatment 0 weeks", "none","pretreatment", "vehicle", "control virus"],
			"type":["genotype/variation","disease state", "protocol","growth protocol","dose", "agent", "description", "time"],
			"gene_names":"IDENTIFIER",
			"no_column_control":["MCF7/BUS, 0 pM","baseline", "minus", "Resistant", "untreated", "Untreated", "Control", "control", "nonresponder", "resistant", "vehicle", "saline","0 pM E2","DMS", "AO ", " AG", "_Veh_","N1+2","N3+4","N5+6"],
			"no_column_title": "title",
			"no_column_accession": "geo_accession"
		}
	},
```
It is usually the case that you will not modify this config in a normal run so you can paste it from above, removing comments.
To generate and store the experiment data you just run the below code, using whatever file-name your config is to replace 'config.json'.
```
from genebench.geoimporter import GEOImporter
importer = GEOImporter('config.json')
importer.importGEOData()
```

### Importing validation data
To be able to test if the methods you create/benchmark are better on the set of data, we need to have some validation data. To do this we provided multiple validation sources and a way to import them. 
The config section needed for this feature is **"ValidationDataImporter"**

```
	"ValidationDataImporter":{
		"base_path":"../data/validation_sources",
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
```

And the coded to import the data:

```
from genebench.validationdataimporter import ValidationDataImporter
importer = ValidationDataImporter('config.json')
importer.importValidationData()
```

### Running benchmarks
After you have all your validation and experiment data in your storage providers (general format) you will now configure the benchmarks you want to run via the **"AccuracyMetrics"**, **"Methods"** and **"Benchmark"** config sections:
- In **AccuracyMetrics** we configure the metrics we want to generate for all benchmarks. Below we provide an example for using the default provided ones
```
"AccuracyMetrics":{
    "output_folder":"../data/metrics" ,
    "metrics":{
        "Kolmogorov": {
            "name": "Kolmogorov",
            "module_name": "evaluationmetrics.kolmogorov",
            "class_name":"Kolmogorov",
            "params":{
            }
        },
        "ROC": {
            "name": "ROC",
            "module_name": "evaluationmetrics.roc",
            "class_name":"ROC",
            "params":{}
        },
        "F1": {
            "name": "F1",
            "module_name": "evaluationmetrics.f1",
            "class_name":"F1",
            "params":{}
        }
    }
},
```
-  In **Methods** we add all the methods we want to use for benchmarking: Each entry has: a key by which the program will address the code, module_name: the module in which the method relies, class_name: the actual method class name and a custom config field so you can use parameters in your own method implementation. Below we provide a config that accounts for all our currently supported methods:
```
	"Methods":{
		"TTest":{
			"module_name":"diffmethods.ttest",
			"class_name":"TTest",
			"config":{}
		},
		"LIMMA":{
			"module_name":"diffmethods.limma",
			"class_name":"LIMMA",
			"config":{}
		},
		"SAM":{
			"module_name":"diffmethods.sam",
			"class_name":"SAM",
			"config":{}
		},
		"Random":{
			"module_name":"diffmethods.random",
			"class_name":"Random",
			"config": {}
		},
		"Characteristic direction py":{
			"module_name":"diffmethods.chdirpy",
			"class_name":"ChDirPy",
			"config":{}
		},
		"MIDGET Neural[n1]":{
			"module_name":"diffmethods.MIDGET.neural",
			"class_name":"MIDGETNeural",
			"config": {
				"data_split": 0.25,
				"method_name":"",
				"batch_size": 40000,
				"test_size": 40000,
				"number_of_epochs": 1000,
				"save_point_num_batches": 10,
				"model_name": "n1",
				"output_folder":"../data/MIDGETNeural" ,
				"feature_file":"features.pkl"
			}
		},
		"MIDGET Neural[n2]":{
			"module_name":"diffmethods.MIDGET.neural",
			"class_name":"MIDGETNeural",
			"config": {
				"data_split": 0.25,
				"method_name":"",
				"batch_size": 40000,
				"test_size": 40000,
				"number_of_epochs": 1000,
				"save_point_num_batches": 10,
				"model_name": "n2",
				"output_folder":"../data/MIDGETNeural",
				"feature_file":"features.pkl"
			}
		},
		"MIDGET Neural[n3]":{
			"module_name":"diffmethods.MIDGET.neural",
			"class_name":"MIDGETNeural",
			"config": {
				"data_split": 0.25,
				"method_name":"",
				"batch_size": 40000,
				"test_size": 40000,
				"number_of_epochs": 1000,
				"save_point_num_batches": 10,
				"model_name": "n3",
				"output_folder":"../data/MIDGETNeural",
				"feature_file":"features.pkl"
			}
		},
		"MIDGET Neural[n4]":{
			"module_name":"diffmethods.MIDGET.neural",
			"class_name":"MIDGETNeural",
			"config": {
				"data_split": 0.25,
				"method_name":"",
				"batch_size": 40000,
				"test_size": 40000,
				"number_of_epochs": 1000,
				"save_point_num_batches": 10,
				"model_name": "n4",
				"output_folder":"../data/MIDGETNeural",
				"feature_file":"features.pkl"
			}
		},
		"MIDGET Neural[n5]":{
			"module_name":"diffmethods.MIDGET.neural",
			"class_name":"MIDGETNeural",
			"config": {
				"data_split": 0.25,
				"method_name":"",
				"batch_size": 40000,
				"test_size": 40000,
				"number_of_epochs": 1000,
				"save_point_num_batches": 10,
				"model_name": "n5",
				"output_folder":"../data/MIDGETNeural",
				"feature_file":"features.pkl"
			}
		},
		"MIDGET Neural[n6]":{
			"module_name":"diffmethods.MIDGET.neural",
			"class_name":"MIDGETNeural",
			"config": {
				"data_split": 0.25,
				"method_name":"",
				"batch_size": 40000,
				"test_size": 40000,
				"number_of_epochs": 1000,
				"save_point_num_batches": 10,
				"model_name": "n6",
				"output_folder":"../data/MIDGETNeural",
				"feature_file":"features.pkl"
			}
		},
		"MIDGET XGB[xgb1]":{
			"module_name":"diffmethods.MIDGET.xgb",
			"class_name":"MIDGETXgBoost",
			"config": {
				"feature_file":"features.pkl",
				"data_split": 0.25,
				"num_round": 1000,
				"model_name": "xgb1",
				"model_folder": "../data/MIDGETXGB",
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
			"module_name":"diffmethods.MIDGET.xgb",
			"class_name":"MIDGETXgBoost",
			"config": {
				"feature_file":"features.pkl",
				"data_split": 0.25,
				"num_round": 1000,
				"model_name": "xgb2",
				"model_folder": "../data/MIDGETXGB",
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
```
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
	},
```

### Providing custom implementations
- **Custom Method**
To use your own method with this framework you need to inherit from DiffMethod located in the diffmethods.base module, or have the following functions in your class body:

```
def setup(self, config):
    #setup code, only done once per benchmark
    pass

def train(self, input: GeneDiffInput):
    #code used for training, if necessary
    pass

def run(self, input: GeneDiffInput) -> GeneMethodResult:
    #actual method run that receives GeneDiffInput and outputs GeneMethodResult
    pass
  ```
Then either add method in config specifying the correct module and class name.

- **Custom Metric**
To use your own method with this framework you need to inherit from Metric located in the evaluationmetrics.base module, or have the following functions in your class body:

```
def __init__(self, output_folder):
    pass

def add(self,
        pf,
        method_name,
        validation: GeneDiffValidation,
        result: GeneMethodResult):
    #function in which you collect method results
    pass

def evaluate(self, group_name):
    #function in which you run the evaluation
    pass
```

Then add metric in config specifying the correct module and class name.
Feel free to check the codebase in the evaluationmetrics folder to see example implementations

- **Custom silico generator**
To use your own silico generator with this framework you need to inherit from SilicoDataGenerator located in the silico.generators.base module, or have the following functions in your class body:

```
def __init__(self, config):
    pass
def generate_experiments(self,
                            validation_data,
                            num_experiments,
                            num_genes) -> GeneData:
    pass
```

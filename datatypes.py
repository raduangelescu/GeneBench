from json import JSONEncoder
import numpy as np
import random


class GeneDiffValidation:
    def __init__(self, json=None):
        if json:
            self.data = json['data']
            self.source = json['source']
        else:
            self.data = []
            self.source = ""

    def to_dict(self):
        return {
            "data": self.data,
            "source": self.source
        }


class GeneListEntry:
    def __init__(self, gene_name, score):
        self.gene_name = gene_name
        self.score = score


class GeneMethodResult:
    def __init__(self, config=None):
        if config is None:
            self.result = []
            self.valid = False
        else:
            self.valid = True
            result_list = config['results']
            self.result = []
            for result in result_list:
                self.result.append(GeneListEntry(
                    **result
                ))

    def to_dict(self):
        array = []
        for x in self.result:
            array.append(x.__dict__)

        ret_dict = {'results': array}
        return ret_dict

    @staticmethod
    def from_pair_list(pair_list):
        new_object = GeneMethodResult()
        for x in pair_list:
            entry = GeneListEntry(x[0], x[1])
            new_object.result.append(entry)
        new_object.valid = True
        return new_object

    @staticmethod
    def from_separate_lists(genes, scores):
        new_object = GeneMethodResult()
        for idx, gene in enumerate(genes):
            entry = GeneListEntry(gene, scores[idx])
            new_object.result.append(entry)
        new_object.valid = True
        return new_object


class GeneData:
    def __init__(self):
        self.gene_input = GeneDiffInput()
        self.gene_validation = GeneDiffValidation()
        self.name = ""


class CustomDataEncoder(JSONEncoder):
    def default(self, o):
        return o.__dict__


class GeoData:
    def __init__(self,
                 config):
        self.name = config['name']
        self.perturbed_series_names = config['perturbed_series_names']
        self.control_series_names = config['control_series_names']
        self.extra_info = config['extra_info']
        self.perturbed_array = config['perturbed_array']
        self.control_array = config['control_array']
        self.source = config['source']
        self.genes = config['genes']
        self.pf = config['pf']

    def get_as_dict(self):
        data = self.get_meta_data()
        data.update(self.get_big_data())
        return data

    def get_big_data(self):
        return {
            'perturbed_array': self.perturbed_array,
            'control_array': self.control_array
        }

    def get_meta_data(self):
        return {
             'name': self.name,
             'perturbed_series_names': self.perturbed_series_names,
             'control_series_names': self.control_series_names,
             'extra_info': self.extra_info,
             'source': self.source,
             'genes': self.genes,
             'pf': self.pf
        }


class GeneDiffInput:
    def __init__(self):
        self.control = np.array([])
        self.perturbed = np.array([])
        self.genes = set([])

    @staticmethod
    def from_geo_data(data: GeoData):
        input = GeneDiffInput()
        # limit data to most important because we don't have
        # enough memory to run the R methods with full data 
        if len(data.control_array) < len(data.genes):
            # some experiments need to be transposed
            control = np.array(data.control_array).T.tolist()
            perturbed = np.array(data.perturbed_array).T.tolist()
        else:
            control = data.control_array
            perturbed = data.perturbed_array
        control = [x[:6] for x in control]
        perturbed = [x[-6:] for x in perturbed]

        input.control = np.array(control)
        input.perturbed = np.array(perturbed)
        input.genes = data.genes

        for idx, gene_name in enumerate(input.genes):
            if not isinstance(gene_name, str):
                input.genes[idx] = f"fake_{random.randint(0,9999)}"

        return input

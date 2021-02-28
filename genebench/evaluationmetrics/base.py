from genebench.datatypes import GeneDiffValidation, GeneMethodResult


class Metric():
    def __init__(self, output_folder):
        pass

    def add(self,
            pf,
            method_name,
            validation: GeneDiffValidation,
            result: GeneMethodResult):

        pass

    def evaluate(self, group_name):
        pass

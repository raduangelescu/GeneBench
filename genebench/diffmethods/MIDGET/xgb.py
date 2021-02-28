import os
import xgboost as xgb
from sklearn.metrics import accuracy_score
import random
from genebench.utils import Utils
from genebench.datatypes import GeneDiffInput, GeneMethodResult
from genebench.diffmethods.MIDGET.common import MIDGET

class MIDGETXgBoostConfig:
    def __init__(self,
                 feature_file,
                 data_split,
                 num_round,
                 param,
                 model_name,
                 model_folder):
        self.num_round = num_round
        self.model_name = model_name
        self.model_folder = model_folder
        self.param = param
        self.feature_file = feature_file
        self.data_split = data_split


class MIDGETXgBoost(MIDGET):

    def setup(self, config):
        self.config = MIDGETXgBoostConfig(**config)
        logger_name = f"MIDGETXgBoost[{self.config.model_name}]"
        self.logger = Utils.get_logger(logger_name)   
        model_path = os.path.join(self.config.model_folder,
                                  self.config.model_name,
                                  'model.json')
        self.logger.info(f"pir: Loading model: {model_path}")
        if os.path.isfile(model_path):
            bst = xgb.Booster(self.config.param)
            bst.load_model(model_path)
            self.model = bst
        else:
            self.logger.warning(f"No model located in {model_path}")
        pass

    def build_model(self):
        pass

    def train(self):
        features_file_name = self.config.feature_file
        self.logger.info("started xgb training ")
        data = self.load_training(features_file_name)
        random.shuffle(data)
        split_point = int(len(data) * self.config.data_split)
        test_data_split = data[:split_point]
        train_data_split = data[split_point:]
        num_round = self.config.num_round
        param = self.config.param
        self.logger.info(f"params for training {param}, num_rounds: {num_round}")
        X_train_split, y_train_split = self.get_data(train_data_split)
        train_data = xgb.DMatrix(data=X_train_split,
                                 label=y_train_split)
        X_test_split, y_test_split = self.get_data(test_data_split)
        test_data = xgb.DMatrix(data=X_test_split,
                                label=y_test_split)
        self.logger.info("Finished preparing data, starting training")
        bst = xgb.train(param, train_data, num_round)
        self.logger.info("Finished training, starting prediction eval")
        y_pred = bst.predict(test_data)
        predictions = [round(value) for value in y_pred]
        accuracy = accuracy_score(y_test_split, predictions)
        self.logger.info(f'Prediction accuracy: {accuracy}, saving')
        model_path = os.path.join(self.config.model_folder,
                                  self.config.model_name)
        Utils.create_folder_if_not_exist(model_path)
        bst.save_model(os.path.join(model_path, "model.json"))

    def run(self, input: GeneDiffInput) -> GeneMethodResult:
        A = input.control
        B = input.perturbed
        genes = input.genes
        self.logger.info('preparing feature vectors')
        X = self.do_run_feature_vectors(A, B, genes)
        self.logger.info('doing predictions')
        X = xgb.DMatrix(X)
        scores = self.model.predict(X)
        genes, scores = self.sort(genes, scores)
        scores = [x.item() for x in scores]
        return GeneMethodResult.from_separate_lists(genes, scores)

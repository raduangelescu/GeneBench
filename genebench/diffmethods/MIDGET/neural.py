import numpy as np
import os
import tensorflow as tf
import random
from genebench.datatypes import GeneDiffInput, GeneMethodResult
from genebench.diffmethods.MIDGET.common import MIDGET
from genebench.utils import Utils

class MIDGETNeuralConfig:
    def __init__(self,
                 feature_file,
                 method_name,
                 data_split,
                 batch_size,
                 test_size,
                 number_of_epochs,
                 save_point_num_batches,
                 model_name,
                 output_folder):
        self.method_name = method_name
        self.batch_size = batch_size
        self.test_size = test_size
        self.number_of_epochs = number_of_epochs
        self.save_point_num_batches = save_point_num_batches
        self.model_name = model_name
        self.output_folder = output_folder
        self.feature_file = feature_file
        self.data_split = data_split


class MIDGETNeural(MIDGET):

    def setup(self, config):
        self.config = MIDGETNeuralConfig(**config)
        logger_name = f"MIDGETNeural[{self.config.model_name}]"
        self.logger = Utils.get_logger(logger_name)

        model_path = os.path.join(self.config.output_folder,
                                  self.config.model_name,
                                  "model")
        self.logger.info(f"Loading model: {model_path}")
        self.model = tf.keras.models.load_model(model_path)

    def build_model(self):
        pass

    def train(self):
        model_path = os.path.join(self.config.output_folder,
                                  self.config.model_name,
                                  "structure.json")

        with open(model_path) as f:
            self.model = tf.keras.models.model_from_json(f.read())
            self.model.compile(loss='binary_crossentropy',
                               optimizer='adam',
                               metrics=['accuracy'])
        feature_file = self.config.feature_file
        self.logger.info(f"loading training data form {feature_file}")
        data = self.load_training(feature_file)
        random.shuffle(data)
        self.logger.info('Building model')
        batch_size = self.config.batch_size
        number_of_epochs = self.config.number_of_epochs
        split_point = int(len(data) * self.config.data_split)
        test_data = data[:split_point]
        train_data = data[split_point:]

        X, y = self.get_data(train_data)

        self.model.fit(X,
                       y,
                       epochs=number_of_epochs,
                       batch_size=batch_size,
                       verbose=1)
        msg = f"accuracy on random data with {len(test_data)} dimension"
        self.logger.info(msg)
        X, y = self.get_data(test_data)
        results = self.model.evaluate(X,
                                      y,
                                      batch_size=batch_size)
        self.logger.info(f'evaluate loss and accuracy {results}')
        self.model.save(self.config.model_name)
        self.logger.info("Done training")

    def run(self, input: GeneDiffInput) -> GeneMethodResult:
        A = input.control
        B = input.perturbed
        genes = input.genes
        self.logger.info("preparing feature vectors")
        X = self.do_run_feature_vectors(A, B, genes)
        self.logger.info("doing predictions")
        shape = self.model.predict(X)
        shape = np.sum(shape, axis=1)
        genes, scores = self.sort(genes, shape)
        scores = [x.item() for x in scores]
        return GeneMethodResult.from_separate_lists(genes, scores)

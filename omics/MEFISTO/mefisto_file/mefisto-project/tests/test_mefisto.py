import unittest
import pandas as pd
from src.data.data_loader import load_data
from src.models.mefisto import MEFISTO

class TestMEFISTO(unittest.TestCase):

    def setUp(self):
        self.data = load_data('data/evodevo.csv')
        self.model = MEFISTO(input_dim=self.data.shape[1], latent_dim=2)

    def test_model_initialization(self):
        self.assertIsNotNone(self.model)
        self.assertEqual(self.model.latent_dim, 2)
        self.assertEqual(self.model.input_dim, self.data.shape[1])

    def test_forward_pass(self):
        output = self.model.forward(self.data)
        self.assertEqual(output.shape, (self.data.shape[0], self.model.latent_dim))

    def test_loss_function(self):
        loss = self.model.loss(self.data)
        self.assertIsInstance(loss, float)

if __name__ == '__main__':
    unittest.main()
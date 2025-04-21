import pandas as pd
import torch
from torch import nn, optim
from src.data.data_loader import load_data
from src.models.mefisto import MEFISTO
from src.training.loss import total_loss

class Trainer:
    def __init__(self, model, data, epochs=100, learning_rate=0.001):
        self.model = model
        self.data = data
        self.epochs = epochs
        self.learning_rate = learning_rate
        self.optimizer = optim.Adam(self.model.parameters(), lr=self.learning_rate)
        self.loss_history = []

    def train(self):
        for epoch in range(self.epochs):
            self.model.train()
            self.optimizer.zero_grad()

            # Forward pass
            reconstructed, latent_factors = self.model(self.data)

            # Compute loss
            loss = total_loss(reconstructed, self.data, latent_factors)
            self.loss_history.append(loss.item())

            # Backward pass and optimization
            loss.backward()
            self.optimizer.step()

            if epoch % 10 == 0:
                print(f'Epoch [{epoch}/{self.epochs}], Loss: {loss.item():.4f}')

        return self.loss_history

def main():
    # Load data
    data = load_data('data/evodevo.csv')
    
    # Initialize model
    model = MEFISTO(input_dim=data.shape[1], latent_dim=2)

    # Create trainer
    trainer = Trainer(model, data)

    # Train the model
    loss_history = trainer.train()

if __name__ == "__main__":
    main()
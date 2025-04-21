import pandas as pd
import torch
from src.data.data_loader import load_data
from src.models.mefisto import MEFISTO
from src.training.trainer import Trainer
from src.visualization.plots import plot_results

def main():
    # Load the dataset
    data = load_data('data/evodevo.csv')
    
    # Initialize the MEFISTO model
    model = MEFISTO(input_dim=data.shape[1])
    
    # Set up the trainer
    trainer = Trainer(model=model, data=data)
    
    # Train the model
    loss_history, recon_history, gp_history = trainer.train(epochs=1000)
    
    # Visualize the results
    plot_results(loss_history, recon_history, gp_history)

if __name__ == "__main__":
    main()
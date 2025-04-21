import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def plot_loss_convergence(loss_history, recon_history, gp_history):
    plt.figure(figsize=(12, 5))
    plt.subplot(1, 2, 1)
    plt.plot(loss_history, label='Total Loss')
    plt.plot(recon_history, label='Reconstruction Loss')
    plt.plot(gp_history, label='GP Prior Loss')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.legend()
    plt.title('Loss Convergence')

def plot_latent_factors(t, Z_learned):
    plt.subplot(1, 2, 2)
    plt.scatter(t, Z_learned[:, 0], color='blue', label='Latent Factor 1')
    plt.scatter(t, Z_learned[:, 1], color='red', label='Latent Factor 2')
    plt.xlabel('Time')
    plt.ylabel('Learned Latent Factor Values')
    plt.legend()
    plt.title('Learned Latent Factors vs. Time')

def visualize_results(loss_history, recon_history, gp_history, t, Z):
    Z_learned = Z.detach().cpu().numpy()
    plot_loss_convergence(loss_history, recon_history, gp_history)
    plot_latent_factors(t, Z_learned)
    plt.tight_layout()
    plt.show()
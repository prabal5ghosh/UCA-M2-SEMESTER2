import numpy as np
import torch
import torch.nn as nn

class GaussianProcessPrior:
    def __init__(self, length_scale=1.0, variance=1.0):
        self.length_scale = length_scale
        self.variance = variance

    def kernel(self, x1, x2):
        """Squared Exponential Kernel"""
        sqdist = torch.sum(x1**2, 1).view(-1, 1) + torch.sum(x2**2, 1) - 2 * torch.mm(x1, x2.t())
        return self.variance * torch.exp(-0.5 / self.length_scale**2 * sqdist)

    def compute_prior(self, Z):
        """Compute the Gaussian Process prior for latent factors Z"""
        K = self.kernel(Z, Z) + 1e-6 * torch.eye(Z.size(0))  # Add small noise for numerical stability
        return K

    def compute_prior_loss(self, Z):
        """Compute the GP prior loss"""
        K = self.compute_prior(Z)
        K_inv = torch.inverse(K)
        prior_loss = 0.5 * torch.mm(Z.t(), torch.mm(K_inv, Z)) + 0.5 * torch.logdet(K) + Z.size(0) * 0.5 * np.log(2 * np.pi)
        return prior_loss.squeeze()
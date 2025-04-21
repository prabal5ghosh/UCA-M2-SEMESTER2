def total_loss(recon_loss, gp_prior_loss, lambda_gp=1.0):
    return recon_loss + lambda_gp * gp_prior_loss

def reconstruction_loss(y_true, y_pred):
    return ((y_true - y_pred) ** 2).mean()

def gp_prior_loss(Z, kernel_func):
    # Assuming Z is the latent representation and kernel_func is a callable kernel function
    K = kernel_func(Z)
    return -torch.logdet(K + 1e-6 * torch.eye(K.size(0)))  # Adding a small value for numerical stability
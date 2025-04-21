def squared_exponential_kernel(x1, x2, length_scale=1.0, variance=1.0):
    """
    Computes the squared exponential (RBF) kernel between two inputs.

    Parameters:
    - x1: First input (numpy array).
    - x2: Second input (numpy array).
    - length_scale: Length scale parameter of the kernel.
    - variance: Variance parameter of the kernel.

    Returns:
    - Kernel value between x1 and x2.
    """
    sqdist = np.sum(x1**2, 1).reshape(-1, 1) + np.sum(x2**2, 1) - 2 * np.dot(x1, x2.T)
    return variance * np.exp(-0.5 / length_scale**2 * sqdist)

def kernel_matrix(X, length_scale=1.0, variance=1.0):
    """
    Computes the kernel matrix for a given dataset.

    Parameters:
    - X: Input dataset (numpy array).
    - length_scale: Length scale parameter of the kernel.
    - variance: Variance parameter of the kernel.

    Returns:
    - Kernel matrix (numpy array).
    """
    return squared_exponential_kernel(X, X, length_scale, variance)
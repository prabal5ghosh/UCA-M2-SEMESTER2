import pandas as pd

def load_data(filepath):
    """
    Load the evodevo dataset from a CSV file.

    Parameters:
    filepath (str): Path to the CSV file.

    Returns:
    pd.DataFrame: Loaded dataset as a pandas DataFrame.
    """
    data = pd.read_csv(filepath)
    return data

def preprocess_data(data):
    """
    Preprocess the dataset for training the MEFISTO model.

    Parameters:
    data (pd.DataFrame): The raw dataset.

    Returns:
    pd.DataFrame: Preprocessed dataset.
    """
    # Example preprocessing steps (customize as needed)
    # Remove rows with missing values
    data = data.dropna()
    
    # Normalize or scale features if necessary
    # data = (data - data.mean()) / data.std()

    return data

def get_data(filepath):
    """
    Load and preprocess the dataset.

    Parameters:
    filepath (str): Path to the CSV file.

    Returns:
    pd.DataFrame: Preprocessed dataset ready for training.
    """
    raw_data = load_data(filepath)
    processed_data = preprocess_data(raw_data)
    return processed_data
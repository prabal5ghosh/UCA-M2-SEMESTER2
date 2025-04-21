import pytest
import pandas as pd
from src.data.data_loader import load_data

def test_load_data():
    # Test loading the dataset
    data = load_data('data/evodevo.csv')
    assert isinstance(data, pd.DataFrame), "Loaded data should be a pandas DataFrame"
    assert not data.empty, "Loaded data should not be empty"
    assert 'column_name' in data.columns, "DataFrame should contain the expected column"

def test_data_shape():
    # Test the shape of the loaded data
    data = load_data('data/evodevo.csv')
    expected_shape = (100, 5)  # Replace with the actual expected shape
    assert data.shape == expected_shape, f"Expected shape {expected_shape}, but got {data.shape}"
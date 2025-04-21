# MEFISTO Project

This project implements the MEFISTO algorithm from scratch using the `evodevo.csv` dataset. The MEFISTO model is designed for analyzing and modeling complex biological data through latent factor analysis.

## Project Structure

```
mefisto-project
├── data
│   └── evodevo.csv
├── src
│   ├── __init__.py
│   ├── data
│   │   ├── __init__.py
│   │   └── data_loader.py
│   ├── models
│   │   ├── __init__.py
│   │   ├── gp_prior.py
│   │   ├── kernels.py
│   │   └── mefisto.py
│   ├── training
│   │   ├── __init__.py
│   │   ├── loss.py
│   │   └── trainer.py
│   └── visualization
│       ├── __init__.py
│       └── plots.py
├── notebooks
│   └── mefisto_analysis.ipynb
├── scripts
│   └── train_mefisto.py
├── tests
│   ├── __init__.py
│   ├── test_data_loader.py
│   ├── test_gp_prior.py
│   └── test_mefisto.py
├── requirements.txt
├── setup.py
└── README.md
```

## Dataset

The dataset used for training the MEFISTO algorithm is located in the `data/evodevo.csv` file. This dataset contains biological data that will be processed and analyzed using the MEFISTO model.

## Installation

To set up the project, clone the repository and install the required dependencies:

```bash
git clone <repository-url>
cd mefisto-project
pip install -r requirements.txt
```

## Usage

To train the MEFISTO model, run the following script:

```bash
python scripts/train_mefisto.py
```

This script will load the dataset, initialize the MEFISTO model, and start the training process. The results will be saved and can be visualized using the provided Jupyter notebook.

## Visualization

After training, you can analyze the results using the Jupyter notebook located in the `notebooks` directory:

```bash
jupyter notebook notebooks/mefisto_analysis.ipynb
```

## Testing

Unit tests are provided to ensure the functionality of the data loading, GP prior, and MEFISTO model implementations. To run the tests, use:

```bash
pytest tests/
```

## Contributing

Contributions are welcome! Please open an issue or submit a pull request for any improvements or bug fixes.

## License

This project is licensed under the MIT License. See the LICENSE file for more details.
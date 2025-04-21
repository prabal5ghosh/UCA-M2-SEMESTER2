from setuptools import setup, find_packages

setup(
    name='mefisto-project',
    version='0.1.0',
    author='Your Name',
    author_email='your.email@example.com',
    description='Implementation of the MEFISTO algorithm using the evodevo dataset.',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    install_requires=[
        'numpy',
        'pandas',
        'matplotlib',
        'scikit-learn',
        'torch',
        'gpytorch',
        'jupyter'
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)
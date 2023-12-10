# Systematically comparing feature-based representations of intra-regional and inter-regional brain dynamics

This repository contains code to accompany our preprint, "Systematically comparing feature-based representations of intra-regional and inter-regional brain dynamics".
Users may follow this repo to reproduce all analyses and visualizations contained in the preprint -- broadly, this includes extracting time-series features from functional magnetic resonance imaging (fMRI) data to serve as the basis for a series of linear support vector machine (SVM) classifiers for case--control comparisons.

All code is a mixture of R and python, with shell scripts used to execute operations.
Please note that the repository is structured such that some steps are designed to be run on a high-performance computing (HPC) cluster for parallelization with a PBS job scheduler.

## Prerequisites

Before running the feature extraction and classification analysis, make sure you have the following packages installed:

Python
- pyspi
- sklearn
- numpy
- pandas
- pyarrow

R
- tidyverse
- theft
- correctR
- R.matlab

The following packages are required to reproduce figures from the preprint:

Python
- nibabel
- nilearn

R
- glue
- icesTAF
- patchwork
- ggseg
- ggsegHO
- ggsegDefaultExtra
- LaCroixColoR
- ggraph
- igraph
- dendextend
- cowplot
- ggsignif
- ggpubr
- feather
- reticulate
- see
- colorspace
- broom
- scales
- ggridges


## Installation

1. Clone this repository to your local machine:

    ```bash
    git clone https://github.com/DynamicsAndNeuralSystems/fMRI_FeaturesDisorders.git
    ```

2. Install the required Python packages:

    ```bash
    pip install -r requirements.txt
    ```

## Usage

1. Place your BOLD fMRI data files in the `data` directory.

2. Modify the configuration file `config.yaml` to specify the case-control comparisons and other parameters.

3. Run the analysis script:

    ```bash
    python analysis.py
    ```

4. The extracted time-series features will be saved in the `output` directory.

## Results

The results of the analysis will be saved in the `results` directory. Include any relevant information about the results here.

## Contributing

If you would like to contribute to this project, please follow the guidelines in [CONTRIBUTING.md](CONTRIBUTING.md).

## License

This project is licensed under the [MIT License](LICENSE).

## Contact

If you have any questions or need further assistance, please contact [annie.bryant@sydney.edu.au](mailto:annie.bryant@sydney.edu.aum).

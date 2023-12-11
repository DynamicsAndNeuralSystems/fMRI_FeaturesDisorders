# Systematically comparing feature-based representations of intra-regional and inter-regional brain dynamics

This repository contains code to accompany our preprint, "Systematically comparing feature-based representations of intra-regional and inter-regional brain dynamics".
Users may follow this repo to reproduce all analyses and visualizations contained in the preprint -- broadly, this includes extracting time-series features from functional magnetic resonance imaging (fMRI) data to serve as the basis for a series of linear support vector machine (SVM) classifiers for case--control comparisons.

All code is a mixture of R and python, with shell scripts used to execute operations.
Please note that the repository is structured such that some steps are designed to be run on a high-performance computing (HPC) cluster for parallelization with a PBS job scheduler.



## Installation

First, clone this repository to your local machine:

    ```bash
    git clone https://github.com/DynamicsAndNeuralSystems/fMRI_FeaturesDisorders.git
    ```

Parts of this project are in R (v4.3.0), and parts are in Python (v3.9.0).
We recommend that users create a conda environment for this project, which can be accomplished as follows:

    ```bash
    conda create -n fMRI_FeaturesDisorders python=3.9.0
    conda activate fMRI_FeaturesDisorders
    ```

From this conda environment, install the required Python packages:

    ```bash
    pip install -r python_requirements.txt
    ```

R packages required for feature extraction and classification analysis can be installed from `R_requirements.txt` as follows:

    ```bash
    while IFS=" " read -r package version; 
    do 
        Rscript -e "devtools::install_version('"$package"', version='"$version"')"; 
    done < "R_requirements.txt"
    ```

### Data visualization

Additional R packages will be needed to reproduce visualizations for the manuscript in `generate_all_figures.ipynb`, which can be installed with the following:

    ```bash
    while IFS=" " read -r package version; 
    do 
        Rscript -e "devtools::install_version('"$package"', version='"$version"')"; 
    done < "R_requirements.txt"
    ```

# Usage

## Preparing data

## Extracting time-series features

## Performing case--control classification

## Data visualization

All figures included in the preprint can be recreated using the Jupyter notebook `generate_all_figures.ipynb`.


## Results

The results of the analysis will be saved in the `results` directory. Include any relevant information about the results here.

## Contributing

If you would like to contribute to this project, please follow the guidelines in [CONTRIBUTING.md](CONTRIBUTING.md).

## License

This project is licensed under the [MIT License](LICENSE).


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
- glue
- R.matlab
- icesTAF

The following packages are required to reproduce figures from the preprint:

Python
- nibabel
- nilearn

R
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


## Contact

If you have any questions or need further assistance, please contact [annie.bryant@sydney.edu.au](mailto:annie.bryant@sydney.edu.aum).

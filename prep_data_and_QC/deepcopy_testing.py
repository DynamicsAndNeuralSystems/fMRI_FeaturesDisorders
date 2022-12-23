# Import the rest of the modules
from pyspi.calculator import Calculator
from pyspi.data import Data
from copy import deepcopy

# Config file
config_file="/headnode1/abry4213/github/fMRI_FeaturesDisorders/prep_data_and_QC/pyspi14_config.yaml"

# Define basecalcs
basecalc = Calculator()
# basecalc_config = Calculator(configfile=config_file)

# Create deepcopies
print("Creating default deepcopy")
basecalc_copy = deepcopy(basecalc)
# print("Creating custom config deepcopy")
# basecalc_config_copy = deepcopy(basecalc_config)
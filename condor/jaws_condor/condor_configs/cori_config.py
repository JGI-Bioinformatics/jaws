# Commands and parameters for Cori
from jaws_condor.config import Config

CORI_NORMAL_C = "#SBATCH -C haswell"
CORI_EXVIVO_C = "#SBATCH -C skylake"

CONFIG_CORI = Config()

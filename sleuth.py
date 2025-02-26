# Imports
from imports import *
from constants import *

class Sleuth:
    """
    Abstract class which contains all the required functions to run the Sleuth Rscript. 
    """
    @abstractmethod
    def run_sleuth():
        """
        Runs the sleuth.R script. 
        """
        os.system(f"Rscript sleuth.R >> processed_data/output_log/output.txt 2>&1")
    
    @abstractmethod
    def print_log():
        """
        Prints the required output for this step in the pipeline (class).
        """
        # Simply takes the output file written by sleuth.R and copies it to the log file. 
        with open(f"{Constants.LOG_FOLDER}/{Constants.LOG_FILE}", "a") as log_file:
            with open(f"{Constants.PROCESSED_FOLDER}/sleuth/sleuth_output", "r") as sleuth_output:
                log_file.write(f"--Sleuth Differentially Expressed Genes--\n{sleuth_output.read()}")
            log_file.write("\n")
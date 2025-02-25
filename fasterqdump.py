# Imports
from imports import *
from constants import *

class FasterqDump:
    """
    Abstract class which can perform fasterq-dump on a directory.
    """
    @abstractmethod
    def fasterq_dump(input_dir, output_dir):
        """
        Runs fasterq_dump on a group of files in a directory.
        
        :parameter input_dir (str): The directory containing our input files.
        :parameter output_dir (str): The directory where we would like to place our output file. 
        """
        files = os.listdir(input_dir)
        for file in files:
            print(f"{input_dir}/{file}")
            os.system(f"fasterq-dump {input_dir}/{file} -o {output_dir}/{file} >> processed_data/output_log/output.txt 2>&1")
            
    
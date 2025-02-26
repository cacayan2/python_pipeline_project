# Imports
from imports import *
from constants import *

class Helpers:
    """An abstract class containing a number of helper methods."""
    @abstractmethod
    def check_arg(args=None):
        """
        Function to parse command line arguments.
        """
        parser = argparse.ArgumentParser(description="ADD TITLE OF SCRIPT HERE (shows on help)")
        parser.add_argument("-r", "--reference",
        help="reference transcriptome",
        required=True)
        parser.add_argument("-c", "--condensed",
        help="condensed pipeline", 
        required=True)
        return parser.parse_args(args)

    @abstractmethod
    def clean(string):
        """
        Function to clean a string being read into python.

        :parameter string (str): The string to be cleaned. 
        """
        return string.replace("\n", "").strip()
    
    @abstractmethod
    def refresh_directory():
        """
        Function which creates a new empty processed_data folder.
        """
        # Deleting the processed_data folder and its contents if it exists
        # to prevent pollution of the data. 
        if os.path.isdir("processed_data"):
            os.system("rm -r processed_data")
        
        # Performing the same above to create the pipeline log file.
        if os.path.isdir(f"{Constants.LOG_FOLDER}"):
            os.system(f"rm -r {Constants.LOG_FOLDER}")

        # Creating and moving into the parent directory
        os.mkdir("processed_data")
        os.chdir("processed_data")
        
        # Creating the subdirectories.
        for dir in Constants.SUBDIRECTORIES:
            os.mkdir(dir)

        # Moving back into the parent directory and populating the pipeline log directory.
        os.chdir("..")
        os.mkdir(Constants.LOG_FOLDER)
        os.chdir(Constants.LOG_FOLDER)
        log = open(Constants.LOG_FILE, "x")
        log.close()
        os.chdir("..")
    
    @abstractmethod 
    def obtain_SRR(SRR: str):
        """
        Downloads the normalized SRR file given an SRR accession number.

        :parameter SRR: The SRR accession number to download. 
        """
        os.system(f"wget -q --show-progress -P {Constants.PROCESSED_FOLDER}/input_data https://sra-pub-run-odp.s3.amazonaws.com/sra/{SRR}/{SRR}")

    @abstractmethod
    def med_list(list):
        """
        Calculates the median value in a list.

        :parameter list (list): The list to calculate the median value out of. 
        :return float: The median of the list. 
        """
        sorted_list = sorted(list)
        n = len(sorted_list)
        if n % 2 == 0:
            mid1 = sorted_list[n // 2 - 1]
            mid2 = sorted_list[n // 2]
            return (mid1 + mid2) / 2
        else:
            return sorted_list[n // 2]
           
    @abstractmethod
    def mean_list(list):
        """
        Calculates the mean of a list.

        :parameter list (list): The list to calculate the mean from.
        :return float: The mean of the list. 
        """
        if not list:
            return 0
        total = sum(list)
        count = len(list)
        return total / count


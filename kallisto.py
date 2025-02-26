# Imports
from imports import *
from constants import *
from helpers import *

class Kallisto:
    """
    Class which contains all the functionality required to create a kallisto index and run TPM quantification.

    :attribute _input_dir (str): The directory which contains the fastq files to be passed to kallisto.
    :attribute _srr_list (list): The list of SRR accession numbers to be used in this analysis.
    :attribute _files (list): The list of files contained in _input_dir.
    :attribute _conditions (dict): A dictionary containing each SRR number and their respective experimental conditions. 
    """
    def __init__(self, srr_list, input_dir):
        """
        Instantiator function for a Kallisto object.
        :parameter srr_list (list): The list of SRR accession numbers to be used in this analysis. 
        :parameter input_dir (str): The directory which contains the fastq files to be passed to kallisto. 
        """
        # Setting the object variables. 
        self._input_dir = input_dir
        self._srr_list = srr_list
        self._files = os.listdir(input_dir)
        self._conditions = dict()
        
        # For each SRR number, we assign a condition (if the SRR is associated with the study of this assignment, they will be given a special condition,
        # otherwise they will be assigned a generic condition).
        for i, file in enumerate(self._srr_list):
            if (file == "SRR5660030" or file == "SRR5660030.fastq"):
                self._conditions[file] = {"Patient": "Donor 1", "Condition": "2dpi"}
                continue
            elif (file == "SRR5660033" or file == "SRR5660033.fastq"):
                self._conditions[file] = {"Patient": "Donor 1", "Condition": "6dpi"}
                continue
            elif (file == "SRR5660044" or file == "SRR5660044.fastq"):
                self._conditions[file] = {"Patient": "Donor 3", "Condition": "2dpi"}
                continue
            elif (file == "SRR5660045" or file == "SRR5660045.fastq"):
                self._conditions[file] = {"Patient": "Donor 3", "Condition": "6dpi"}
                continue
            else:
                self._conditions[file] = {"Condition": f"Condition {i}"}
    
    def quantification(self):
        """
        Runs kallisto quantification using the object variables for each SRR. 
        """
        for srr in self._srr_list:
            os.system(f"kallisto quant -i {Constants.PROCESSED_FOLDER}/kallisto/index.idx -o {Constants.PROCESSED_FOLDER}/kallisto/{srr} -b 30 -t 2 {self._input_dir}/{srr}_1.fastq {self._input_dir}/{srr}_2.fastq >> processed_data/output_log/output.txt 2>&1")

    def print_log(self):
        """
        Prints the required output to the log folder for this step in the pipeline (class).
        """
        # Write the header section of this print_log().
        with open(f"{Constants.LOG_FOLDER}/{Constants.LOG_FILE}", "a") as log:
            log.write(f"--Kallisto Quantification Statistics--\nsample\tcondition\tmin_tpm\tmed_tpm\tmean_tpm\tmax_tpm\n")

        # For each SRR number, calculate each of the required metrics and write them to the log file. 
        for srr in self._srr_list:
            with open(f"{Constants.PROCESSED_FOLDER}/kallisto/{srr}/abundance.tsv", "r") as abundance:
                reader = csv.reader(abundance, delimiter = "\t")
                next(reader)
                condition = self._conditions[srr]["Condition"]
                sample = srr
                tpm = [float(row[4]) for row in reader if row]
                min_tpm = min(tpm)
                med_tpm = Helpers.med_list(tpm)
                mean_tpm = Helpers.mean_list(tpm)
                max_tpm = max(tpm)

                with open(f"{Constants.LOG_FOLDER}/{Constants.LOG_FILE}", "a") as log:
                    log.write(f"{sample}\t{condition}\t{min_tpm}\t{med_tpm}\t{mean_tpm}\t{max_tpm}\n")
        with open(f"{Constants.LOG_FOLDER}/{Constants.LOG_FILE}", "a") as log:
            log.write("\n")

    def create_sleuth_input(self):
        """
        Creates the input files required for Sleuth analysis. 
        """
        # Opens the sleuth input file. 
        with open(f"{Constants.PROCESSED_FOLDER}/sleuth/sleuth_input", "w") as sleuth_input:
            # Writes the header to the log file. 
            sleuth_input.write("sample\tcondition\tpath\n")

            # Goes through each of the SRR's and writes the sleuth data to an output file. 
            for srr in self._conditions:
                path = f"{Constants.PROCESSED_FOLDER}/kallisto/{srr}"
                sleuth_input.write(f"{str(srr)}\t{self._conditions[srr]['Condition']}\t{str(path)}\n")



                

                        
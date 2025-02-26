# Imports
from imports import *
from constants import *

class Spades:
    """
    Class which defines the functionality to perform SPAdes assembly.

    :attribute _srr_list (str): The list of SRR's that are going to be assembled using SPAdes.
    :attribute _spades_cmd (str): The command that will be sent to the SPAdes assembler.
    :attribute _conditions (dict): The dictionary containing the conditions associated with a given SRR entry. 
    """
    def __init__(self, srr_list):
        """
        Instantiator method for a Spades object.

        :param srr_list (list): The list of SRR entries to be assembled using Spades. 
        """
        
        # Presetting the class variables. 
        self._srr_list = srr_list
        self._spades_cmd = ""
        self._conditions = dict()

        # Iterating through the SRR list and setting the condition.
        # Because this is done on a pre-annotated study, conditions for specific
        # SRR entries are hard-coded into this logic.
        # If the entered SRR entries are not the ones associated with the study, 
        # then a generic condition will be assigned. 
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

    def run_spades(self):
        """
        Function which runs the SPAdes assembler on the reads.
        """
        # This SPAdes command is the beginning of any of the SPAdes commands we will send to os.system().
        spades_cmd = "spades.py -k 77 -t 2 --only-assembler"
        spades_cmd1 = ""
        spades_cmd2 = ""

        # We set the path to the reads using path_to_srr (this contains the reads mapped to only by bowtie2).
        # We also set a count for donor1 and donor3 - this helps when generating the SPAdes command. 
        path_to_srr = f"{Constants.PROCESSED_FOLDER}/bowtie2/mapped_unzipped"
        donor1_count = 0
        donor3_count = 0

        # We iterate through each of the SRR's in the list. 
        for i, srr in enumerate(self._srr_list):
            # If any of the SRR's are the ones we are studying (HCMV study), we run a special command to ensure that the paired-end reads are paired correctly. 
            if (srr == "SRR5660030" 
            or srr == "SRR5660033"
            or srr == "SRR5660044"
            or srr == "SRR5660045"):
                if self._conditions[srr]["Patient"] == "Donor 1":
                    donor1_count += 1
                    spades_cmd1 = f"{spades_cmd1} --pe-1 {donor1_count} {path_to_srr}/{srr}_mapped_1.fq --pe-2 {donor1_count} {path_to_srr}/{srr}_mapped_2.fq"
                else:
                    donor3_count += 1
                    spades_cmd2 = f"{spades_cmd2} --pe-1 {donor3_count} {path_to_srr}/{srr}_mapped_1.fq --pe-2 {donor3_count} {path_to_srr}/{srr}_mapped_2.fq"
            else:     
                spades_cmd = f"{spades_cmd} --pe-1 {i + 1} {path_to_srr}/{srr}_mapped_1.fq --pe-2 {i + 1} {path_to_srr}/{srr}_mapped_2.fq"
        
        # This checks whether we found that our set of SRR's were the ones of interest mentioned previously. 
        if spades_cmd == "spades.py -k 77 -t 2 --only-assembler":
            spades_cmd1 = f"{spades_cmd1} -o {Constants.PROCESSED_FOLDER}/spades/Donor1"
            spades_cmd2 = f"{spades_cmd2} -o {Constants.PROCESSED_FOLDER}/spades/Donor3"
            self._spades_cmd = f"{spades_cmd}{spades_cmd1}\n{spades_cmd}{spades_cmd2}\n"
            os.system(f"{spades_cmd}{spades_cmd1} >> processed_data/output_log/output.txt")
            os.system(f"{spades_cmd}{spades_cmd2} >> processed_data/output_log/output.txt")
        else:         
            spades_cmd = f"{spades_cmd} -o {Constants.PROCESSED_FOLDER}/spades >> processed_data/output_log/output.txt 2>&1"
            self._spades_cmd = spades_cmd
            os.system(spades_cmd)

    def print_log(self):
        """
        Prints the appropriate output for this step in the pipeline (class).
        """
        with open(f"{Constants.LOG_FOLDER}/{Constants.LOG_FILE}", "a") as log:
            log.write(f"--Assembly Using SPAdes--\n{self._spades_cmd}\n")
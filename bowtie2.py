# Imports
from imports import *
from constants import *
from fasterqdump import *

class Bowtie2:
    """
    Class which contains all the functionality required to perform Bowtie2 alignment.

    :attribute _ref_accession (str): The accession number for the reference transcriptome.
    :attribute _ref_assembly (str): The assembly number for the reference transcriptome.
    :attribute _ref_assembly_name (str): The project associated with the reference transcriptome.
    :attribute _ref_assembly_path (str): The path to the assembly .fna file.
    :attribute _ref_organism_name (str): The name of the organism.
    :attribute _srr_list (str): The list of SRR accession numbers to be used in the analysis.
    :attribute _fasterq_dir (str): The directory of fastq files to be used.
    :attribute _conditions (list): A dictionary of conditions associated with an SRR. 
    """
    def __init__(self, ref_accession, srr_list, fasterq_dir):
        """
        Instantiator function for Bowtie2 objects. 

        :parameter ref_accession (str): The accession number for the reference transcriptome.
        :parameter srr_list (str): The list of SRR accession numbers to be used in the analysis.
        :parameter fasterq_dir (str): The directory of fastq files to be used. 
        """
        # Setting class variables. 
        self._ref_accession = ref_accession
        self._ref_assembly = ""
        self._ref_assemnbly_name = ""
        self._ref_assembly_path = ""
        self._ref_organism_name = ""
        self._srr_list = srr_list
        self._fasterq_dir = fasterq_dir

        # Iterating through the SRR list and setting the condition.
        # Because this is done on a pre-annotated study, conditions for specific
        # SRR entries are hard-coded into this logic.
        # If the entered SRR entries are not the ones associated with the study, 
        # then a generic condition will be assigned. 
        self._conditions = dict()
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

    def retrieve_reference(self):
        """
        Downloads the reference genome assembly from NCBI. 
        """
        # This may seem like a lot, but all this does is:
        #   1. Sets an appropriate name for output files.
        #   2. Takes the reference genome number, obtains the assembly accession number, downloads it, then unzips it. 
        Entrez.email = "ecacayan@luc.edu"
        handle = Entrez.esearch(db = "assembly", term = self._ref_accession)
        record = Entrez.read(handle)
        handle = Entrez.esummary(db = "assembly", id = record["IdList"])
        summary = Entrez.read(handle)
        handle.close()
        self._ref_assembly = summary["DocumentSummarySet"]["DocumentSummary"][0]["AssemblyAccession"]
        self._ref_assembly_name = summary["DocumentSummarySet"]["DocumentSummary"][0]["AssemblyName"]
        os.system(f"datasets download genome accession {self._ref_assembly} --include gff3,rna,cds,protein,genome,seq-report --filename {Constants.PROCESSED_FOLDER}/bowtie2/{self._ref_accession}_dataset.zip >> processed_data/output_log/output.txt 2>&1")
        os.system(f"unzip {Constants.PROCESSED_FOLDER}/bowtie2/{self._ref_accession}_dataset.zip -d {Constants.PROCESSED_FOLDER}/bowtie2 >> processed_data/output_log/output.txt 2>&1")
        self._ref_assembly_path = f"{Constants.PROCESSED_FOLDER}/bowtie2/ncbi_dataset/data/{self._ref_assembly}/{self._ref_assembly}_{self._ref_assembly_name}_genomic.fna"
    
    def bowtie2_build(self):
        """
        Builds the Bowtie2 index. 
        """
        # First we obtain the organism name from NCBI in order to create a name for our output file. 
        Entrez.email = "ecacayan@luc.edu"
        handle = Entrez.esearch(db = "assembly", term = self._ref_accession)
        record = Entrez.read(handle)
        handle = Entrez.esummary(db = "assembly", id = record["IdList"])
        summary = Entrez.read(handle)
        self._ref_organism_name = f"{summary['DocumentSummarySet']['DocumentSummary'][0]['Organism'].replace(' ', '').replace('(', '').replace(')', '')}"

        # We then take the assembly that we downloaded and use it to build our Bowtie2 genome index. 
        os.system(f"bowtie2-build {self._ref_assembly_path} ./{Constants.PROCESSED_FOLDER}/bowtie2/{self._ref_organism_name} >> processed_data/output_log/output.txt 2>&1")

    def run_bowtie2(self):
        """
        Runs Bowtie2 mapping. 
        """
        # First we create all the directories we need to store our output. 
        os.system(f"mkdir -p {Constants.PROCESSED_FOLDER}/bowtie2/output")
        os.system(f"mkdir -p {Constants.PROCESSED_FOLDER}/bowtie2/mapped_compressed")
        os.system(f"mkdir -p {Constants.PROCESSED_FOLDER}/bowtie2/mapped_unzipped")

        # For each SRR, we then run bowtie2 mapping with our index using the paired-end reads and unzip the file. 
        for srr in self._srr_list:
            os.system(f"bowtie2 -x {Constants.PROCESSED_FOLDER}/bowtie2/{self._ref_organism_name} -1 {self._fasterq_dir}/{srr}_1.fastq -2 {self._fasterq_dir}/{srr}_2.fastq -S {Constants.PROCESSED_FOLDER}/bowtie2/output/{srr}.sam --al-conc-gz {Constants.PROCESSED_FOLDER}/bowtie2/mapped_compressed/{srr}_mapped_%.fq.gz >> processed_data/output_log/output.txt 2>&1")
        for file in os.listdir(f"{Constants.PROCESSED_FOLDER}/bowtie2/mapped_compressed"):
            os.system(f"gunzip -c {Constants.PROCESSED_FOLDER}/bowtie2/mapped_compressed/{file} > {Constants.PROCESSED_FOLDER}/bowtie2/mapped_unzipped/{file.replace('.gz', '')}")

    
    def print_log(self):
        """
        Writes the appropriate output for this step in the pipeline (class) to the log file. 
        """
        # First we set the title. 
        with open(f"{Constants.LOG_FOLDER}/{Constants.LOG_FILE}", "a") as log:
            log.write("--Bowtie2 Filtering--\n")
        
        # We then iterate through each of the SRR's and obtain the condition. 
        for srr in self._srr_list:
            current_condition = ""
            if (srr == "SRR5660030" 
                or srr == "SRR5660033"
                or srr == "SRR5660044"
                or srr == "SRR5660045"):
                current_condition = f"{self._conditions[srr]['Patient']} ({self._conditions[srr]['Condition']})"
            else:
                current_condition = f"{self._conditions[srr]['Condition']}"
            
            # We then count how many reads are contained before mapping to the reference genome. 
            before_count = 0
            before_directory = self._fasterq_dir
            for file in os.listdir(before_directory):
                if srr in file:
                    with open(f"{before_directory}/{file}", "r") as current_srr:
                        for line in current_srr.readlines():
                            before_count += 1
            before_count = int(before_count / 8)

            # We then count how many reads are contained after mapping to the reference genome.
            after_count = 0
            after_directory = f"{Constants.PROCESSED_FOLDER}/bowtie2/mapped_unzipped"
            for file in os.listdir(after_directory):
                if srr in file:
                    with open(f"{after_directory}/{file}", "r") as current_srr:
                        for line in current_srr.readlines():
                            after_count += 1
            after_count = int(after_count / 8)

            # For each condition, we write how many read pairs are present before and after filtering. 
            with open(f"{Constants.LOG_FOLDER}/{Constants.LOG_FILE}", "a") as log:
                log.write(f"{current_condition} had {str(before_count)} read pairs before Bowtie2 filtering and {str(after_count)} read pairs after.\n")
        with open(f"{Constants.LOG_FOLDER}/{Constants.LOG_FILE}", "a") as log:
            log.write("\n")

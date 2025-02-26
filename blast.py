# Imports
from imports import *
from constants import *

class Blast:
    """
    Class which contains all the functionality required to run BLAST.

    :attribute _subfamily (str): The subfamily of the reference genome.
    :attribute _srr_list (str): The list of SRR's to use in the analysis.
    :attribute _conditions (str): The conditions associated with each SRR. 
    """
    def __init__(self, srr_list):
        """
        Instantiator method for a Blast object.

        :param srr_list: The list of SRR's to use in the analysis. 
        """
        # Setting class variables. 
        self._subfamily = ""
        self._srr_list = srr_list

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
        
    def makeblastdb(self, accession):
        """
        Creates the local database for blast.

        :parameter accession (str): The accesion number for the reference genome. 
        """
        # Obtain the taxonomy id from NCBI for the reference genome. 
        Entrez.email = "ecacayan@luc.edu"
        handle = Entrez.esearch(db = "nucleotide", term = accession)
        record = Entrez.read(handle)
        handle = Entrez.esummary(db = "nucleotide", id = record["IdList"])
        summary = Entrez.read(handle)
        tax_id = int(summary[0]["TaxId"])
        
        # Obtain the subfamily of the reference genome using the taxonomy id. 
        handle = Entrez.efetch(db = "taxonomy", id = tax_id, retmode = "xml")
        records = Entrez.read(handle)
        handle.close()
        for record in records[0]["LineageEx"]:
            if record["Rank"] == "subfamily":
                self._subfamily = record["ScientificName"].lower()
                break
        
        # Creating a database containing the subfamily of the reference genome. 
        os.system(f"datasets download virus genome taxon {self._subfamily} --filename {Constants.PROCESSED_FOLDER}/blast/ncbi_dataset.zip --include genome >> processed_data/output_log/output.txt 2>&1")
        os.system(f"unzip {Constants.PROCESSED_FOLDER}/blast/ncbi_dataset.zip -d {Constants.PROCESSED_FOLDER}/blast >> processed_data/output_log/output.txt 2>&1")
        os.system(f"makeblastdb -in {Constants.PROCESSED_FOLDER}/blast/ncbi_dataset/data/genomic.fna -out {Constants.PROCESSED_FOLDER}/blast/{self._subfamily} -title {self._subfamily} -dbtype nucl >> processed_data/output_log/output.txt 2>&1")

    def blast(self):
        """
        Runs blast on the query sequences against the local database. 
        """
        # Creates a new directory to contain the results from the blast. 
        os.system(f"mkdir {Constants.PROCESSED_FOLDER}/blast/results")

        # Creates a temporary file containing the longest contig in each assembly and runs blastn on that contig. 
        for dir in os.listdir(f"{Constants.PROCESSED_FOLDER}/spades"):
            file_path = f"{Constants.PROCESSED_FOLDER}/spades/{dir}/contigs.fasta"
            record = next(SeqIO.parse(file_path, "fasta"))
            longest_contig = f">{record.id}\n{record.seq}"

            with open(f"{Constants.PROCESSED_FOLDER}/blast/fasta_tmp.fasta", "w") as fasta_tmp:
                fasta_tmp.write(longest_contig)   
            os.system(f"blastn -max_hsps 1 -query {Constants.PROCESSED_FOLDER}/blast/fasta_tmp.fasta -db {Constants.PROCESSED_FOLDER}/blast/{self._subfamily} -out {Constants.PROCESSED_FOLDER}/blast/results/results{dir}.tsv -outfmt '6 sacc pident length qstart qend sstart send bitscore evalue stitle' >> processed_data/output_log/output.txt 2>&1")

    def print_log(self):
        """
        Prints the appropriate output for this step in the pipeline (class) to the pipeline log file. 
        """
        # We open the log file and create a title. 
        with open(f"{Constants.LOG_FOLDER}/{Constants.LOG_FILE}", "a") as log:
            log.write("--BLAST Alignment--\n")
            
            # For each of the subdirectories we find in the spades folder (which are named after the conditions), we write the top 10 results from the .tsv file generated from the blast. 
            for dir in sorted(os.listdir(f"{Constants.PROCESSED_FOLDER}/spades")):
                log.write(f"{dir}:\nsacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n")
                with open (f"{Constants.PROCESSED_FOLDER}/blast/results/results{dir}.tsv", "r") as curr_result:
                    lines = curr_result.readlines()
                    
                    # This checks the length of the file, if the length is less than 10 the whole file is included in the log file, otherwise we iterate through each line and count
                    # how many lines we pass, and if we reach 10 then we stop writing to the log file. 
                    if len(lines) <= 10:
                        for line in lines:
                            log.write(line)
                    else:
                        counter = 0
                        for line in lines:
                            if counter == 10:
                                break
                            log.write(line)
                            counter += 1
                log.write("\n")
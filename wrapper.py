# Imports
from imports import *
from helpers import *
from transcriptomeindex import *
from constants import *
from kallisto import *
from fasterqdump import *
from sleuth import *
from bowtie2 import *
from spades import *
from blast import *

def main():
    """
    Controls the execution flow of the pipeline. Please note that most of the os.system() calls in the functions below
    have their output silenced and instead written to the file processed_data/output_log (the processed_data folder
    is created the moment that wrapper.py is run).
    """
    # Reading the arguments into python. 
    arguments = Helpers.check_arg(sys.argv[1:])
    ref_accession = arguments.reference
    condensed = arguments.condensed
    
    # input_dir refers to the location of where the SRR files will be downloaded.
    # fasterq_dir refers to the location of where the fasterq-dump files are located.
    # fasterq_dir is more important as which directory is selected determines
    # whether the full dataset or the sample dataset will be used. 
    input_dir = ""
    fasterq_dir = ""
    srr_list = list()
    print("Output for all functions can be found in processed_data/output_log/output.txt")

    # Deletes and recreates the processed_data directory, where all the data files for the pipeline will be stored.
    print("Cleaning directory...")
    Helpers.refresh_directory()

    # Logic which determines from the user input whether the sample dataset will be used
    # or the full dataset. This is the branch of the logic which
    # allows the full dataset to be used. 
    if (condensed == "false" 
        or condensed == "False" 
        or condensed == "f" 
        or condensed == "F"):
        
        # Setting the appropriate directories for the complete dataset. 
        input_dir = f"{Constants.PROCESSED_FOLDER}/input_data"
        fasterq_dir = f"{Constants.PROCESSED_FOLDER}/fasterqdump"

        # Iterating through the input_srr file, extracting the SRR's to download from NCBI. 
        print("Downloading SRR files...")
        with open("input_srr", "r") as input_srr:
            lines = input_srr.readlines()
            for line in lines:
                srr_list.append(Helpers.clean(line))
                Helpers.obtain_SRR(Helpers.clean(line))

        # Unpacking the files using fasterq-dump. 
        print("Unpacking files using fasterq-dump...")
        FasterqDump.fasterq_dump(input_dir, fasterq_dir)
        
        # Obtaining the transcriptome index for kallisto.
        print("Generating transcriptome index...")
        transcriptome_index = TranscriptomeIndex()
        transcriptome_index.obtain_transcriptome(ref_accession)
        transcriptome_index.build_idx()
        transcriptome_index.print_log()

        # Building the index and quantifying TPM with kallisto. 
        print("Quantifying using Kallisto...")
        kallisto = Kallisto(srr_list, fasterq_dir)
        kallisto.quantification()
        kallisto.print_log()
        kallisto.create_sleuth_input()

        # Determining differentially expressed genes with Sleuth. 
        print("Analyzing differentially expressed genes with Sleuth...")
        Sleuth.run_sleuth()
        Sleuth.print_log()

        # Filtering reads that align to the reference transcriptome with Bowtie2.
        print("Filtering with Bowtie2...")
        bowtie2 = Bowtie2(ref_accession, srr_list, fasterq_dir)
        bowtie2.retrieve_reference()
        bowtie2.bowtie2_build()
        bowtie2.run_bowtie2()
        bowtie2.print_log()

        # Creating the assembly with SPAdes. 
        print("Generating assemblies with SPAdes...")
        spades = Spades(srr_list)
        spades.run_spades()
        spades.print_log()

        # Finding alignments with BLAST (generating the local database then running blastn).
        print("Generating alignments with BLAST...")
        blast = Blast(srr_list)
        blast.makeblastdb(ref_accession)
        blast.blast()
        blast.print_log()

        print("Finished pipeline!")

    # Logic which determines from the user input whether the sample dataset will be used
    # or the full dataset. This is the branch of the logic which
    # allows the partial dataset to be used. 
    elif (condensed == "true"
          or condensed == "True"
          or condensed == "t"
          or condensed == "T"):
        fasterq_dir = "sample_data"
        
        # Obtaining the transcriptome index for kallisto.
        print("Generating transcriptome index...")
        transcriptome_index = TranscriptomeIndex()
        transcriptome_index.obtain_transcriptome(ref_accession)
        transcriptome_index.build_idx()
        transcriptome_index.print_log()

        # Because the test dataset is pre-curated, the srr_list is hard-coded into python. 
        srr_list = ["SRR5660030", "SRR5660033", "SRR5660044", "SRR5660045"]

        # Building the index and quantifying TPM with kallisto. 
        print("Quantifying using Kallisto...")
        kallisto = Kallisto(srr_list, fasterq_dir)
        kallisto.quantification()
        kallisto.print_log()
        kallisto.create_sleuth_input()

        # Determining differentially expressed genes with Sleuth. 
        print("Analyzing differentially expressed genes with Sleuth...")
        Sleuth.run_sleuth()
        Sleuth.print_log()

        # Filtering reads that align to the reference transcriptome with Bowtie2.
        print("Filtering with Bowtie2...")
        bowtie2 = Bowtie2(ref_accession, srr_list, fasterq_dir)
        bowtie2.retrieve_reference()
        bowtie2.bowtie2_build()
        bowtie2.run_bowtie2()
        bowtie2.print_log()

        # Creating the assembly with SPAdes. 
        print("Generating assemblies with SPAdes...")
        spades = Spades(srr_list)
        spades.run_spades()
        spades.print_log()

        # Finding alignments with BLAST (generating the local database then running blastn).
        print("Generating alignments with BLAST...")
        blast = Blast(srr_list)
        blast.makeblastdb(ref_accession)
        blast.blast()
        blast.print_log()

        print("Finished pipeline!")

    # If an unexpected argument is detected for -c,
    # execution stops. 
    else:
        raise SystemExit("Not a correct argument for condensed - please enter true/false (t/f)!")

# Executes main().
if __name__ == "__main__":
    main()
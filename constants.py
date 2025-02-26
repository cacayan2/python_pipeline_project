class Constants:
    """
    A class of constants (mainly to store important directories/subdirectories).
    """
    # The folder to store the log file. 
    LOG_FOLDER = "PipelineProject"

    # The name of the log file. 
    LOG_FILE = "PipelineProject.log"

    # The name of the folder which contains all the through data handled by the pipeline. 
    PROCESSED_FOLDER = "processed_data"

    # The subdirectories contained in PROCESSED_FOLDER. 
    SUBDIRECTORIES =  subdirectories = ["blast", "bowtie2", "kallisto", "sleuth", "transcriptomeindex", "fasterqdump", "input_data", "spades", "output_log"]
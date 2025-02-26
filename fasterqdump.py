from imports import *
from constants import *

class FasterqDump:
    @abstractmethod
    def fasterq_dump(input_dir, output_dir):
        files = os.listdir(input_dir)
        for file in files:
            print(f"{input_dir}/{file}")
            os.system(f"fasterq-dump {input_dir}/{file} -o {output_dir}/{file} >> processed_data/output_log/output.txt 2>&1")
            
    
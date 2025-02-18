import os
import argparse
import sys
from abc import ABC, abstractmethod

class Helpers:
    """An abstract class containing a number of helper methods."""
    @abstractmethod
    def check_arg(args=None):
        """Function to parse command line arguments."""
        parser = argparse.ArgumentParser(description="ADD TITLE OF SCRIPT HERE (shows on help)")
        parser.add_argument("-i", "--input",
        help="input directory",
        required=True)
        return parser.parse_args(args)

    @abstractmethod
    def clean(string):
        """Function to clean a string being read into python."""
        return string.replace("\n", "").strip()

# Setting nothing for the input
arguments = Helpers.check_arg(sys.argv[1:])
input_data_dir = arguments.input

# Checking if the test_data or input_data directory exists.
if not os.path.isdir(input_data_dir):
    raise SystemExit("Input data directory not found: please specify a directory that exists. Stopping execution.")

# Clearing any directories that might be present in the 
## utils.py
import gzip
from os import PathLike
import pandas as pd
import logging

# Set up basic configuration for logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Define the columns in the GFF3 file
GFF3_COLS=["Sequence ID","source","Feature Type","Feature Start","Feature End","Score","Strand","Phase","INFO"]
# Define the columns in the input file
INPUT_COLS=["Gene Symbol","Transcript ID","Genome Build","Protein Domains"]

class Utils:
    @staticmethod
    def read_input_file(file_path: str) -> pd.DataFrame:
        """
        Reads the input file containing protein data and returns it as a pandas DataFrame.

        Parameters:
        file_path (str): The path to the input file.

        Returns:
        pd.DataFrame: The protein data as a pandas DataFrame.
        """
        try:
            protein_data = pd.read_csv(file_path, sep='\t')
            for col in INPUT_COLS:
                if col not in protein_data.columns:
                    raise ValueError(f"The input file {file_path} does not contain the required columns {INPUT_COLS}.")
            return protein_data
        except Exception as e:
            logging.error(f"file reading error for {file_path}")
            raise e

    @staticmethod
    def read_feature_file(file_path:  PathLike ) -> pd.DataFrame:
        """
        Reads the parquet feature file and returns its contents as a dictionary.

        Parameters:
        file_path (str): The path to the feature file.

        Returns:
        dict: The features data as a dictionary.
        """
        try:
            feature_data = pd.read_parquet(file_path)
            return feature_data
        except FileNotFoundError:
            logging.error(f"Error: The file {file_path} was not found. Please check the file path.")
            raise FileNotFoundError
        except Exception as e:
            logging.error(f"An error occurred while reading the feature file: {e}")
            raise e

    @staticmethod
    def read_gff3_to_df(file_path: PathLike) -> pd.DataFrame:
        """
        Read a GFF3 file and convert it into a pandas DataFrame.

        Args:
            file_path (PathLike): The path to the GFF3 file.

        Raises:
            ValueError: If the file is not a recognized file type.

        Returns:
            pd.DataFrame: A DataFrame containing the parsed GFF3 data.
        """
       
        try:
            # Check if the file is gzipped or not
            if file_path.suffix==".gz":
                FILE=gzip.open(file_path, 'r')
            elif file_path.suffix==".gff3" or file_path.suffix==".gff":
                FILE=open(file_path, 'r')
            else:
                raise ValueError(f"The file {file_path} is not a recognized file type. Please provide a .gz, .gff3 or .gff file.")
                
            with FILE as file:
                outlist=list()
                for line in file.readlines():
                    line=line.decode().strip()
                    if not line.startswith("#"):
                        # Parse the line and create a dictionary
                        parsed_dict=dict(zip(GFF3_COLS,line.split("\t")))
                        # Check if the number of columns in the line matches the number of columns in GFF3_COLS
                        if len(parsed_dict) != len(GFF3_COLS):
                            raise ValueError(f"The number of columns in the line does not match the number of columns in GFF3_COLS.")
                        # Parse the INFO field and create a dictionary
                        INFO_dict=dict([item.split("=") for item in parsed_dict["INFO"].split(";")])
                        # Update the parsed dictionary with the INFO dictionary
                        parsed_dict.update(INFO_dict)
                        parsed_dict.pop("INFO")
                        outlist.append(parsed_dict)

            return pd.DataFrame(outlist)
        except FileNotFoundError:
            logging.error(f"Error: The file {file_path} was not found. Please check the file path.")
            return pd.DataFrame()
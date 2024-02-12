## main.py
import argparse
import logging

import pandas as pd

from importlib_resources import files
from pathlib import Path
from os import PathLike

from protein2genome.converter import Converter
from protein2genome.utils import Utils

def parse_arguments():
    parser = argparse.ArgumentParser(description='Convert protein coordinates to genomic coordinates.')

    try:
        default_feature_file_path = Path(files('protein2genome').joinpath('resources').resolve()) / "gencode_CDS_DataFrame.parquet"
    except (ModuleNotFoundError, FileNotFoundError):
        default_feature_file_path = Path('protein2genome/resources/genecode_CDS_DataFrame.parquet')
        logging.info(f"Default feature file not found by importlib. Using {default_feature_file_path} as default.")

    parser.add_argument('-f', '--feature_file_path',  help='Path to the feature file containing genomic data. '
                        'It should be a parquet file. Default: resources/genecode_CDS_DataFrame.parquet.',
                        default=default_feature_file_path)
    parser.add_argument('-o', '--output_file_path', help='Path to the output file. If not provided, the output will be '
                        'written to output.tsv.', default="output.tsv")
    parser.add_argument('input_file_path',
                        help="Path to the input file containing protein data. It should be a tab-separated file and "
                             "contains the following columns: Gene Symbol, Transcript ID, Genome Build, Protein Domains.")

    args = parser.parse_args()
    return args

def protein2genome(protein_data: pd.DataFrame, feature_file_path: PathLike):
    """
    function to convert protein coordinates to genomic coordinates.

    Args:
        input_file_path (str): Path to the input file.
        feature_file_path (str): Path to the feature file.
        output_file_path (str, optional): Path to the output file. Defaults to "output.tsv".
    """
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    # Convert protein coordinates to genomic coordinates
    covnerter = Converter(feature_file_path)
    logging.info(f"Feature file {feature_file_path} read successfully.")
    
    logging.info("Converting protein coordinates to genomic coordinates...")
    output_df=covnerter.convert(protein_data)

    # Write output file
    return output_df
    
    
def main():
    args=parse_arguments()
    # Read input file
    protein_data = Utils.read_input_file(args.input_file_path)
    logging.info(f"Input file {args.input_file_path} read successfully.")
    output_df=protein2genome(protein_data, args.feature_file_path)
    output_df.to_csv(args.output_file_path, sep='\t', index=False)
    logging.info(f"Output written to { args.output_file_path}.")
    
    
if __name__ == "__main__":
    main()


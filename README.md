# MS Basecalling

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

## Description

This is the demo python library for converting protein coordinates to genomic coordinates.
It's only a demo for my coding skill, and not supposed to be used in real world research and production.

## Scripts
- 'converter.py' is the main Converter Class which has the core function to convert the amino acid position to genomic postion.
- 'utils.py' is the Class which contains the static methodes for IO purpose, it also contains a read_gff3_to_df function which can read gff3 files as pandas DataFrame.
- '\__main__.py' is the main entry point of the program.
- 'prototype.ipynb' is the scratch notebook I used to build up the core function and prepare the feature file.

## Resource File
The resource files are stored in `protein2genome/resources`.  The `gencode_CDS_DataFrame.parquet` is the feature file of CDS dataframes extracted from the [gencode gff3 files](https://www.gencodegenes.org/human/).
Code can be found in prototype.ipynb.

## Installation

Run `python setup.py install` to install this package.

## Test

Run `python .\tests\test_protein2genome.py` to run the unit test.

## Usage
It can also run as `python -m protein2genome -o example/test_out.tsv example/teststest_input.tsv` without installation.
Or you can run `protein2genome` command after installtion of the package.

## License

This project is licensed under the [MIT License](LICENSE).

## Contact

- Email: jiawku17@gmail.com
- LinkedIn: [Jiawei Gu](https://www.linkedin.com/in/jiawei-gu/)

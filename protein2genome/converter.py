## converter.py
import logging
import pandas as pd

from os import PathLike
from tqdm import tqdm
from protein2genome.utils import Utils

class Converter:
    def __init__(self, feature_file_path: PathLike):
        """
        Initializes the Converter object.

        Args:
            feature_file_path (PathLike): The path to the feature file.

        Returns:
            None
        """
        self.cds_df = Utils.read_feature_file(feature_file_path)

    def convert(self, protein_data: pd.DataFrame) -> pd.DataFrame:
        """
        Converts protein data to genomic coordinates.

        Args:
            protein_data (str): Protein data to be converted.

        Returns:
            pd.DataFrame: DataFrame containing the converted genomic coordinates.
        """

        output_list=list()
        for _, row in tqdm(protein_data.iterrows(), total=protein_data.shape[0]):
            for item in row["Protein Domains"].split(";"):
                try:
                    domain_dict=dict()
                    domain_dict["Domain Name"],domain_dict["protein_pos"]=item.split(":")
                    domain_dict["protein_start"],domain_dict["protein_end"]=list(domain_dict["protein_pos"].replace("—","-").strip().split("-"))
                    domain_dict.pop("protein_pos")

                    # Map protein region to genomic coordinates
                    genomic_loc=self.map_protein_to_genomic(row["Transcript ID"],int(domain_dict["protein_start"]),int(domain_dict["protein_end"]),row["Genome Build"])
                    
                    domain_dict["Gene Symbol"]=row["Gene Symbol"]
                    domain_dict["Genome_Build"]=row["Genome Build"]
                    domain_dict["Chrom"]=genomic_loc[0]
                    domain_dict["Domain Coordinates"]=domain_dict["protein_start"]+"-"+domain_dict["protein_end"]
                    
                    # Calculate the length of the protein and nucleotide sequences, as gffs format is 1-based, add 1 for the length
                    domain_dict["AA Length"]=int(domain_dict["protein_end"])-int(domain_dict["protein_start"])+1
                    domain_dict["Genomic Coordiantes"]=f"{genomic_loc[1]}-{genomic_loc[2]}"
                    domain_dict["NUC Length"]=int(genomic_loc[2])-int(genomic_loc[1])+1
            
                    output_list.append(domain_dict)
                except (ValueError, KeyError) as e:
                    # Log error if there's an exception
                    logging.error(f"Error on {row['Transcript ID']} {item} {row['Genome Build']}")
                    logging.exception(e)
                    continue
        return pd.DataFrame(output_list)[["Gene Symbol", "Genome_Build", "Chrom", "Domain Name", "Domain Coordinates", "AA Length", "Genomic Coordiantes", "NUC Length"]]

    def map_protein_to_genomic(self, transcript_id: str, protein_start: int, protein_end: int, genome_build: str) -> tuple:
        """
        Maps a protein region to its corresponding genomic coordinates.

        Args:
            transcript_id (str): The transcript ID.
            protein_start (int): The start position of the protein region.
            protein_end (int): The end position of the protein region.
            genome_build (str): The genome build version (e.g., "hg19", "hg38").

        Returns:
            tuple: A tuple containing the chromosome, genomic start position, and genomic end position.

        Raises:
            ValueError: If no transcript is found.
            ValueError: If no CDS (Coding DNA Sequence) is found for the given protein region.

        """
        # Convert genome build to the appropriate format
        if genome_build.lower() == "hg19" or genome_build.lower() == "grch37":
            genome_build = "GRCh37"
        elif genome_build.lower() == "hg38" or genome_build.lower() == "grch38":
            genome_build = "GRCh38"
        else:
            raise ValueError(f"Genome build {genome_build} not supported")
        
        # Filter the transcript dataframe based on transcript ID and genome build
        transcript_df=self.cds_df[self.cds_df["transcript_id"].str.startswith(transcript_id) 
                                  & (self.cds_df["Genome_Build"]== genome_build)].copy()

        if transcript_df.empty:
            raise ValueError("No transcript found")
        elif transcript_df.iloc[0]["Strand"]=="+":
            transcript_df.sort_values("Feature Start",inplace=True)
        elif transcript_df.iloc[0]["Strand"]=="-":
            # sort the CDS records in descending order for the reverse strand
            transcript_df.sort_values("Feature Start",inplace=True,ascending=False)

        # Calculate the length of each feature as gffs format is 1-based, add 1 for the length
        transcript_df["length"]=transcript_df["Feature End"].astype(int)-transcript_df["Feature Start"].astype(int)+1
        # Calculate the start and end positions of each CDS records in the transcript
        transcript_df['CDS_start'] = (transcript_df.groupby('protein_id')['length'].cumsum() - transcript_df['length'])
        transcript_df['CDS_end'] = transcript_df.groupby('protein_id')['length'].cumsum()
        
        # Convert protein start and end positions to CDS start and end positions
        CDS_start, CDS_end = protein_start*3, protein_end*3

        # Find the overlapping CDS regions
        overlap_CDS = transcript_df[(transcript_df["CDS_end"] >= CDS_start) & (transcript_df["CDS_start"] <= CDS_end)]

        if overlap_CDS.empty:
            raise ValueError(f"No CDS found for {transcript_id} {protein_start} {protein_end} {genome_build}")
        elif overlap_CDS.iloc[0]["Strand"] == "-":
            # Calculate the genomic start and end positions for the reverse strand, add offsets for the stop and start codons
            genomic_start = int(overlap_CDS.iloc[-1]["Feature End"]) - (CDS_end - overlap_CDS.iloc[-1]["CDS_start"])+1
            genomic_end = int(overlap_CDS.iloc[0]["Feature End"]) - (CDS_start - overlap_CDS.iloc[0]["CDS_start"])+3
        elif overlap_CDS.iloc[0]["Strand"] == "+":
            # Calculate the genomic start and end positions for the forward strand, add offsets for the stop and start codons
            genomic_start = int(overlap_CDS.iloc[0]["Feature Start"]) + (CDS_start - overlap_CDS.iloc[0]["CDS_start"])-3
            genomic_end = int(overlap_CDS.iloc[-1]["Feature Start"]) + (CDS_end - overlap_CDS.iloc[-1]["CDS_start"])-1
        else:
            raise ValueError(f"Format Error: No strand info found in feature file for {transcript_id}")
        
        chromosome = overlap_CDS.iloc[0]["Sequence ID"]
        
        return (chromosome, genomic_start, genomic_end)

# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 21:43:05 2021

@author: James
"""
import re
import math
import os 

path_to_program_directory=os.path.join((os.path.expanduser("~")),"Documents","Mutation_Scanner")

program_directory_to_reference_genomes=os.path.join(path_to_program_directory,"Input","Reference Genome")

def chromosome_fasta_sequence(chromosome): 
    """
    This program simply converts a fasta file into a string,
    with each character representing a base of the chromosome of interest.
    The header line (starting with '>') is removed.
    """
    file_name=f"chr{chromosome}.fa"
    file_path=os.path.join(program_directory_to_reference_genomes,file_name)
    try:
        with open(file_path,"r") as chromosome_fasta_file:
            chromosome_data = chromosome_fasta_file.read()
            # Remove the header line (if present)
            chromosome_data_formatted = chromosome_data.replace(">chr{}".format(str(chromosome)), "")

            # Remove newline characters                                   
            chromosome_data_formatted=re.sub("\n","",chromosome_data_formatted)
        return chromosome_data_formatted
    except FileNotFoundError:
        print(f"File {file_path} not found.")
        return None
    
#chromosome_fasta_sequence(1, program_file_location)



# It reviews how much of a fragments sequence is 
def fragment_coverage(FragmentNumber,FragmentSize,ChromosomeSequence):

    FragmentSequence=(ChromosomeSequence[((FragmentNumber*FragmentSize)+1):((FragmentNumber*FragmentSize)+FragmentSize)])
    Bases_known=str(((FragmentSize-(FragmentSequence.count("N")+FragmentSequence.count("n")))/FragmentSize)*100)+"%"

    return Bases_known


def all_fragments_chromosome_base_coverage(FragmentSize,ChromosomeSequence):
    NumberOfFragments=math.ceil(len(ChromosomeSequence)/FragmentSize)
    list=[fragment_coverage(NumberOfFragments,FragmentSize,ChromosomeSequence)for NumberOfFragments in range(0,(NumberOfFragments),1)] 
    return list

def all_chromosomes_fragments_base_coverage(FragmentSize):
    Initial_list=[all_fragments_chromosome_base_coverage(FragmentSize,chromosome_fasta_sequence(chromosome)) for chromosome in range(1,23)]
    List_of_fragment_mutated_percentages=[item for sublist in Initial_list for item in sublist]
    
    return List_of_fragment_mutated_percentages


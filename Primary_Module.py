# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 16:16:39 2019

@author: jshar
"""
import os,glob,errno,shutil,math,csv,pandas as pd,time,re,numpy as np,collections,getpass,sys,os.path
from time import sleep
from multiprocessing import Pool
from Tertiary_Module import directory_creation
from decimal import *
from Secondary_Module import all_chromosomes_fragments_base_coverage
import decimal 

#Below is a list of all functions
def chromosome_number_fragment_list(desired_fragment_size, chromosome_number):

    # Open the chromosome file from the reference genome directory and read its data.
    with open(os.path.join(program_directory_to_reference_genomes, f"chr{chromosome_number}.fa")) as chromosome_file:
        
        chromosome_data = chromosome_file.read()

    # Replace all newline characters in the file data to create a continuous sequence string.
    formatted_chromosome_data = re.sub("\n", "", chromosome_data)

    # Determine the length of the chromosome data, accounting for the 'Chr' header.
    # If the chromosome number is 10 or higher, an additional digit is included in the header.
    if chromosome_number < 10:

        chromosome_size = len(formatted_chromosome_data[5:])

    else:

        chromosome_size = len(formatted_chromosome_data[6:])

    # Calculate the number of fragments by rounding up the division of chromosome length by fragment size.
    number_of_fragments = math.ceil(chromosome_size / desired_fragment_size)
    
    # Initialize an empty list to hold fragment descriptors.
    chromosome_fragments = []

    # Generate and append fragment descriptors to the list.
    for fragment_index in range(0, number_of_fragments):

        fragment_descriptor = f"Chromosome {chromosome_number} Fragment {fragment_index}, Base pairs {fragment_index * desired_fragment_size + 1} - {fragment_index * desired_fragment_size + desired_fragment_size}"
        
        chromosome_fragments.append(fragment_descriptor)

    # Return the list of fragment descriptors.
    return chromosome_fragments

# This function generates a list of fragment descriptors for all 22 autosomal human chromosomes.
# It utilizes the chromosome_number_fragment_list function to generate fragment descriptors for each chromosome.
# Finally, it flattens these into a single list that represents fragment descriptors across the entire genome.
def genome_multi_chromosome_fragment_generator(selected_fragment_size):
    
    # Generate a list of lists, where each inner list contains fragment descriptors for a specific chromosome.
    # Looping through chromosomes 1 to 22.
    genome_matrix = [chromosome_number_fragment_list(selected_fragment_size, chromosome) for chromosome in range(1, 23)]

    # Flatten the 2D list to create a 1D list of fragment descriptors for all chromosomes.
    # This is done by iterating through each sublist (chromosome_fragments) in genome_matrix.
    genome_list = [fragment_descriptor for chromosome_fragments in genome_matrix for fragment_descriptor in chromosome_fragments]

    # Return the 1D list of all fragment descriptors across all chromosomes.
    return genome_list

#Flagged: Can we remove the initial filler from the variable titled chromosome_mutation_list
# This function is designed to extract positions of mutations that occur on a specific chromosome.
# It takes in a chromosome number and a list of mutation data as its parameters.
# The function returns a list of integer positions where mutations are found on the given chromosome.
def extract_mutation_positions(chromosome_number, mutation_data):
    
    # Initialize an empty list with a filler value "0-1" to avoid empty list issues.
    # The filler value is not representative of any mutation and will be filtered out later.
    chromosome_mutation_list = ["0-1"]
    
    # Loop through each entry in the provided mutation data.
    for mutation_position in mutation_data:
        
        # Check if the current mutation position starts with the specified chromosome number followed by a colon.
        # If so, this implies the mutation is on the desired chromosome.
        if mutation_position.startswith(str(chromosome_number) + ":"):
            
            # Find the index positions for ":" and "-" in the string to extract the position information.
            start_index = mutation_position.find(":")

            end_index = mutation_position.find("-")
            
            # Append only the extracted position to the chromosome-specific mutation list.
            chromosome_mutation_list.append(mutation_position[start_index+1:end_index])
            
    # The list now contains integer positions where mutations have occurred on the specified chromosome.
    return chromosome_mutation_list

# Creates a Dictionary which has mutation count as one key and chromosome fragment location as the other.
def dictionary_of_fragment_mutation(selected_fragment_size,formatted_mutation_data,chromosome_number,chromosome_sizes):
    
    chromosome_mutations={}

    # Adjust for the off-by-one difference between chromosome numbering and list indexing
    adjusted_chromosome_index = chromosome_number - 1
    
    # Use the adjusted index to fetch the chromosome size
    chromosome_size = chromosome_sizes[adjusted_chromosome_index]

    number_of_fragments=math.ceil(chromosome_size/selected_fragment_size)
    
    for fragment_index in range(0,(number_of_fragments)):

        start_fragment = fragment_index * selected_fragment_size + 1

        end_fragment = fragment_index * selected_fragment_size + selected_fragment_size

        fragment_descriptor=(f"Chromosome {chromosome_number} Fragment {fragment_index}, Basepairs {start_fragment}-{end_fragment}")
        
        chromosome_mutations.update({fragment_descriptor:0})

        for mutation in formatted_mutation_data:

            if mutation=="0-1":
    
                continue
                        
            elif start_fragment <= int(mutation) < end_fragment:
                
                # Update or initialize the count for the given fragment_descriptor
                current_count = chromosome_mutations.get(fragment_descriptor, 0)
    
                chromosome_mutations[fragment_descriptor] = current_count + 1

    return(chromosome_mutations)

#Converts chromosome dictionary into chromosome mutation fragment list for dataframe addition. 
def dictionary_to_list(chromosome_dictionary):

    chromosome_list=[*chromosome_dictionary.values()]

    return(chromosome_list)

#Calculate mutation locations for all chromosomes based on selected fragment size.
def all_chromosome_mutation_locations(selected_fragment_size, original_mutation_data, chromosome_sizes):
    """
    Parameters:
    - selected_fragment_size: Size of each fragment for which mutations are to be determined.
    - original_mutation_data: Original data containing mutation information.
    - chromosome_sizes: Dictionary or list containing sizes of each chromosome.
    
    Returns:
    - A list containing mutation data for all chromosome fragments.
    """
    
    # Initialize an empty list to store mutation data for all chromosome fragments.
    all_chromosome_fragments_mutation_data = []
    
    # Iterate over all 22 chromosomes.
    for chromosome_number in range(1, 23):
        
        # Extract mutation positions for the current chromosome from the original data.
        formatted_mutation_data = extract_mutation_positions(chromosome_number, original_mutation_data)
        
        # Generate a dictionary containing fragment mutation data for the current chromosome.
        AllChromosomeFragmentMutationDictionary = dictionary_of_fragment_mutation(selected_fragment_size, formatted_mutation_data, chromosome_number, chromosome_sizes)
        
        # Convert the dictionary to a list and append it to the main list.
        all_chromosome_fragments_mutation_data += (dictionary_to_list(AllChromosomeFragmentMutationDictionary))
    
    # Return the combined list of mutation data for all chromosome fragments.
    return all_chromosome_fragments_mutation_data

def coding_mutation_data(df, sample_list, sample_dataframe_positions, selected_fragment_size, chromosome_sizes):
    
    samples_mutation_points = []
    
    # Filter the dataframe to get mutation positions
    #Sample dataframe positions is a dictionary containing key value pairs like {index:sample_number_position}
    #mutation_positions = df[df['Mutation Position'].notna()]['Mutation Position'].tolist()
    for sample in sample_list:
        
        sample_rows = [row_index  for row_index, sample_value in sample_dataframe_positions.items() if int(sample_value) == sample]
        #Extracts all of a single samples mutation data

        mutation_genome_locations = df.iloc[sample_rows]['Mutation Position'].tolist()
        
        mutation_genome_locations.sort()

        samples_mutation_points.append(all_chromosome_mutation_locations(selected_fragment_size, mutation_genome_locations, chromosome_sizes))

    return samples_mutation_points

def non_coding_mutation_data(df, sample_list, sample_dataframe_positions, selected_fragment_size, chromosome_sizes):
    # Since the logic is the same as coding_mutation_data, you can just call that function
    return coding_mutation_data(df, sample_list, sample_dataframe_positions, selected_fragment_size, chromosome_sizes)


def exome_mutation_frequency_dataframe_creator(mutation_data, sample_list, selected_fragment_size, autosome_chromosomes_fragment_names, fragment_bases_known, cancer_type):
    # Initialize the original DataFrame
    data = {"GenomeFragments": autosome_chromosomes_fragment_names, "Encoded Base Percentage": fragment_bases_known}
    dataframe = pd.DataFrame(data, index=autosome_chromosomes_fragment_names)
    
    # Create a new DataFrame to hold the new columns
    new_columns = pd.DataFrame({str(sample): mutation_data[i] for i, sample in enumerate(sample_list)}, index=dataframe.index)

    print("concatting")

    # Concatenate the original DataFrame with the new DataFrame containing all the new columns
    dataframe = pd.concat([dataframe, new_columns], axis=1)

    file_name=f"{cancer_type} Individual Samples Exome Mutation Frequency.xlsx"

    file_path=os.path.join(program_directory_to_cancer_genome_mutation_frequency,file_name)

    dataframe.to_excel (os.path.join(file_path), index = None, header=True)

#This is a critical function which is responsible for joining the intron mutation frequency data 
#with the exome mutation frequency data to create a whole genome mutation frequency data 
#it ensures that if a sample is already present then it data is combined 
#with the current mutation frequency score
#If it is not present a new column is added

def cancer_exome_intron_mutation_dataframe_creator(intron_mutation_frequency_data,intron_sample_ids,cancer_type):#Non Coding sample mutation data and coding sample mutation data combined into new excel table.
    
    file_name=f"{cancer_type} Individual Samples Exome Mutation Frequency.xlsx"

    file_path=os.path.join(program_directory_to_cancer_genome_mutation_frequency,file_name)
    
    exome_dataframe=pd.read_excel(file_path, sheet_name= "Sheet1")

    samples_with_exome_mutation_data=(exome_dataframe.columns.values.tolist())

    for intron_sample_index in range(0,len(intron_sample_ids)):

        intron_sample=intron_sample_ids[intron_sample_index]

        if str(intron_sample) in samples_with_exome_mutation_data:
            
            # We begin at position two because the first two columns are not sample data, we are iterating through each relevant sample column.
            for exome_sample_index in range(2,(len(samples_with_exome_mutation_data)-2)):
                    
                    if int(samples_with_exome_mutation_data[exome_sample_index])== int(intron_sample): 

                        for exome_sample_match_row_index in range(0,len(exome_dataframe.index)):

                            updated_exome_sample=int((intron_mutation_frequency_data[intron_sample])[exome_sample_match_row_index])+int(exome_dataframe.iat[(exome_sample_match_row_index),(exome_sample_index+2)])#Problem Code
                            
                            exome_dataframe.at[exome_sample_match_row_index,samples_with_exome_mutation_data[(exome_sample_index+2)]]=int(updated_exome_sample)#Works correctly    
        
        else:

            new_intron_sample=intron_mutation_frequency_data[intron_sample]

            exome_dataframe[str(intron_sample_ids[intron_sample])] = new_intron_sample

    file_name=f"{cancer_type} Individual Samples Whole Genome Mutation Frequency Distribution.xlsx"

    file_path=os.path.join(program_directory_to_cancer_genome_mutation_frequency,file_name)

    exome_dataframe.to_excel (file_path, index = None, header=True)

def fragment_mean(dataframe):

    sums=0

    for i in range(0,len(dataframe)):

        sums=sums+dataframe.iat[i,2]

    newsum=sums/(len(dataframe))

    return(newsum)

def pattern(dataframe,mean):

    pattern_yes_no=[]

    for i in range(0,len(dataframe)):

        value=dataframe.iat[i,2]

        if value <mean:

            pattern_yes_no.append(0)

        else:

            pattern_yes_no.append(1)

    return(pattern_yes_no)

def biomarker_list(dataframe):

    return dataframe[dataframe.iloc[:, 2] == 1].iloc[:, 1].tolist()

def mutation_rate_list(dataframe):

    return dataframe[dataframe.iloc[:, 2] == 1].iloc[:, 0].tolist()

def compute_continuous_mean(dataframe):

    return Decimal(dataframe['MutationRate'].mean())

def compute_continuous_pattern(dataframe, mean_value):

    return dataframe['MutationRate'].apply(lambda x: 0 if x < mean_value else 1)

def compute_continuous_mean(dataframe):

    return Decimal(dataframe['MutationRate'].mean())

def compute_continuous_pattern(dataframe, mean_value):

    return dataframe['MutationRate'].apply(lambda x: 0 if x < mean_value else 1)

def mean_table(cancer_type, iteration):

    file_name=(f"{cancer_type} Biomarkers with mutation frequency data highest 50% iteration {iteration}.xlsx")
    
    file_path=(program_directory_to_cancer_biomarker_candidates,file_name)
    
    dataframe= pd.read_excel(file_path)
    
    mean_value = compute_continuous_mean(dataframe)
    dataframe['Pattern'] = compute_continuous_pattern(dataframe, mean_value)
    
    file_name=(f"{cancer_type} Biomarkers with mutation frequency data highest 50% iteration {iteration+0.5}.xlsx")

    file_path=(program_directory_to_cancer_biomarker_candidates,file_name)

    dataframe.to_excel(file_path, index=None, header=True)
    
    # Assuming 'biomarker_list' and 'mutation_rate_list' are functions you've defined elsewhere
    dataframe['NewBiomarkers'] = biomarker_list(dataframe)

    dataframe['NewMutationRates'] = mutation_rate_list(dataframe)
    
    file_name=(f"{cancer_type} Biomarkers with mutation frequency data highest 50% iteration {iteration+1}.xlsx")

    file_path=(program_directory_to_cancer_biomarker_candidates,file_name)
    
    dataframe[['NewBiomarkers', 'NewMutationRates']].to_excel(file_path, index=None, header=True)

#  This performs the setup of the programs file structure.
directory_creation() 

# This is a list of file paths for different part of the program to save or load data

path_to_program_directory=os.path.join((os.path.expanduser("~")),"Documents","Mutation_Scanner")

program_directory_to_formatted_data=os.path.join(path_to_program_directory,"Input","Samples","Formatted")

program_directory_to_crude_data=os.path.join(path_to_program_directory,"Input","Samples","Crude data")

program_directory_to_reference_genomes=os.path.join(path_to_program_directory,"Input","Reference Genome")

program_directory_to_cancer_genome_mutation_frequency=os.path.join(path_to_program_directory,"Output","Cancer Genomes Mutation Frequency Distribution")

program_directory_to_cancer_biomarker_candidates=os.path.join(path_to_program_directory,"Output","Cancer Biomarker Candidates")

#Important variables used in numerous functions


#selected_fragment_size=int(input("Enter genome fragment size."))

selected_fragment_size=1000000

print(f"{selected_fragment_size}bp selected as fragment size.")

exome_data_mutation_location_column=25

exome_data_sample_column=5

intron_data_mutation_location_column=15

intron_data_sample_column=1

cancer_type="Bone"

#cancer_type=input("Enter Cancer you wish to process.")

print(f"{cancer_type} selected as cancer type.")

#This is a list of chromosome and the number of base pairs they have in there sequences.
chromosome_sizes=[]   

for chromosome_number in range(1,23):

    with open(os.path.join(program_directory_to_reference_genomes,f"chr{chromosome_number}.fa")) as chromosome_file:
            chromosome_data = chromosome_file.read()

            formatted_chromosome_data=re.sub("\n","",chromosome_data)

            if chromosome_number<10:

                chromosome_sizes.append(len(formatted_chromosome_data[5:]))

            else:
                chromosome_sizes.append(len(formatted_chromosome_data[6:]))

print("1.Chromosome Sizes Calculated.")                             

# These are lists,sets and variable used to store cancer patient sample information, cancer patient mutation information, row position,
# and a unique set of pairs of sample and mutation data.
exome_samples=[]

exome_mutation_genome_positions=[]

exome_row_counter=0

exome_mutation_sample_unique_pairs=set()

# Open a cancer mutation file and stores the data in a variable
#The function below would be the same as coding file would be the same as non coding, 
#the reason it not because the dataframe has a different structure.
with open(os.path.join(program_directory_to_crude_data,f"{cancer_type}_Coding.csv")) as exome_cancer_mutation__data:
    #This open the cancer types mutation file as a csv file
    csv_exome_cancer_mutation_data=csv.reader(exome_cancer_mutation__data)
    # This is a for loop which will iterate through all rows present in the csv file
    for row in csv_exome_cancer_mutation_data:

        #Skips the first row, as this contains the columns titles.
        if exome_row_counter==0:
            exome_row_counter+=1
        else:
            #If the mutation data in a row is missing then it is skipped.
            if row[exome_data_mutation_location_column]=="null":
                continue
            else:
                #Selects columns 6 and columns 26 of each row in the table and prevents duplicate mutations in the same sample.
                if (str(row[exome_data_sample_column])+str(row[exome_data_mutation_location_column])) not in exome_mutation_sample_unique_pairs:

                    # Each row represents a specific mutation identified within a particular patient cancer tissue sample
                    # The rows patient sample number is added to a list called exome_samples 
                    exome_samples.append(row[exome_data_sample_column])
                    #The rows mutation is added to a list called exome_mutation_genome_positions
                    exome_mutation_genome_positions.append(row[exome_data_mutation_location_column])
                    #The specific mutation patient pair is added to a set of unique pairs to prevent duplicates being added.
                    exome_mutation_sample_unique_pairs.add(str(row[exome_data_sample_column])+str(row[exome_data_mutation_location_column]))

#This determines the number of times an exome file has to be split for all data to be viewable in excel 
number_of_cancer_exome_mutation_chunks=math.ceil((len(exome_samples))/(selected_fragment_size))

#This process splits up the large exome cancer file into  smaller exome files which can be viewed.
for data_chunk_index in range(0,number_of_cancer_exome_mutation_chunks):

    initial_base_position=(selected_fragment_size*data_chunk_index)

    end_base_position=(selected_fragment_size*(data_chunk_index+1))

    #Creates a dictionary called file_data. 
    exome_file_data= {"sample Id":exome_samples[initial_base_position:end_base_position],"Mutation Position":exome_mutation_genome_positions[initial_base_position:end_base_position]}

    # Dictionary file_data used to make a dataframe
    exome_mutation_chunk_dataframe=pd.DataFrame(exome_file_data,index=(exome_samples[initial_base_position:end_base_position]))
    
    exome_chunks_file_save_name=f"{cancer_type} cancer exome data chunk {data_chunk_index}.xlsx"

    exome_chunk_file_path=os.path.join(program_directory_to_formatted_data,exome_chunks_file_save_name)

    exome_mutation_chunk_dataframe.to_excel(exome_chunk_file_path, index = None, header=True)

intron_samples=[]

intron_mutation_genome_positions=[]

intron_row_counter=0 

intron_mutation_sample_unique_pairs=set()

with open (os.path.join(program_directory_to_crude_data,f"{cancer_type}_Non_Coding.csv")) as intron_cancer_mutation__data:
    
    csv_intron_cancer_mutation_data=csv.reader(intron_cancer_mutation__data)

    for row in csv_intron_cancer_mutation_data:
        if intron_row_counter==0:
            intron_row_counter+=1
        else:
            if row[intron_data_mutation_location_column]=="null":
                continue
            else:
                if (str(row[intron_data_sample_column])+str(row[intron_data_mutation_location_column])) not in intron_mutation_sample_unique_pairs:

                    intron_samples.append(row[intron_data_sample_column])

                    intron_mutation_genome_positions.append(row[intron_data_mutation_location_column])

                    intron_mutation_sample_unique_pairs.add(str(row[intron_data_sample_column])+str(row[intron_data_mutation_location_column]))

number_of_cancer_intron_mutation_chunks=math.ceil((len(intron_samples))/(selected_fragment_size))

for data_chunk_index in range(0,number_of_cancer_intron_mutation_chunks):

    initial_base_position=(selected_fragment_size*data_chunk_index)

    end_base_position=(selected_fragment_size*(data_chunk_index+1))

    intron_file_Data= {"sample Id":intron_samples[initial_base_position:end_base_position],"Mutation Position":intron_mutation_genome_positions[initial_base_position:end_base_position]}

    intron_mutation_chunk_dataframe=pd.DataFrame(intron_file_Data,index=(intron_samples[initial_base_position:end_base_position]))

    intron_chunk_file_save_name=f"{cancer_type} cancer intron data chunk {data_chunk_index}.xlsx"

    intron_chunk_file_path=os.path.join(program_directory_to_formatted_data,intron_chunk_file_save_name)
    
    intron_mutation_chunk_dataframe.to_excel(intron_chunk_file_path, index = None, header=True)

print("2.data processed into accessible excel files.")

#Looks through all files in the directory location listed and counts the number of files

names_of_formated_mutation_files=([x[2] for x in os.walk(program_directory_to_formatted_data)])[0]

number_of_processed_coding_files=0

for file_name in names_of_formated_mutation_files:

    if f"{cancer_type} cancer exome data chunk" in file_name:

        number_of_processed_coding_files+=1

number_of_processed_non_coding_files=0

for file_name in names_of_formated_mutation_files:

    if f"{cancer_type} cancer intron data chunk" in file_name:

        number_of_processed_non_coding_files+=1
        
#Commence Refactoring
for data_chunk_index in range(0,number_of_processed_coding_files):

    file_of_interest=f"{cancer_type} cancer exome data chunk {data_chunk_index}.xlsx"

    file_path=os.path.join(program_directory_to_formatted_data,file_of_interest)

    #This opens the coding mutation data excel file.
    exome_mutation_file_data=pd.read_excel(file_path)
    
    print("3.Coding Excel Input File opened and Cordinates for Column containing ID_SAMPLE discovered.")     
   
    #This is a list containing coding sample duplicates, it does this by adding each row value 
    # underneath the column sample id.
    exome_cancer_cohort_samples_with_repeats = exome_mutation_file_data.iloc[:, 0].tolist()
    
    print("List of Coding samples containing duplicates created.")
    
    # Using a set to get unique samples and then converting back to a list
    exome_cancer_cohort_samples = list(set(int(sample) for sample in exome_cancer_cohort_samples_with_repeats))
    
    print("4.List of Coding samples without duplicates created.")
    
    exome_sample_dataframe_positions = {i: int(sample) for i, sample in enumerate(exome_cancer_cohort_samples_with_repeats)}
    
    print("5.Dictionary of Coding samples with corresponding row coordinates created.")    
    
    all_chromosome_fragments=genome_multi_chromosome_fragment_generator(selected_fragment_size)#This generates a list of all genome fragments.
    
    print("6.Completed genome fragment list generation.")
    
    fragment_sequence_percentage=all_chromosomes_fragments_base_coverage(selected_fragment_size)#This generates a list of all genome fragments base coverage.
    
    print("7.Completed chromosomes fragments base coverage percentage list.")
    
    exome_fragment_mutation_frequency_data=coding_mutation_data(exome_mutation_file_data,exome_cancer_cohort_samples,exome_sample_dataframe_positions,selected_fragment_size,chromosome_sizes)#This is the coding mutation data ready to be added to a dataframe.
    
    print("8.Completed creation of Coding data.")
    
    exome_mutation_frequency_dataframe_creator(exome_fragment_mutation_frequency_data,exome_cancer_cohort_samples,selected_fragment_size,all_chromosome_fragments,fragment_sequence_percentage,cancer_type)#This processes the coding mutation data into a dataframe.
    

# Complete work on a file within a directory which has non coding mutation data for a cohort.

for data_chunk_index in range(0,number_of_processed_coding_files):

    # input file
    file_of_interest=f"{cancer_type} cancer intron data chunk {data_chunk_index}.xlsx"
    
    file_path=os.path.join(program_directory_to_formatted_data,file_of_interest)

    intron_mutation_file_data=pd.read_excel(file_path) #This opens the coding mutation data excel file.
    

    print("10.Non Coding Excel Input file opened, and Coordinates for Column containing ID_Sample discovered.")     
    
    #This is a list containing Non Coding sample duplicates, it does this by adding each row value underneath the column sample id.
    intron_cancer_cohort_samples_with_repeats = intron_mutation_file_data.iloc[:, 0].tolist()
    
    print("11.List of Non Coding samples containing duplicates created.")
    
    # Directly convert the list with potential repeats to a set to get unique values.
    intron_cancer_cohort_samples = list(set(map(int, intron_cancer_cohort_samples_with_repeats)))

    print("12.List of Non Coding samples without duplicates created.")
    
    intron_sample_dataframe_positions = {i: int(sample) for i, sample in enumerate(intron_cancer_cohort_samples_with_repeats)}

    print("13.Dictionary of Non Coding samples with corresponding row coordinates created.")       
    
    intron_fragment_mutation_frequency_data=non_coding_mutation_data(intron_mutation_file_data,intron_cancer_cohort_samples,intron_sample_dataframe_positions,selected_fragment_size,chromosome_sizes)#This is the non coding mutation data ready to be added to a dataframe.
    
    print("14.Completed creation of Non Coding data.")
    
    #This critical stage combines exome and intron data to truly understand mutation frequency hotspots 
    #for certain types of cancer
    cancer_exome_intron_mutation_dataframe_creator(intron_fragment_mutation_frequency_data,intron_cancer_cohort_samples,cancer_type)#This processes the coding and non coding mutation data into a dataframe.


print("15.Completed creation of Coding and Non Coding Excel sample data File.")


print("16.Commencing sample mutation average excel file creation.")

#I need to refactor the code by doing list comprehenesions 

#The mutation frequencies for each individual patient are known, however we are interested in the overall average cancer cohorts
# As a result we need to determine the average mutation frequency for every fragments created so we can see the most common 
# cancer mutation frequency hotspots for each cancer type

#File is read
file_of_interest = f"{cancer_type} Individual Samples Whole Genome Mutation Frequency Distribution.xlsx"

file_path = os.path.join(program_directory_to_cancer_genome_mutation_frequency, file_of_interest)

cancer_samples_mutation_freq_distribution_data = pd.read_excel(file_path)

multi_chromosome_fragment_mutation_freq = []

# Instead of using nrows and ncols, use shape attribute from pandas.
num_rows, num_cols = cancer_samples_mutation_freq_distribution_data.shape

# Note: Adjusted the range to use num_rows and num_cols
for fragment in range(1, num_rows):
    
    fragment_mean_mutation_frequency = 0
    
    # Again, adjusted the range for pandas dataframe
    for cancer_sample in range(2, num_cols):
        
        fragment_mean_mutation_frequency += int(cancer_samples_mutation_freq_distribution_data.iloc[fragment, cancer_sample])

    # Instead of ncols, use num_cols
    multi_chromosome_fragment_mutation_freq.append(fragment_mean_mutation_frequency / (num_cols - 2))

cancer_cohort_average_mutation_frequency_distribution = {
    
    "Chromosome Fragments": all_chromosome_fragments,
    
    "Encoded Base Percentage": fragment_sequence_percentage,
    
    "Samples mean mutation rate per fragment": multi_chromosome_fragment_mutation_freq
}

cancer_cohort_average_mutation_frequency_distribution_dataframe = pd.DataFrame(
    
    cancer_cohort_average_mutation_frequency_distribution,
    
    index=all_chromosome_fragments
)

# File is saved

file_to_save=f"{cancer_type} cohort average mutation frequency distribution.xlsx"

file_path=os.path.join(program_directory_to_cancer_genome_mutation_frequency,file_to_save)

cancer_cohort_average_mutation_frequency_distribution_dataframe.to_excel (file_path, index = None, header=True) 

print("17.Completed Creation of Mean sample Whole Genome Fragment Mutation Rate Excel Table.")

print("18.Biomarker Addition beginning.")

# Importing the genomic data for blood proteins from an Excel file
blood_protein_genomic_file = os.path.join(program_directory_to_crude_data, "Blood Proteins Genomic data.xlsx")

blood_protein_genomic_file_data = pd.read_excel(blood_protein_genomic_file, sheet_name="Sheet1")

# Initialize an empty list to hold proteins ordered by their fragment range
ordered_proteins_by_fragment_range = []

# Loop through each fragment in the average mutation frequency dataframe
for fragment in range(0, len(cancer_cohort_average_mutation_frequency_distribution_dataframe )):
    # Extract chromosome number, start position, and end position from the fragment
    chromosome_number = int((re.search('Chromosome(.*)Fragment', cancer_cohort_average_mutation_frequency_distribution_dataframe .iat[fragment, 0])).group(1))

    start_position = int((re.search('Basepairs(.*)-', cancer_cohort_average_mutation_frequency_distribution_dataframe .iat[fragment, 0])).group(1))

    end_position = int((re.search('-(.*)', cancer_cohort_average_mutation_frequency_distribution_dataframe .iat[fragment, 0])).group(1))
    
    # Initialize an empty list to hold proteins found in the current fragment's genomic location
    proteins_in_fragment_genomic_location = []
    
    # Define the condition for selecting proteins that lie within the current fragment
    condition = (blood_protein_genomic_file_data['Genomic Position'] >= start_position) & \
                (blood_protein_genomic_file_data['Genomic Position'] <= end_position) & \
                (blood_protein_genomic_file_data['Chromosome'] == int(chromosome_number))
    
    # Extract the names of proteins that satisfy the condition
    proteins_in_fragment = blood_protein_genomic_file_data.loc[condition, 'Protein Name'].tolist()
    
    # Append the list of proteins found in the current fragment to the ordered list
    ordered_proteins_by_fragment_range.append(proteins_in_fragment)

# Add the ordered list of proteins as a new column in the average mutation frequency dataframe
cancer_cohort_average_mutation_frequency_distribution_dataframe ["Blood soluble proteins coded in fragment"] = ordered_proteins_by_fragment_range

# Calculate the mean using the fragment_mean function
mean = fragment_mean(cancer_cohort_average_mutation_frequency_distribution_dataframe )

cancer_cohort_average_mutation_frequency_distribution_dataframe ["pattern_yes_no"]=pattern(cancer_cohort_average_mutation_frequency_distribution_dataframe ,mean)

# Creating appropriate filename
excel_filename=f"AvgMutFreq_{cancer_type}_BloodSolProtein.xlsx"

# Full path for saving the Excel file
output_path=os.path.join(program_directory_to_cancer_biomarker_candidates,excel_filename)
                         
# Save the DataFrame to Excel
cancer_cohort_average_mutation_frequency_distribution_dataframe .to_excel(output_path, index = None, header=True)

# Read the previously saved DataFrame containing average mutation frequencies and blood proteins
avg_mut_freq_blood_sol_protein_file = output_path

avg_mut_freq_blood_sol_protein_data = pd.read_excel(avg_mut_freq_blood_sol_protein_file, sheet_name="Sheet1")

# Initialize empty lists to store Blood Soluble Proteins and Mutation Rates
blood_soluble_proteins = []

mutation_rate_data = []

# Loop through the DataFrame to filter rows based on specific conditions
for row in range(len(avg_mut_freq_blood_sol_protein_data)):
    # Check if the column at index 4 is 1 and the list of blood proteins at index 3 is not empty
    if avg_mut_freq_blood_sol_protein_data.iat[row, 4] == 1 and (avg_mut_freq_blood_sol_protein_data.iat[row, 3] != "[]"):

        # Append relevant Blood Soluble Proteins and Mutation Rates to the lists
        blood_soluble_proteins.append(avg_mut_freq_blood_sol_protein_data.iat[row, 3])

        mutation_rate_data.append(avg_mut_freq_blood_sol_protein_data.iat[row, 2])

# Create a simplified DataFrame containing only the rows with Blood Soluble Proteins
blood_proteins_only_fragments_data = {"MutationRate": mutation_rate_data, "Blood soluble proteins coded in fragment": blood_soluble_proteins}

blood_proteins_only_fragments_dataframe = pd.DataFrame(blood_proteins_only_fragments_data, index=mutation_rate_data)

# Save this filtered DataFrame to an Excel file
biomarker_excel_filename = (f"{cancer_type} Biomarkers with mutation frequency data highest 50% iteration 1.xlsx")

biomarker_output_path = os.path.join(program_directory_to_cancer_biomarker_candidates, biomarker_excel_filename)

blood_proteins_only_fragments_dataframe.to_excel(biomarker_output_path, index=None, header=True)

#This part calculates the mean frequency for each subsequenty group which meets the relevant criteria,
# I can definitely streamline this process by selecting the upper 5% means

mean_table(cancer_type,1)

mean_table(cancer_type,2)

mean_table(cancer_type,3)

mean_table(cancer_type,4)

print("19.Program completed.")
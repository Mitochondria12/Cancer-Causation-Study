# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 17:53:28 2023

@author: James
"""
import os
import pandas as pd
import re
import decimal as Decimal

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

path_to_program_directory=os.path.join((os.path.expanduser("~")),"Documents","Mutation_Scanner")

program_directory_to_formatted_data=os.path.join(path_to_program_directory,"Input","Samples","Formatted")

program_directory_to_crude_data=os.path.join(path_to_program_directory,"Input","Samples","Crude data")

program_directory_to_reference_genomes=os.path.join(path_to_program_directory,"Input","Reference Genome")

program_directory_to_cancer_genome_mutation_frequency=os.path.join(path_to_program_directory,"Output","Cancer Genomes Mutation Frequency Distribution")

program_directory_to_cancer_biomarker_candidates=os.path.join(path_to_program_directory,"Output","Cancer Biomarker Candidates")

file_to_load=f"Bone cohort average mutation frequency distribution.xlsx"

file_path=os.path.join(program_directory_to_cancer_genome_mutation_frequency,file_to_load)

cancer_cohort_average_mutation_frequency_distribution_dataframe=pd.read_excel(file_path) 


blood_protein_genomic_file = os.path.join(program_directory_to_crude_data)

blood_protein_genomic_file_data = pd.read_excel(blood_protein_genomic_file)

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
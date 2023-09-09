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

#  This performs the setup of the programs file structure.
directory_creation() 

# This is a list of file paths for different part of the program to save or load data
path_to_program_directory=os.path.join((os.path.expanduser("~")),"Documents","Mutation_Scanner")
program_directory_to_formatted_data=os.path.join(path_to_program_directory,"Input","Samples","Formatted")
program_directory_to_crude_data=os.path.join(path_to_program_directory,"Input","Samples","Crude data")
program_directory_to_reference_genomes=os.path.join(path_to_program_directory,"Input","Reference Genome")
program_directory_to_cancer_genome_mutation_frequency=os.path.join(path_to_program_directory,"Output","Cancer Genomes Mutation Frequency Distribution")
program_directory_to_cancer_biomarker_candidates=os.path.join(path_to_program_directory,"Output","Cancer Biomarker Candidates")
# This function generates a list of fragment descriptors for a specified chromosome.
# The list is used to create a table representing mutation frequency across the chromosome.

def chromosome_number_fragment_list(desired_fragment_size, chromosome_number):

    # Open the chromosome file from the reference genome directory and read its data.
    with open(os.path.join(program_directory_to_reference_genomes, f"Chr{chromosome_number}.fa"), "r") as chromosome_file:
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
    chromosomes={}

    # Adjust for the off-by-one difference between chromosome numbering and list indexing
    adjusted_chromosome_index = chromosome_number - 1
    
    # Use the adjusted index to fetch the chromosome size
    chromosome_size = chromosome_sizes[adjusted_chromosome_index]

    number_of_fragments=math.ceil(chromosome_size/selected_fragment_size)
    
    for fragment_index in range(0,(number_of_fragments)):
        fragment_descriptor=(f"Chromosome {chromosome_number} Fragment {fragment_index}, Basepairs {fragment_index*selected_fragment_size+1}-{fragment_index*selected_fragment_size+selected_fragment_size}")
        chromosomes.update({fragment_descriptor:0})

    for mutation in formatted_mutation_data:
        if mutation==0:
            continue
        else:
            i=math.floor(mutation/1000000)
            fragment_descriptor=(f"Chromosome {chromosome_number} Fragment {fragment_index}, Basepairs {fragment_index*selected_fragment_size+1}-{fragment_index*selected_fragment_size+selected_fragment_size}")
            
            # Update or initialize the count for the given fragment_descriptor
            current_count = chromosomes.get(fragment_descriptor, 0)
            chromosomes[fragment_descriptor] = current_count + 1

    return(chromosomes)

#Converts chromosome dictionary into chromosome mutation fragment list for dataframe addition. 
def dictionary_to_list(chromosome_dictionary):
    chromosome_list=[*chromosome_dictionary.values()]
    return(chromosome_list)
#Determines the number of mutation found within a chromosome fragment location, assigning this mapping to a dictionary,
#Creates a list from the previous dictionary containing number of mutations per specific chromosome mutation fragment list 
def all_chromosome_mutation_locations(selected_fragment_size,original_mutation_data,chromosome_sizes):
    AllChromosomeFragmentMutationList=[]
    for chromosome_number in range(1,23,1):
        formatted_mutation_data=extract_mutation_positions(chromosome_number,original_mutation_data)
        AllChromosomeFragmentMutationDictionary=dictionary_of_fragment_mutation(selected_fragment_size,formatted_mutation_data,chromosome_number,chromosome_sizes) 
        AllChromosomeFragmentMutationList+=(dictionary_to_list(AllChromosomeFragmentMutationDictionary))
    return(AllChromosomeFragmentMutationList)    
   
def coding_mutation_data(excelsheet,sample_list,dictA,selected_fragment_size,chromosome_sizes):#Each sample has its coding mutation data processed and added to a giant list of all sample mutation data.
    Samplesandmutationpoints=[]
    ExcelGenomePositionColumnCoordinates=0
    for i in range(0,excelsheet.ncols):
        if excelsheet.cell_value(0,i)=="MutationPosition":
            ExcelGenomePositionColumnCoordinates+=(int(i))
    
    for i in range(0,len(sample_list)):
        #Begin=time.time()
        SampleData=[SampleTablerowLocation for SampleTablerowLocation, Sample in dictA.items() if int(Sample)==sample_list[i]]
        MutationGenomeLocationList=[excelsheet.cell_value(Samplerows,ExcelGenomePositionColumnCoordinates) for Samplerows in SampleData]
        MutationGenomeLocationList.sort()
        Samplesandmutationpoints.append(all_chromosome_mutation_locations(selected_fragment_size,MutationGenomeLocationList,chromosome_sizes))    
        #Complete=time.time()
        #print(Complete-Begin)
    return(Samplesandmutationpoints)    

def non_coding_mutation_data(excelsheet,sample_list,dictA,selected_fragment_size,chromosome_sizes):#Each sample has its non coding mutation data processed and added to a giant list of all sample mutation data.
    Samplesandmutationpoints=[]
    ExcelGenomePositionColumnCoordinates=0
    for i in range(0,excelsheet.ncols):
        if excelsheet.cell_value(0,i)=="MutationPosition":
            ExcelGenomePositionColumnCoordinates+=(int(i))
    for i in range(0,len(sample_list)):
        SampleData=[SampleTablerowLocation for SampleTablerowLocation, Sample in dictA.items() if int(Sample)==sample_list[i]]
        MutationGenomeLocationList=[excelsheet.cell_value(Samplerows,ExcelGenomePositionColumnCoordinates) for Samplerows in SampleData]
        MutationGenomeLocationList.sort()
        Samplesandmutationpoints.append(all_chromosome_mutation_locations(selected_fragment_size,MutationGenomeLocationList,chromosome_sizes))   
    return(Samplesandmutationpoints)    

def IndividualSamples_ExcelCodingCreator(mutation_data,sample_list,selected_fragment_size,autosome_chromsomes_fragment_names,BasePercentageColumnData,cancer_type,window_directory):# Coding sample mutation data presented into excel table with corresponding chromosome fragment location and base coverage percentages.
    data = {"GenomeFragments":autosome_chromsomes_fragment_names,"Encoded Base Percentage":BasePercentageColumnData}
    df = pd.DataFrame(data,index=autosome_chromsomes_fragment_names)

    for i in range(0,len(sample_list),1):
        value=mutation_data[i]
        df[str(sample_list[i])] = value 
    df.to_excel (r""+str(window_directory)+"/Output/Cancer Genomes Mutation Frequency Distribution/"+cancer_type+"IndividualSamplesWholeGenomeFragmentMutationTable.xlsx", index = None, header=True)
    
def IndividualSamples_ExcelCodingAndNonCodingCreator(NonCodingMutationData,NonCodingSampleNameList,cancer_type,window_directory):#Non Coding sample mutation data and coding sample mutation data combined into new excel table.
    file_name='{name}/Output/Cancer Genomes Mutation Frequency Distribution/{errno}IndividualSamplesWholeGenomeFragmentMutationTable.xlsx'.format(name=window_directory, errno=cancer_type)#(r""+str(window_directory)+("/Output/")+str(cancer_type)+("IndividualSamplesExomeGenomeFragmentMutationTable.xlsx"))
    df=pd.read_excel(io=file_name, sheet_name= "Sheet1")
    SampleListCodingData=(df.columns.values.tolist())#list
    for NonCodingSample in range(0,len(NonCodingSampleNameList)):
        Individual_Sample=NonCodingSampleNameList[NonCodingSample]
        if str(Individual_Sample) in SampleListCodingData:
            for CodingSampleColumn in range(2,(len(SampleListCodingData)-2)):
                    if int(SampleListCodingData[CodingSampleColumn])== int(Individual_Sample): 
                        for CodingSamplerow in range(0,len(df.index)):
                            UpdatedValue=int((NonCodingMutationData[NonCodingSample])[CodingSamplerow])+int(df.iat[(CodingSamplerow),(CodingSampleColumn+2)])#Problem Code
                            df.at[CodingSamplerow,SampleListCodingData[(CodingSampleColumn+2)]]=int(UpdatedValue)#Works correctly    
        else:
            Newvalue=NonCodingMutationData[NonCodingSample]
            df[str(NonCodingSampleNameList[NonCodingSample])] = Newvalue
    df.to_excel (r""+str(window_directory)+"/Output/Cancer Genomes Mutation Frequency Distribution/"+str(cancer_type)+"IndividualSamplesWholeGenomeFragmentMutationTable.xlsx", index = None, header=True)

def fragmentmean(df):
    sums=0
    for i in range(0,len(df)):
        sums=sums+df.iat[i,2]
    newsum=sums/(len(df))
    return(newsum)

def pattern(df,mean):
    pattern_yes_no=[]
    for i in range(0,len(df)):
        value=df.iat[i,2]
        if value <mean:
            pattern_yes_no.append(0)
        else:
            pattern_yes_no.append(1)
    return(pattern_yes_no)
def BiomarkerList(df1):
        A=[]
        for row in range(0,len(df1)):
            if df1.iat[row,2]==1:
                A.append(df1.iat[row,1])
        return(A)

def mutation_rate_list(random_dataframe):
    row_value_one_positions=[]
    # Iterates through each row in a dataframe and if one of the rows has a value of 1 at column 2 then the row value is added to a list [1,1,1,1,1,1,1,1,]
    for row in range(0,len(random_dataframe)):
        if random_dataframe.iat[row,2]==1:
            row_value_one_positions.append(random_dataframe.iat[row,0])
    return(row_value_one_positions)

def Continousfragmentmean(df1):
    sums=0
    for i in range(0,len(df1)):
        sums=sums+Decimal(df1.iat[i,0])
    newsum=sums/(len(df1))
    return(newsum)

def Continouspattern(df1,mean):
    pattern_yes_no=[]
    for i in range(0,len(df1)):
        value=df1.iat[i,0]
        if value <mean:
            pattern_yes_no.append(0)
        else:
            pattern_yes_no.append(1)
    return(pattern_yes_no)



def mean_table(window_directory,cancer_type,iteration):
    file_name4=(r""+str(window_directory)+"/Output/BiomarkerMutationTables/"+str(cancer_type)+"BiomarkerTypeAndBiomarkerLocationMutationTableI"+str(iteration)+".xlsx")
    df1 = pd.read_excel(io=file_name4, sheet_name= "Sheet1")
    meangenomefragmentmutationrate=Continousfragmentmean(df1)
    Biomarkers=[]
    MutationRates=[]
    for row in range(0,len(df1)):
        Biomarkers.append(df1.iat[row,1])
        MutationRates.append(df1.iat[row,0])
    FirstData= {"MutationRate":MutationRates,"Biomarker":Biomarkers,"Pattern":Continouspattern(df1,meangenomefragmentmutationrate)}
    FirstDataExcelTable = pd.DataFrame(FirstData,index=MutationRates)
    FirstDataExcelTable.to_excel(r""+str(window_directory)+"/Output/BiomarkerMutationTables/"+str(cancer_type)+"BiomarkerTypeAndBiomarkerLocationMutationTableI"+str(int(iteration)+0.5)+".xlsx", index = None, header=True)
    file_name=(r""+str(window_directory)+"/Output/BiomarkerMutationTables/"+str(cancer_type)+"BiomarkerTypeAndBiomarkerLocationMutationTableI"+str(int(iteration)+0.5)+".xlsx")
    df2=pd.read_excel(io=file_name, sheet_name= "Sheet1")

    df2 = pd.DataFrame(FirstData,index=MutationRates)
    NewBiomarkers=BiomarkerList(df2)
    NewMutationRates=mutation_rate_list(df2)
    SecondData= {"MutationRate":NewMutationRates,"Biomarker":NewBiomarkers}
    df2 = pd.DataFrame(SecondData,index=NewMutationRates)
    df2.to_excel (r""+str(window_directory)+"/Output/BiomarkerMutationTables/"+str(cancer_type)+"BiomarkerTypeAndBiomarkerLocationMutationTableI"+str(int(iteration)+1)+".xlsx", index = None, header=True)

selected_fragment_size=int(input("Enter genome fragment size."))
print(f"{selected_fragment_size}bp selected as fragment size.")

cancer_type=input("Enter Cancer you wish to process.")
print(f"{cancer_type} selected as cancer type.")


#This is a list of chromosome and the number of base pairs they have in there sequences.
chromosome_sizes=[]   
for chromosome_number in range(1,23):
    with open(os.path.join(program_directory_to_reference_genomes,f"Chr{chromosome_number}.fa","r")) as chromosome_file:
            chromosome_data = chromosome_file.read()
            formatted_chromosome_data=re.sub("\n","",chromosome_data)
            if chromosome_number<10:
                chromosome_sizes.append(len(formatted_chromosome_data[5:]))
            else:
                chromosome_sizes.append(len(formatted_chromosome_data[6:]))

print("1.Chromosome Sizes Calculated.This took ")                             

# These are lists,sets and variable used to store cancer patient sample information, cancer patient mutation information, row position,
# and a unique set of pairs of sample and mutation data.
sample_id_coding=[]
mutation_genome_positions_coding=[]
row_number_coding=0
mutation_sample_unique_pairs_coding=set()

# Open a cancer mutation file and place its data in a variable called CSV FileC
#The function below would be the same as coding file would be the same as non coding, the reason it not because the dataframe has a different structure
with open(os.path.join(program_directory_to_samples,"Crude data",f"{cancer_type}_Coding.csv"),"r") as cancer_mutation_coding_data:
    #This open the cancer types mutation file as a csv file
    csv_cancer_mutation_data=csv.reader(cancer_mutation_coding_data)
    # This is a for loop which will iterate through all rows present in the csv file
    for row in csv_cancer_mutation_data:

        #Skips the first row, as this contains the columns titles.
        if row_number_coding==0:
            row_number_coding+=1
        else:
            #If the mutation data in a row is missing then it is skipped.
            if row[25]=="null":
                continue
            else:
                #Selects columns 6 and columns 26 of each row in the table and prevents duplicate mutations in the same sample.
                if (str(row[5])+str(row[25])) not in mutation_sample_unique_pairs_coding:

                    # Each row represents a specific mutation identified within a particular patient cancer tissue sample
                    # The rows patient sample number is added to a list called sample_id_coding 
                    sample_id_coding.append(row[5])
                    #The rows mutation is added to a list called mutation_genome_positions_coding
                    mutation_genome_positions_coding.append(row[25])
                    #The specific mutation patient pair is added to a set of unique pairs to prevent duplicates being added.
                    mutation_sample_unique_pairs_coding.add(str(row[5])+str(row[25]))

# This variable number_of_files_coding is a number which represents the number of samples divided by one million
number_of_files_coding=((len(sample_id_coding))/(1000000))
#Using the above variable undergo a for loop in range of the number_of_files_coding variable rounded up. 
#This process will split up the large cancer file into subsequent small files which can be processed.
for IterationC in range(0,math.ceil(number_of_files_coding)):
    XC=(1000000*IterationC)
    YC=(1000000*(IterationC+1))
    #Creates a dictionary called file_data. 
    file_data= {"SampleId":sample_id_coding[XC:YC],"MutationPosition":mutation_genome_positions_coding[XC:YC]}
    # Dictionary file_data used to make a dataframe
    DFC=pd.DataFrame(file_data,index=(sample_id_coding[XC:YC]))
    DFC.to_excel((r""+str(window_directory)+"/Input/Samples/SmallFiles/"+str(cancer_type)+"_Coding_Sub"+str(IterationC+1)+".xlsx"), index = None, header=True)


#(str(window_directory)+"/Input/Samples/"+(str(cancer_type)+"_Non_Coding.csv))
SampleIdNC=[]
MutationPositionNC=[]
rowCounterNC=0 
DuplicantFreeNC=set()

with open (r""+str(window_directory)+"/Input/Samples/Crude data/"+(str(cancer_type)+"_Non_Coding.csv")) as InformationNC:
    CSV_FileNC=csv.reader(InformationNC,delimiter=",")
    for row in CSV_FileNC:
        if rowCounterNC==0:
            rowCounterNC+=1
        else:
            if row[15]=="null":
                continue
            else:
                if (str(row[1])+str(row[15])) not in DuplicantFreeNC:
                    SampleIdNC.append(row[1])
                    MutationPositionNC.append(row[15])
                    DuplicantFreeNC.add(str(row[1])+str(row[15]))

NolistsNC=((len(SampleIdNC))/(1000000))

for IterationNC in range(0,math.ceil(NolistsNC)):#round(Nolists)+1)):
    XNC=(1000000*IterationNC)
    YNC=(1000000*(IterationNC+1))
    FileDataNC= {"SampleId":SampleIdNC[XNC:YNC],"MutationPosition":MutationPositionNC[XNC:YNC]}
    DFNC=pd.DataFrame(FileDataNC,index=(SampleIdNC[XNC:YNC]))
    DFNC.to_excel((r""+str(window_directory)+"/Input/Samples/SmallFiles/"+(str(cancer_type)+"_Non_Coding_Sub"+str(IterationNC+1)+".xlsx")), index = None, header=True)
print("2.data processed into accessible excel files.This took "+str(time.time()-StartTime)+" Seconds.")

#Looks through all files in the directory location listed and counts the number of files
names_of_formated_mutation_files=([x[2] for x in os.walk(r""+str(window_directory)+"/Input/Samples/SmallFiles/")])[0]
number_of_processed_coding_files=0
for file_name in names_of_formated_mutation_files:
    if str(cancer_type+"_Coding_Sub") in file_name:
        number_of_processed_coding_files+=1

number_of_processed_non_coding_files=0
for file_name in names_of_formated_mutation_files:
    if str(cancer_type+"_Non_Coding_Sub") in file_name:
        number_of_processed_non_coding_files+=1

        
totalsamples=len(SampleIdNC)+len(sample_id_coding)
if(totalsamples)>16380:
    print("Sheet Size not large enough.")
else:
    print("Sheet size large enough.")
for i in range(0,number_of_processed_coding_files):
    if i ==0:
        CancerCodingSampleData=((str(window_directory)+"/Input/Samples/SmallFiles/"+(str(cancer_type)+"_Coding_Sub1.xlsx")))
        CWB=xlrd.open_workbook(CancerCodingSampleData) #This opens the coding mutation data excel file.
        CodingExcelSheet=CWB.sheet_by_index(0) 
        CodingExcelSheet.cell_value(0, 0) 
        
        print("3.Coding Excel Input File opened and Cordinates for Column containing ID_SAMPLE discovered.This took "+str(time.time()-StartTime)+" Seconds.")     
        SamplesWithDuplicates=[CodingExcelSheet.cell_value(i,0) for i in range(1,CodingExcelSheet.nrows)] #This is a list containing coding sample duplicates, it does this by adding each row value underneath the column sample id.
        print("List of Coding samples containing duplicates created.This took "+str(time.time()-StartTime)+" Seconds.")
        SampleListC=[]#Might be able to use a list comprehension to shorten this code.
        for Sample in SamplesWithDuplicates:
            if int(Sample) not in SampleListC:
                SampleListC.append(int(Sample))
        
        print("4.List of Coding samples without duplicates created.This took "+str(time.time()-StartTime)+" Seconds.")
        DictCSamplerows={}#CodingDictionary ,This is a dictionary of Coding Sample with there row excel locations.
        CounterC=1    
        for Sample in SamplesWithDuplicates:
            DictCSamplerows.update({CounterC:int(Sample)})
            CounterC=CounterC+1
        print("5.Dictionary of Coding samples with corresponding row coordinates created.This took "+str(time.time()-StartTime)+" Seconds.")    
        
        all_chromosome_fragments=genome_multi_chromosome_fragment_generator(selected_fragment_size)#This generates a list of all genome fragments.
        print("6.Completed genome fragment list generation.This took"+str(time.time()-StartTime)+"Seconds.")
        fragment_sequence_percentage=all_chromosomes_fragments_base_coverage(window_directory,selected_fragment_size)#This generates a list of all genome fragments base coverage.
        print("7.Completed chromosomes fragments base coverage percentage list.This took "+str(time.time()-StartTime)+" Seconds.")
        
        CodingData=coding_mutation_data(CodingExcelSheet,SampleListC,DictCSamplerows,selected_fragment_size,chromosome_sizes)#This is the coding mutation data ready to be added to a dataframe.
        print("8.Completed creation of Coding data.This took "+str(time.time()-StartTime)+" Seconds.")
        IndividualSamples_ExcelCodingCreator(CodingData,SampleListC,selected_fragment_size,all_chromosome_fragments,fragment_sequence_percentage,cancer_type,window_directory)#This processes the coding mutation data into a dataframe.
    else:
        #input Samples small files
        CancerNonCodingSampleData=((str(window_directory)+"/Input/Samples/SmallFiles/"+(str(cancer_type)+"_Coding_Sub"+str(i+1)+".xlsx")))
        NCWB=xlrd.open_workbook(CancerNonCodingSampleData) #This opens the coding mutation data excel file.
        NonCodingExcelSheet=NCWB.sheet_by_index(0) 
        NonCodingExcelSheet.cell_value(0, 0) 
        ExcelSampleIdColumnNC=0#This value equals the Non Coding excel table location of the column ID sample, the below code is how it is discovered.
        for i in range(0,NonCodingExcelSheet.ncols):
            if NonCodingExcelSheet.cell_value(0,i)=="SampleId":
                ExcelSampleIdColumnNC+=(int(i))
        print("10.Non Coding Excel Input file opened, and Coordinates for Column containing ID_Sample discovered.This took "+str(time.time()-StartTime)+" Seconds.")     
        SamplesWithDuplicatesNC=[NonCodingExcelSheet.cell_value(i,ExcelSampleIdColumnNC) for i in range(1,NonCodingExcelSheet.nrows)]#This is a list containing Non Coding sample duplicates, it does this by adding each row value underneath the column sample id.
        print("11.List of Non Coding samples containing duplicates created.This took "+str(time.time()-StartTime)+" Seconds.")
        SampleListNCD=set()#Might be able to use a list comprehension to shorten this code.
        for Sample in SamplesWithDuplicatesNC:
            if Sample not in SampleListNCD:
                SampleListNCD.add(int(Sample))
        SampleListNC=[]
        for a in SampleListNCD:
            SampleListNC.append(int(a))
        print("12.List of Non Coding samples without duplicates created.This took "+str(time.time()-StartTime)+" Seconds.")
        DictNCSamplerows={}#NonCoding Dictionary ,This is a dictionary of Non Coding Sample with there row excel locations.
        CounterNC=1    
        for Sample in SamplesWithDuplicatesNC:
            DictNCSamplerows.update({CounterNC:int(Sample)})
            CounterNC=CounterNC+1
        print("13.Dictionary of Non Coding samples with corresponding row coordinates created.This took "+str(time.time()-StartTime)+" Seconds.")       
        NonCodingData=non_coding_mutation_data(NonCodingExcelSheet,SampleListNC,DictNCSamplerows,selected_fragment_size,chromosome_sizes)#This is the non coding mutation data ready to be added to a dataframe.
        print("14.Completed creation of Non Coding data.This took "+str(time.time()-StartTime)+" Seconds.")
        IndividualSamples_ExcelCodingAndNonCodingCreator(NonCodingData,SampleListNC,cancer_type,window_directory)
    print("Processed coding files")
for i in range(0,number_of_processed_coding_files):

    # input file
    CancerNonCodingSampleData=((str(window_directory)+"/Input/Samples/SmallFiles/"+(str(cancer_type)+"_Non_Coding_Sub"+str(i+1)+".xlsx")))
    NCWB=xlrd.open_workbook(CancerNonCodingSampleData) #This opens the coding mutation data excel file.
    NonCodingExcelSheet=NCWB.sheet_by_index(0) 
    NonCodingExcelSheet.cell_value(0, 0) 
    ExcelSampleIdColumnNC=0#This value equals the Non Coding excel table location of the column ID sample, the below code is how it is discovered.
    for i in range(0,NonCodingExcelSheet.ncols):
        if NonCodingExcelSheet.cell_value(0,i)=="SampleId":
            ExcelSampleIdColumnNC+=(int(i))
    print("10.Non Coding Excel Input file opened, and Coordinates for Column containing ID_Sample discovered.This took "+str(time.time()-StartTime)+" Seconds.")     
    SamplesWithDuplicatesNC=[NonCodingExcelSheet.cell_value(i,ExcelSampleIdColumnNC) for i in range(1,NonCodingExcelSheet.nrows)]#This is a list containing Non Coding sample duplicates, it does this by adding each row value underneath the column sample id.
    print("11.List of Non Coding samples containing duplicates created.This took "+str(time.time()-StartTime)+" Seconds.")
    SampleListNCD=set()#Might be able to use a list comprehension to shorten this code.
    for Sample in SamplesWithDuplicatesNC:
        if Sample not in SampleListNCD:
            SampleListNCD.add(int(Sample))
    SampleListNC=[]
    for a in SampleListNCD:
        SampleListNC.append(int(a))
    print("12.List of Non Coding samples without duplicates created.This took "+str(time.time()-StartTime)+" Seconds.")
    DictNCSamplerows={}#NonCoding Dictionary ,This is a dictionary of Non Coding Sample with there row excel locations.
    CounterNC=1    
    for Sample in SamplesWithDuplicatesNC:
        DictNCSamplerows.update({CounterNC:int(Sample)})
        CounterNC=CounterNC+1
    print("13.Dictionary of Non Coding samples with corresponding row coordinates created.This took "+str(time.time()-StartTime)+" Seconds.")       
    NonCodingData=non_coding_mutation_data(NonCodingExcelSheet,SampleListNC,DictNCSamplerows,selected_fragment_size,chromosome_sizes)#This is the non coding mutation data ready to be added to a dataframe.
    print("14.Completed creation of Non Coding data.This took "+str(time.time()-StartTime)+" Seconds.")
    IndividualSamples_ExcelCodingAndNonCodingCreator(NonCodingData,SampleListNC,cancer_type,window_directory)#This processes the coding and non coding mutation data into a dataframe.
print("15.Completed creation of Coding and Non Coding Excel Sample data File.This took "+str(time.time()-StartTime)+" Seconds.")






print("16.Commencing sample mutation average excel file creation.")

#I need to refactor the code by doing list comprehenesions 

#The mutation frequencies for each individual patient are known, however we are interested in the overall average cancer cohorts
# As a result we need to determine the average mutation frequency for every fragments created so we can see the most common 
# cancer mutation frequency hotspots for each cancer type

#File is read
file_of_interest=f"{cancer_type}IndividualSamplesWholeGenomeFragmentMutationTable.xlsx"

file_path=os.path.join(program_directory_to_cancer_genome_mutation_frequency,file_of_interest)

cancer_samples_mutation_freq_distribution_file=xlrd.open_workbook(file_path) 

cancer_samples_mutation_freq_distribution_data=cancer_samples_mutation_freq_distribution_file.sheet_by_index(0) 

multi_chromosome_fragment_mutation_freq=[]

for fragment in range(1,cancer_samples_mutation_freq_distribution_data.nrows):

    fragment_mean_mutation_frequency=0

    for cancer_sample in range(2,cancer_samples_mutation_freq_distribution_data.ncols):

        fragment_mean_mutation_frequency=fragment_mean_mutation_frequency+int(cancer_samples_mutation_freq_distribution_data.cell_value(fragment,cancer_sample))
    
    multi_chromosome_fragment_mutation_freq.append(fragment_mean_mutation_frequency/(cancer_samples_mutation_freq_distribution_data.ncols-2))

cancer_cohort_average_mutation_frequency_distribution = {"Chromosome Fragments":all_chromosome_fragments,"Encoded Base Percentage":fragment_sequence_percentage,"Samples mean mutatione rate per fragment":multi_chromosome_fragment_mutation_freq}

cancer_cohort_average_mutation_frequency_distribution_dataframe = pd.DataFrame(cancer_cohort_average_mutation_frequency_distribution,index=all_chromosome_fragments)

# File is saved

file_to_save=f"{cancer_type} cohort average mutation frequency distribution.xlsx"

file_path=os.path.join(program_directory_to_cancer_genome_mutation_frequency,file_to_save)

cancer_cohort_average_mutation_frequency_distribution_dataframe.to_excel (file_path, index = None, header=True) 

print("17.Completed Creation of Mean Sample Whole Genome Fragment Mutation Rate Excel Table.This took "+str(time.time()-StartTime)+" Seconds.")

print("18.Biomarker Addition beginning.")

# Importing the genomic data for blood proteins from an Excel file
blood_protein_genomic_file = os.path.join(program_directory_to_crude_data, "Blood Proteins Genomic data.xlsx")

blood_protein_genomic_file_data = pd.read_excel(blood_protein_genomic_file, sheet_name="Sheet1")

# Initialize an empty list to hold proteins ordered by their fragment range
ordered_proteins_by_fragment_range = []

# Loop through each fragment in the average mutation frequency dataframe
for fragment in range(0, len(average_mutation_frequency_dataframe)):
    # Extract chromosome number, start position, and end position from the fragment
    chromosome_number = int((re.search('Chromosome(.*)Fragment', average_mutation_frequency_dataframe.iat[fragment, 0])).group(1))

    start_position = int((re.search('Basepairs(.*)-', average_mutation_frequency_dataframe.iat[fragment, 0])).group(1))

    end_position = int((re.search('-(.*)', average_mutation_frequency_dataframe.iat[fragment, 0])).group(1))
    
    # Initialize an empty list to hold proteins found in the current fragment's genomic location
    proteins_in_fragment_genomic_location = []
    
    # Define the condition for selecting proteins that lie within the current fragment
    condition = (blood_protein_genomic_data['Genomic Position'] >= start_position) & \
                (blood_protein_genomic_data['Genomic Position'] <= end_position) & \
                (blood_protein_genomic_data['Chromosome'] == int(chromosome_number))
    
    # Extract the names of proteins that satisfy the condition
    proteins_in_fragment = blood_protein_genomic_data.loc[condition, 'Protein Name'].tolist()
    
    # Append the list of proteins found in the current fragment to the ordered list
    ordered_proteins_by_fragment_range.append(proteins_in_fragment)

# Add the ordered list of proteins as a new column in the average mutation frequency dataframe
average_mutation_frequency_dataframe["Blood soluble proteins coded in fragment"] = ordered_proteins_by_fragment_range

# Calculate the mean using the fragmentmean function
mean = fragmentmean(average_mutation_frequency_dataframe)

average_mutation_frequency_dataframe["pattern_yes_no"]=pattern(average_mutation_frequency_dataframe,mean)

# Creating appropriate filename
excel_filename=f"AvgMutFreq_{cancer_type}_BloodSolProtein.xlsx"

# Full path for saving the Excel file
output_path=os.path.join(program_directory_to_cancer_biomarker_candidates,excel_filename)
                         
# Save the DataFrame to Excel
average_mutation_frequency_dataframe.to_excel(output_path, index = None, header=True)

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
biomarker_excel_filename = f"{cancer_type}_Biomarker_Mutation_Table_1.xlsx"

biomarker_output_path = os.path.join(program_directory_to_cancer_biomarker_candidates, biomarker_excel_filename)

blood_proteins_only_fragments_dataframe.to_excel(biomarker_output_path, index=None, header=True)

#This part calculates the mean frequency for each subsequenty group which meets the relevant criteria,
# I can definitely streamline this process by selecting the upper 5% means

mean_table(window_directory,cancer_type,1)

mean_table(window_directory,cancer_type,2)

mean_table(window_directory,cancer_type,3)

mean_table(window_directory,cancer_type,4)

print("19.Program completed.")
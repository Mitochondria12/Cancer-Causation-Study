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
program_directory_to_crude_data=os.path.join(path_to_program_directory,"Input","Samples","Crude Data")
program_directory_to_reference_genomes=os.path.join(path_to_program_directory,"Input","Reference Genome")
program_directory_to_output=os.path.join(path_to_program_directory,"Output")

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




def Dictionary_Of_Fragment_mutation(selected_fragment_size,formatted_mutation_data,chromosome_number,chromosome_sizes):# Creates a Dictionary which has mutation count as one key and chromosome fragment location as the other.
    chromosomes={}

    # Adjust for the off-by-one difference between chromosome numbering and list indexing
    adjusted_chromosome_index = chromosome_number - 1
    
    # Use the adjusted index to fetch the chromosome size
    chromosome_size = chromosome_sizes[adjusted_chromosome_index]

    number_of_fragments=math.ceil(chromosome_size/selected_fragment_size)
    
    for fragment_index in range(0,(number_of_fragments)):
        a=(f"Chromosome {chromosome_number} Fragment {fragment_index}, Basepairs {fragment_index*selected_fragment_size+1}-{fragment_index*selected_fragment_size+selected_fragment_size}")
        chromosomes.update({a:0})

    for mutation in formatted_mutation_data:
        if mutation==0:
            continue
        else:
            i=math.floor(mutation/1000000)
            a=("Chromosome"+str(chromosome_number)+"Fragment"+str(i)+",Basepairs "+str(i*selected_fragment_size+1)+"-"+str(i*selected_fragment_size+selected_fragment_size))
            if chromosomes.get(a) ==0:
                    chromosomes.update({a:1})
            else:
                b=((chromosomes.get(a))+1)
                chromosomes.update({a:b})
    return(chromosomes)

def dictionary_to_list(ABC):#Converts chromosome dictionary into chromosome mutation fragment list for dataframe addition. 
    C=[*ABC.values()]
    return(C)

def all_chromosome_mutation_locations(selected_fragment_size,original_mutation_data,chromosome_sizes):#Processes the mutation data for every chromosome,plots each chromosome data mutation point to a chromosome fragment location in  a dictionary,Converts chromosome dictionary into chromosome mutation fragment list for dataframe addition,Adds each chromosome mutation fragment dictionary list conversion into one giant list file.. 
    AllChromosomeFragmentMutationList=[]
    for chromosome_number in range(1,23,1):
        formatted_mutation_data=extract_mutation_positions(chromosome_number,original_mutation_data)
        AllChromosomeFragmentMutationDictionary=Dictionary_Of_Fragment_mutation(selected_fragment_size,formatted_mutation_data,chromosome_number,chromosome_sizes) 
        AllChromosomeFragmentMutationList+=(dictionary_to_list(AllChromosomeFragmentMutationDictionary))
    return(AllChromosomeFragmentMutationList)    
   
def MutationPointDataCoding(excelsheet,SampleList,dictA,selected_fragment_size,chromosome_sizes):#Each sample has its coding mutation data processed and added to a giant list of all sample mutation data.
    Samplesandmutationpoints=[]
    ExcelGenomePositionColumnCoordinates=0
    for i in range(0,excelsheet.ncols):
        if excelsheet.cell_value(0,i)=="MutationPosition":
            ExcelGenomePositionColumnCoordinates+=(int(i))
    
    for i in range(0,len(SampleList)):
        #Begin=time.time()
        SampleData=[SampleTablerowLocation for SampleTablerowLocation, Sample in dictA.items() if int(Sample)==SampleList[i]]
        MutationGenomeLocationList=[excelsheet.cell_value(Samplerows,ExcelGenomePositionColumnCoordinates) for Samplerows in SampleData]
        MutationGenomeLocationList.sort()
        Samplesandmutationpoints.append(all_chromosome_mutation_locations(selected_fragment_size,MutationGenomeLocationList,chromosome_sizes))    
        #Complete=time.time()
        #print(Complete-Begin)
    return(Samplesandmutationpoints)    

def MutationPointDataNonCoding(excelsheet,SampleList,dictA,selected_fragment_size,chromosome_sizes):#Each sample has its non coding mutation data processed and added to a giant list of all sample mutation data.
    Samplesandmutationpoints=[]
    ExcelGenomePositionColumnCoordinates=0
    for i in range(0,excelsheet.ncols):
        if excelsheet.cell_value(0,i)=="MutationPosition":
            ExcelGenomePositionColumnCoordinates+=(int(i))
    for i in range(0,len(SampleList)):
        SampleData=[SampleTablerowLocation for SampleTablerowLocation, Sample in dictA.items() if int(Sample)==SampleList[i]]
        MutationGenomeLocationList=[excelsheet.cell_value(Samplerows,ExcelGenomePositionColumnCoordinates) for Samplerows in SampleData]
        MutationGenomeLocationList.sort()
        Samplesandmutationpoints.append(all_chromosome_mutation_locations(selected_fragment_size,MutationGenomeLocationList,chromosome_sizes))   
    return(Samplesandmutationpoints)    

def IndividualSamples_ExcelCodingCreator(mutationdata,samplelist,selected_fragment_size,GenomeFragmentColumnData,BasePercentageColumnData,cancer_type,Window_Directory):# Coding sample mutation data presented into excel table with corresponding chromosome fragment location and base coverage percentages.
    data = {"GenomeFragments":GenomeFragmentColumnData,"Encoded Base Percentage":BasePercentageColumnData}
    df = pd.DataFrame(data,index=GenomeFragmentColumnData)

    for i in range(0,len(samplelist),1):
        value=mutationdata[i]
        df[str(samplelist[i])] = value 
    df.to_excel (r""+str(Window_Directory)+"/Output/WholeGenomeMutationTables/"+cancer_type+"IndividualSamplesWholeGenomeFragmentMutationTable.xlsx", index = None, header=True)
    
def IndividualSamples_ExcelCodingAndNonCodingCreator(NonCodingMutationData,NonCodingSampleNameList,cancer_type,Window_Directory):#Non Coding sample mutation data and coding sample mutation data combined into new excel table.
    file_name='{name}/Output/WholeGenomeMutationTables/{errno}IndividualSamplesWholeGenomeFragmentMutationTable.xlsx'.format(name=Window_Directory, errno=cancer_type)#(r""+str(Window_Directory)+("/Output/")+str(cancer_type)+("IndividualSamplesExomeGenomeFragmentMutationTable.xlsx"))
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
    df.to_excel (r""+str(Window_Directory)+"/Output/WholeGenomeMutationTables/"+str(cancer_type)+"IndividualSamplesWholeGenomeFragmentMutationTable.xlsx", index = None, header=True)

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

def meantable(Window_Directory,cancer_type,i):
    file_name4=(r""+str(Window_Directory)+"/Output/BiomarkerMutationTables/"+str(cancer_type)+"BiomarkerTypeAndBiomarkerLocationMutationTableI"+str(i)+".xlsx")
    df1 = pd.read_excel(io=file_name4, sheet_name= "Sheet1")
    def BiomarkerList(df1):
        A=[]
        for row in range(0,len(df1)):
            if df1.iat[row,2]==1:
                A.append(df1.iat[row,1])
        return(A)
    def MutationRateList(df1):
        B=[]
        for row in range(0,len(df1)):
            if df1.iat[row,2]==1:
                B.append(df1.iat[row,0])
        return(B)
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
    meangenomefragmentmutationrate=Continousfragmentmean(df1)
    Biomarkers=[]
    MutationRates=[]
    for row in range(0,len(df1)):
        Biomarkers.append(df1.iat[row,1])
        MutationRates.append(df1.iat[row,0])
    FirstData= {"MutationRate":MutationRates,"Biomarker":Biomarkers,"Pattern":Continouspattern(df1,meangenomefragmentmutationrate)}
    FirstDataExcelTable = pd.DataFrame(FirstData,index=MutationRates)
    FirstDataExcelTable.to_excel(r""+str(Window_Directory)+"/Output/BiomarkerMutationTables/"+str(cancer_type)+"BiomarkerTypeAndBiomarkerLocationMutationTableI"+str(int(i)+0.5)+".xlsx", index = None, header=True)
    file_name=(r""+str(Window_Directory)+"/Output/BiomarkerMutationTables/"+str(cancer_type)+"BiomarkerTypeAndBiomarkerLocationMutationTableI"+str(int(i)+0.5)+".xlsx")
    df2=pd.read_excel(io=file_name, sheet_name= "Sheet1")

    df2 = pd.DataFrame(FirstData,index=MutationRates)
    NewBiomarkers=BiomarkerList(df2)
    NewMutationRates=MutationRateList(df2)
    SecondData= {"MutationRate":NewMutationRates,"Biomarker":NewBiomarkers}
    df2 = pd.DataFrame(SecondData,index=NewMutationRates)
    df2.to_excel (r""+str(Window_Directory)+"/Output/BiomarkerMutationTables/"+str(cancer_type)+"BiomarkerTypeAndBiomarkerLocationMutationTableI"+str(int(i)+1)+".xlsx", index = None, header=True)

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
with open(os.path.join(program_directory_to_samples,"Crude Data",f"{cancer_type}_Coding.csv"),"r") as cancer_mutation_coding_data:
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
    DFC.to_excel((r""+str(Window_Directory)+"/Input/Samples/SmallFiles/"+str(cancer_type)+"_Coding_Sub"+str(IterationC+1)+".xlsx"), index = None, header=True)


#(str(Window_Directory)+"/Input/Samples/"+(str(cancer_type)+"_Non_Coding.csv))
SampleIdNC=[]
MutationPositionNC=[]
rowCounterNC=0 
DuplicantFreeNC=set()

with open (r""+str(Window_Directory)+"/Input/Samples/Crude Data/"+(str(cancer_type)+"_Non_Coding.csv")) as InformationNC:
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
    DFNC.to_excel((r""+str(Window_Directory)+"/Input/Samples/SmallFiles/"+(str(cancer_type)+"_Non_Coding_Sub"+str(IterationNC+1)+".xlsx")), index = None, header=True)
print("2.Data processed into accessible excel files.This took "+str(time.time()-StartTime)+" Seconds.")

#Looks through all files in the directory location listed and counts the number of files
names_of_formated_mutation_files=([x[2] for x in os.walk(r""+str(Window_Directory)+"/Input/Samples/SmallFiles/")])[0]
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
        CancerCodingSampleData=((str(Window_Directory)+"/Input/Samples/SmallFiles/"+(str(cancer_type)+"_Coding_Sub1.xlsx")))
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
        
        GenomeFragmentColumn=genome_multi_chromosome_fragment_generator(selected_fragment_size)#This generates a list of all genome fragments.
        print("6.Completed genome fragment list generation.This took"+str(time.time()-StartTime)+"Seconds.")
        BaseCoveragePercentageColumn=all_chromosomes_fragments_base_coverage(Window_Directory,selected_fragment_size)#This generates a list of all genome fragments base coverage.
        print("7.Completed chromosomes fragments base coverage percentage list.This took "+str(time.time()-StartTime)+" Seconds.")
        
        CodingData=MutationPointDataCoding(CodingExcelSheet,SampleListC,DictCSamplerows,selected_fragment_size,chromosome_sizes)#This is the coding mutation data ready to be added to a dataframe.
        print("8.Completed creation of Coding Data.This took "+str(time.time()-StartTime)+" Seconds.")
        IndividualSamples_ExcelCodingCreator(CodingData,SampleListC,selected_fragment_size,GenomeFragmentColumn,BaseCoveragePercentageColumn,cancer_type,Window_Directory)#This processes the coding mutation data into a dataframe.
    else:
        #input Samples small files
        CancerNonCodingSampleData=((str(Window_Directory)+"/Input/Samples/SmallFiles/"+(str(cancer_type)+"_Coding_Sub"+str(i+1)+".xlsx")))
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
        NonCodingData=MutationPointDataNonCoding(NonCodingExcelSheet,SampleListNC,DictNCSamplerows,selected_fragment_size,chromosome_sizes)#This is the non coding mutation data ready to be added to a dataframe.
        print("14.Completed creation of Non Coding Data.This took "+str(time.time()-StartTime)+" Seconds.")
        IndividualSamples_ExcelCodingAndNonCodingCreator(NonCodingData,SampleListNC,cancer_type,Window_Directory)
    print("Processed coding files")
for i in range(0,number_of_processed_coding_files):

    # input file
    CancerNonCodingSampleData=((str(Window_Directory)+"/Input/Samples/SmallFiles/"+(str(cancer_type)+"_Non_Coding_Sub"+str(i+1)+".xlsx")))
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
    NonCodingData=MutationPointDataNonCoding(NonCodingExcelSheet,SampleListNC,DictNCSamplerows,selected_fragment_size,chromosome_sizes)#This is the non coding mutation data ready to be added to a dataframe.
    print("14.Completed creation of Non Coding Data.This took "+str(time.time()-StartTime)+" Seconds.")
    IndividualSamples_ExcelCodingAndNonCodingCreator(NonCodingData,SampleListNC,cancer_type,Window_Directory)#This processes the coding and non coding mutation data into a dataframe.
print("15.Completed creation of Coding and Non Coding Excel Sample Data File.This took "+str(time.time()-StartTime)+" Seconds.")
print("16.Comencing sample mutation average excel file creation.")
#This commences processing a new dataframe containing average mutation rate for all samples per genome fragment.

#Output file
WGWB = xlrd.open_workbook(str(Window_Directory)+"/Output/WholeGenomeMutationTables/"+str(cancer_type)+"IndividualSamplesWholeGenomeFragmentMutationTable.xlsx") 
WholeGenomeExcelSheet=WGWB.sheet_by_index(0) 
WholeGenomeExcelSheet.cell_value(0, 0)
WholeGenomeFragmentMeanMutationList=[]
for Number in range(1,WholeGenomeExcelSheet.nrows):
    Mean=0
    for i in range(2,WholeGenomeExcelSheet.ncols):
        Mean=Mean+int(WholeGenomeExcelSheet.cell_value(Number,i))
    WholeGenomeFragmentMeanMutationList.append(Mean/(WholeGenomeExcelSheet.ncols-2))
print(WholeGenomeExcelSheet.ncols)
Data = {"GenomeFragments":GenomeFragmentColumn,"Encoded Base Percentage":BaseCoveragePercentageColumn,"Samples mean mutatione rate per fragment":WholeGenomeFragmentMeanMutationList}
MeanGenomeExcelTable = pd.DataFrame(Data,index=GenomeFragmentColumn)

# Output File
export_xlsx =MeanGenomeExcelTable.to_excel (r""+str(Window_Directory)+"/Output/WholeGenomeMutationTables/"+str(cancer_type)+"MeanSamplesWholeGenomeFragmentMutationTable.xlsx", index = None, header=True) 
print("17.Completed Creation of Mean Sample Whole Genome Fragment Mutation Rate Excel Table.This took "+str(time.time()-StartTime)+" Seconds.")
print("18.Biomarker Addition beginning.")

Data=MeanGenomeExcelTable.columns
# Input File
file_name2=(r""+str(Window_Directory)+"/Input/Samples/Biomarkers/Genes.xlsx")
PotentialBiomarkersExcelTable=pd.read_excel(io=file_name2, sheet_name="Sheet1")
Biomarker=[]
for fragment in range(0,len(MeanGenomeExcelTable)):
    chromosome=((re.search('Chromosome(.*)Fragment', MeanGenomeExcelTable.iat[fragment,0])).group(1))
    basepair1=((re.search('Basepairs(.*)-', MeanGenomeExcelTable.iat[fragment,0])).group(1))
    basepair2=((re.search('-(.*)', MeanGenomeExcelTable.iat[fragment,0])).group(1))
    James=[]
    for row in range(0,len(PotentialBiomarkersExcelTable)):
        if PotentialBiomarkersExcelTable.iat[row,2]<=(int(basepair2)) and PotentialBiomarkersExcelTable.iat[row,2]>=(int(basepair1))and (int(chromosome)-int(PotentialBiomarkersExcelTable.iat[row,1]))==0:
            James.append(PotentialBiomarkersExcelTable.iat[row,0])
    Biomarker.append(James)
mean=fragmentmean(MeanGenomeExcelTable)    
MeanGenomeExcelTable["Biomarker"]=Biomarker
MeanGenomeExcelTable["pattern_yes_no"]=pattern(MeanGenomeExcelTable,mean)

program_directory_to_output
#output file
MeanGenomeExcelTable.to_excel(r""+str(Window_Directory)+"/Output/BiomarkerMutationTables/"+str(cancer_type)+"MeanSamplesWGFMWithBiomarkersTable.xlsx", index = None, header=True)

#output file but data comes in 
Alpha_File=(r""+str(Window_Directory)+"/Output/BiomarkerMutationTables/"+str(cancer_type)+"MeanSamplesWGFMWithBiomarkersTable.xlsx")
Alpha1=pd.read_excel(io=Alpha_File, sheet_name="Sheet1")
Biomarkers=[]
MutationRates=[]
for row in range(0,len(Alpha1)):
    if Alpha1.iat[row,4]==1 and (Alpha1.iat[row,3] !="[]"):
        Biomarkers.append(Alpha1.iat[row,3])
        MutationRates.append(Alpha1.iat[row,2])
datax = {"MutationRate":MutationRates,"Biomarker":Biomarkers}#"Biomarker":Biomarkers
BiomarkerTable=pd.DataFrame(datax,index=MutationRates)


BiomarkerTable.to_excel (r""+str(Window_Directory)+"/Output/BiomarkerMutationTables/"+str(cancer_type)+"BiomarkerTypeAndBiomarkerLocationMutationTableI1.xlsx", index = None, header=True)
meantable(Window_Directory,cancer_type,1)
meantable(Window_Directory,cancer_type,2)
meantable(Window_Directory,cancer_type,3)
meantable(Window_Directory,cancer_type,4)
print("19.Completed mutation data processing.This took "+str(time.time()-StartTime)+" Seconds.")
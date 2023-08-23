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
print("yes")
path_to_program_directory=os.path.join((os.path.expanduser("~")),"Documents","Mutation_Scanner")
program_directory_to_samples=os.path.join(path_to_program_directory,"Input","Samples")
program_directory_to_reference_genomes=os.path.join(path_to_program_directory,"Input","Reference Genome")
program_directory_to_output=os.path.join(path_to_program_directory,"Output")

def chromosomexfragmentlist(FragmentSize,ChromosomeNumber):#Fragment List Generators for chromosomes.
    with open(os.path.join(program_directory_to_reference_genomes,f"Chr{ChromosomeNumber}.fa","r")) as text_file:
        text_data = text_file.read()
        list=re.sub("\n","",text_data)
        if ChromosomeNumber<10:
            chromosomesize=len(list[5:])
        else:
            chromosomesize=len(list[6:])
    Number_of_fragments=math.ceil(chromosomesize/FragmentSize)
    chromosomefragmentlist=[]
    for i in range(0,(Number_of_fragments),1):
        a=("Chromosome"+str(ChromosomeNumber)+"Fragment"+str(i)+",Basepairs "+str(i*FragmentSize+1)+"-"+str(i*FragmentSize+FragmentSize))
        chromosomefragmentlist.append(str(((a))))
    return(chromosomefragmentlist)

def WholeGenomeFragmentGenerator():#Creates a complete list of all chromosomes corresponding fragment lists.
    list=[chromosomexfragmentlist(FragmentSize,i) for i in range(1,23,1)]
    A=[item for sublist in list for item in sublist]
    return(A) 

def chromosomex_mutation_data(chromosomenumber,mutationlist):#strips irrelavant characters from mutation point.
    chromosomexlist=["0-1"]#0-1 is a filler
    for mutationposition in mutationlist:
        if (mutationposition[0:2] == str(chromosomenumber)):
            chromosomexlist.append(mutationposition[3:])
        elif(mutationposition[0:2] == (str(chromosomenumber)+":")):
            chromosomexlist.append(mutationposition[2:])
        
        else:
            continue
    Puremutationdatapoints=[int(mutationposition.split("-")[0])for mutationposition in chromosomexlist]
    return(Puremutationdatapoints) 

def Dictionary_Of_Fragment_mutation(FragmentSize,PurifiedMutationData,ChromosomeNumber,ChromosomesSizes):# Creates a Dictionary which has mutation count as one key and chromosome fragment location as the other.
    chromosomes={}
    ChromosomeSize=ChromosomesSizes[ChromosomeNumber-1]
    Number_of_fragments=math.ceil(ChromosomeSize/FragmentSize)
    
    for i in range(0,(Number_of_fragments)):
        a=("Chromosome"+str(ChromosomeNumber)+"Fragment"+str(i)+",Basepairs "+str(i*FragmentSize+1)+"-"+str(i*FragmentSize+FragmentSize))
        chromosomes.update({a:0})
    for mutation in PurifiedMutationData:
        if mutation==0:
            continue
        else:
            i=math.floor(mutation/1000000)
            a=("Chromosome"+str(ChromosomeNumber)+"Fragment"+str(i)+",Basepairs "+str(i*FragmentSize+1)+"-"+str(i*FragmentSize+FragmentSize))
            if chromosomes.get(a) ==0:
                    chromosomes.update({a:1})
            else:
                b=((chromosomes.get(a))+1)
                chromosomes.update({a:b})
    return(chromosomes)

def DictionaryToList(ABC):#Converts chromosome dictionary into chromosome mutation fragment list for dataframe addition. 
    C=[*ABC.values()]
    return(C)

def AllChromosomeMutationListPackage(FragmentSize,PureMutationData,ChromosomeSizes):#Processes the mutation data for every chromosome,plots each chromosome data mutation point to a chromosome fragment location in  a dictionary,Converts chromosome dictionary into chromosome mutation fragment list for dataframe addition,Adds each chromosome mutation fragment dictionary list conversion into one giant list file.. 
    AllChromosomeFragmentMutationList=[]
    for i in range(1,23,1):
        PurifiedMutationData=chromosomex_mutation_data(i,PureMutationData)
        AllChromosomeFragmentMutationDictionary=Dictionary_Of_Fragment_mutation(FragmentSize,PurifiedMutationData,i,ChromosomeSizes) 
        AllChromosomeFragmentMutationList+=(DictionaryToList(AllChromosomeFragmentMutationDictionary))
    return(AllChromosomeFragmentMutationList)    
   
def MutationPointDataCoding(excelsheet,SampleList,dictA,FragmentSize,ChromosomeSizes):#Each sample has its coding mutation data processed and added to a giant list of all sample mutation data.
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
        Samplesandmutationpoints.append(AllChromosomeMutationListPackage(FragmentSize,MutationGenomeLocationList,ChromosomeSizes))    
        #Complete=time.time()
        #print(Complete-Begin)
    return(Samplesandmutationpoints)    

def MutationPointDataNonCoding(excelsheet,SampleList,dictA,FragmentSize,ChromosomeSizes):#Each sample has its non coding mutation data processed and added to a giant list of all sample mutation data.
    Samplesandmutationpoints=[]
    ExcelGenomePositionColumnCoordinates=0
    for i in range(0,excelsheet.ncols):
        if excelsheet.cell_value(0,i)=="MutationPosition":
            ExcelGenomePositionColumnCoordinates+=(int(i))
    for i in range(0,len(SampleList)):
        SampleData=[SampleTablerowLocation for SampleTablerowLocation, Sample in dictA.items() if int(Sample)==SampleList[i]]
        MutationGenomeLocationList=[excelsheet.cell_value(Samplerows,ExcelGenomePositionColumnCoordinates) for Samplerows in SampleData]
        MutationGenomeLocationList.sort()
        Samplesandmutationpoints.append(AllChromosomeMutationListPackage(FragmentSize,MutationGenomeLocationList,ChromosomeSizes))   
    return(Samplesandmutationpoints)    

def IndividualSamples_ExcelCodingCreator(mutationdata,samplelist,FragmentSize,GenomeFragmentColumnData,BasePercentageColumnData,cancer_type,Window_Directory):# Coding sample mutation data presented into excel table with corresponding chromosome fragment location and base coverage percentages.
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

StartTime=time.time()   
FragmentSize=int(input("Enter genome fragment size."))
print(str(FragmentSize)+"bp selected as fragment size.")
cancer_type=input("Enter Cancer you wish to process.")
print(cancer_type+" selected as cancer type.")
ChromosomeSizes=[]   #This is a list which contains each chromosome base pair size and the corresponding code used to calculate this value.
for i in range(1,23):
    with open(os.path.join(program_directory_to_reference_genomes,f"Chr{ChromosomeNumber}.fa","r")) as text_file:
            text_data = text_file.read()
            list=re.sub("\n","",text_data)
            if i<10:
                ChromosomeSizes.append(len(list[5:]))
            else:
                ChromosomeSizes.append(len(list[6:]))
print("1.Chromosome Sizes Calculated.This took "+str(time.time()-StartTime)+" Seconds.")                             

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
    csv_cancer_mutation_data=csv.reader(Cancer_Mutation_Coding_Data)
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
# This variable NoFilesC is a number which represents the number of samples divided by one million
NoFilesC=((len(sample_id_coding))/(1000000))
#Using the above variable undergo a for loop in range of the NoFilesC variable rounded up. 
#This process will split up the large cancer file into subsequent small files which can be processed.
for IterationC in range(0,math.ceil(NoFilesC)):
    XC=(1000000*IterationC)
    YC=(1000000*(IterationC+1))
    #Creates a dictionary called FileData. 
    FileData= {"SampleId":sample_id_coding[XC:YC],"MutationPosition":mutation_genome_positions_coding[XC:YC]}
    # Dictionary FileData used to make a dataframe
    DFC=pd.DataFrame(FileData,index=(sample_id_coding[XC:YC]))
    DFC.to_excel((r""+str(Window_Directory)+"/Input/Samples/SmallFiles/"+str(cancer_type)+"_Coding_Sub"+str(IterationC+1)+".xlsx"), index = None, header=True)
print(len(sample_id_coding))
print(len(mutation_sample_unique_pairs_coding))
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

SampleFiles=([x[2] for x in os.walk(r""+str(Window_Directory)+"/Input/Samples/SmallFiles/")])[0]
CodingFileCountC=0
for i in SampleFiles:
    if str(cancer_type+"_Coding_Sub") in i:
        CodingFileCountC+=1

SampleFiles=([x[2] for x in os.walk(r""+str(Window_Directory)+"/Input/Samples/SmallFiles/")])[0]
CodingFileCountNC=0
for i in SampleFiles:
    if str(cancer_type+"_Non_Coding_Sub") in i:
        CodingFileCountNC+=1
print(CodingFileCountC)
totalsamples=len(SampleIdNC)+len(sample_id_coding)
if(totalsamples)>16380:
    print("Sheet Size not large enough.")
else:
    print("Sheet size large enough.")
for i in range(0,CodingFileCountC):
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
        
        GenomeFragmentColumn=WholeGenomeFragmentGenerator()#This generates a list of all genome fragments.
        print("6.Completed genome fragment list generation.This took"+str(time.time()-StartTime)+"Seconds.")
        BaseCoveragePercentageColumn=all_chromosomes_fragments_base_coverage(Window_Directory,FragmentSize)#This generates a list of all genome fragments base coverage.
        print("7.Completed chromosomes fragments base coverage percentage list.This took "+str(time.time()-StartTime)+" Seconds.")
        
        CodingData=MutationPointDataCoding(CodingExcelSheet,SampleListC,DictCSamplerows,FragmentSize,ChromosomeSizes)#This is the coding mutation data ready to be added to a dataframe.
        print("8.Completed creation of Coding Data.This took "+str(time.time()-StartTime)+" Seconds.")
        IndividualSamples_ExcelCodingCreator(CodingData,SampleListC,FragmentSize,GenomeFragmentColumn,BaseCoveragePercentageColumn,cancer_type,Window_Directory)#This processes the coding mutation data into a dataframe.
    else:
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
        NonCodingData=MutationPointDataNonCoding(NonCodingExcelSheet,SampleListNC,DictNCSamplerows,FragmentSize,ChromosomeSizes)#This is the non coding mutation data ready to be added to a dataframe.
        print("14.Completed creation of Non Coding Data.This took "+str(time.time()-StartTime)+" Seconds.")
        IndividualSamples_ExcelCodingAndNonCodingCreator(NonCodingData,SampleListNC,cancer_type,Window_Directory)
    print("Processed coding files")
for i in range(0,CodingFileCountC):
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
    NonCodingData=MutationPointDataNonCoding(NonCodingExcelSheet,SampleListNC,DictNCSamplerows,FragmentSize,ChromosomeSizes)#This is the non coding mutation data ready to be added to a dataframe.
    print("14.Completed creation of Non Coding Data.This took "+str(time.time()-StartTime)+" Seconds.")
    IndividualSamples_ExcelCodingAndNonCodingCreator(NonCodingData,SampleListNC,cancer_type,Window_Directory)#This processes the coding and non coding mutation data into a dataframe.
print("15.Completed creation of Coding and Non Coding Excel Sample Data File.This took "+str(time.time()-StartTime)+" Seconds.")
print("16.Comencing sample mutation average excel file creation.")
#This commences processing a new dataframe containing average mutation rate for all samples per genome fragment.
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
export_xlsx =MeanGenomeExcelTable.to_excel (r""+str(Window_Directory)+"/Output/WholeGenomeMutationTables/"+str(cancer_type)+"MeanSamplesWholeGenomeFragmentMutationTable.xlsx", index = None, header=True) 
print("17.Completed Creation of Mean Sample Whole Genome Fragment Mutation Rate Excel Table.This took "+str(time.time()-StartTime)+" Seconds.")
print("18.Biomarker Addition beginning.")

Data=MeanGenomeExcelTable.columns
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
MeanGenomeExcelTable.to_excel(r""+str(Window_Directory)+"/Output/BiomarkerMutationTables/"+str(cancer_type)+"MeanSamplesWGFMWithBiomarkersTable.xlsx", index = None, header=True)
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
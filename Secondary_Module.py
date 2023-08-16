# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 21:43:05 2021

@author: James
"""
import re
import math
def ChromosomeSequenceData(ChromosomeNumber,WindowDirectory):
    with open(r""+str(WindowDirectory)+"/Input/ReferenceGenome/chr"+str(ChromosomeNumber)+".fa") as text_file:
        text_data = text_file.read()
        listA=re.sub("\n","",text_data)
        if ChromosomeNumber<10:
            ChromosomeSequenceData=(listA[5:])
        else:
            ChromosomeSequenceData=listA[6:]
    return(ChromosomeSequenceData)#return

def FragmentBaseCoverage(FragmentNumber,FragmentSize,ChromosomeSequence):#Part 1 BPList
    FragmentSequence=(ChromosomeSequence[((FragmentNumber*FragmentSize)+1):((FragmentNumber*FragmentSize)+FragmentSize)])
    Bases_known=str(((FragmentSize-(FragmentSequence.count("N")+FragmentSequence.count("n")))/FragmentSize)*100)+"%"
    return(Bases_known)
 
def AllFragmentsOfChromosomeBaseCoverage(FragmentSize,ChromosomeSequence):#Part 1 BPList
    NumberOfFragments=math.ceil(len(ChromosomeSequence)/FragmentSize)
    list=[FragmentBaseCoverage(NumberOfFragments,FragmentSize,ChromosomeSequence)for NumberOfFragments in range(0,(NumberOfFragments),1)] 
    return(list)

def AllChromosomesFragmentsBaseCoverage(WindowDirectory,FragmentSize):
    Initial_list=[AllFragmentsOfChromosomeBaseCoverage(FragmentSize,ChromosomeSequenceData(i,WindowDirectory)) for i in range(1,23,1)]
    List_of_fragment_mutated_percentages=[item for sublist in Initial_list for item in sublist]
    return(List_of_fragment_mutated_percentages)

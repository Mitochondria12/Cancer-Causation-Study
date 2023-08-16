# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 17:06:46 2021

@author: James
"""

def DirectoryCreation(WindowDirectory):
    try:
        os.makedirs(WindowDirectory+'/Input/ReferenceGenome')
        os.makedirs(WindowDirectory+'/Input/Samples')
        os.makedirs(WindowDirectory+'/Output')
        Option_1=1
        print("Mutation Scanner files created,program will restart in twenty seconds. Please now enter your documents ,leave your sample data in the designated program file as as an excel document and add a complete fasta list of all human chromosomes to the reference genome file.")
    except:
        pass
    try:
        if Option_1==1 :
            sleep(20)
            os.execl(sys.executable, sys.executable, *sys.argv)
    except:
        pass
    if input("Please ensure the sample data is in the correct location before analysis begins, if it is press enter,else insert sample data and restart program.") == "":
        print("\n"+"Program starting.")
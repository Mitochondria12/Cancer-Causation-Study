# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 17:06:46 2021

@author: James
"""
import os
import sys
from time import sleep
import os,glob,errno,shutil,math,csv,pandas as pd,time,re,numpy as np,collections,getpass,sys,os.path
from time import sleep
from multiprocessing import Pool
from decimal import *

# Get the username to define the window directory

username =getpass.getuser()
window_directory="C:/Users/"+username+"/Documents/"


def create_directories(window_directory):
    
    # Indicates whether a restart is required
    restart_decision=True
    
    # Define the directories to be created
    directories = [
        window_directory + 'Mutation_Scanner/Input/ReferenceGenome',
        window_directory + 'Mutation_Scanner/Input/Samples',
        window_directory + 'Mutation_Scanner/Output'
    ]
    
    #Count of created subdirectories
    number_of_sub_folders=0
    for directory in directories:
        try:
            os.makedirs(directory)
            number_of_sub_folders+=1
        except FileExistsError: # Directory already exists; pass silently
            pass
        except Exception as e: # Handle any other exception
            print(f"An error occurred while creating {directory}: {e}")
            print("Ensure your main directory is on a windows os.")
            os.execl(sys.executable, sys.executable, *sys.argv)
            
    #This occurs when directories already exist.
    if number_of_sub_folders==0: 
        restart_decision=False
        
    # This occurs when no directories exist.    
    elif number_of_sub_folders==3:
        print("Genome Mutation Frequency Program directories created.")
        # Delays are added to give the user time to read the messages
        sleep(1.5)
        print("Please now add the relevant data to the directories.")
        sleep(1.5)
        print("This program will now restart.")
        sleep(1.5)
        
    #This occurs when some directories exist but not all.
    elif number_of_sub_folders>0 and number_of_sub_folders <3:
        print("Genome Mutation Frequency Program directories initally missing and now readded.")
        sleep(1.5)
        print("Please now add the relevant data to the directories.")
        sleep(1.5)
        print("This program will now restart.")
        sleep(1.5)

    return(restart_decision)

def check_data_presence(window_directory,directory_status=True):
    
    
    if directory_status==True:
        pass
    else:
            
        #Logic Gates
        data_missing=[]
        restart_decision=True
        
        # Define paths to the required data
        reference_genome_path = window_directory + 'Mutation_Scanner/Input/ReferenceGenome'
        samples_path = window_directory + 'Mutation_Scanner/Input/Samples'
    
        # Check if the reference genome files are present
        if not glob.glob(reference_genome_path + '/*.fa'):
            data_missing.append("Genome") #Reference Genome missing
    
        # Check if the sample data files are present (e.g., Excel files)
        if not glob.glob(samples_path + '/*.xlsx'):
            data_missing.append("Sample") # Sample Data missing
        
        #This occurs when your directorie has all the data required.
        if len(data_missing)==0:
            restart_decision=False
            
        #This occurs when your directories only has some of the data required.
        elif len(data_missing)==1:
            print("The directory has no {} data present, please add them to continue.".format(data_missing[0]))
            sleep(1)
            print("This program will now restart.")
            sleep(1)
            
        #This occurs when your directories are empty and need data added.
        elif len(data_missing)==2:
            print("Please now add your data files to the directory.")
            sleep(1)
            print("This program will now restart.")
            sleep(1)
        return restart_decision

def directory_creation(window_directory):
    # Create required directories and get the restart decision
    directories__restart_decision = create_directories(window_directory)

    # Check the presence of required data and get the restart decision
    data_restart_decision = check_data_presence(window_directory,directories__restart_decision)
    
    # Restart decision is not triggered because both data and directories are present.
    if directories__restart_decision == False and data_restart_decision ==False:
        pass 
    
    #Restart decision triggered because either or both data and directories are missing.
    else:
        sleep(5)
        os.execl(sys.executable, sys.executable, *sys.argv) # Restart the program

username = os.getlogin()
window_directory = "C:/Users/" + username + "/Documents/"
directory_creation(window_directory)

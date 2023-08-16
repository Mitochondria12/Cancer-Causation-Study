# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 17:06:46 2021

@author: James
"""
import os
import sys
import glob
from time import sleep
import getpass

def create_single_directory(directory):
    # Create a single directory and return True if created, False if it already exists
    try:
        os.makedirs(directory)
        return True
    except FileExistsError:
        return False
    except Exception as e:
        print(f"An error occurred while creating {directory}: {e}")
        print("Ensure your main directory is on a Windows OS.")
        sys.exit(1)

def print_with_delays(messages, delay=1.5):
    # Print messages with delays between them
    for message in messages:
        print(message)
        sleep(delay)

def create_directories(window_directory):
    # Create required directories and print messages based on their creation status
    directories = [
        window_directory + 'Mutation_Scanner/Input/ReferenceGenome',
        window_directory + 'Mutation_Scanner/Input/Samples',
        window_directory + 'Mutation_Scanner/Output'
    ]
    created_count = sum([create_single_directory(directory) for directory in directories])
    if created_count == 3:
        print_with_delays([
            "Genome Mutation Frequency Program directories created.",
            "Please now add the relevant data to the directories.",
            "This program will now restart."
        ])
        return True
    elif created_count > 0:
        print_with_delays([
            "Genome Mutation Frequency Program directories initially missing and now readded.",
            "Please now add the relevant data to the directories.",
            "This program will now restart."
        ])
        return True
    return False

def check_data_presence(window_directory, directory_restart):
    # Check if required data files are present in the directories
    if not directory_restart:
        reference_genome_path = window_directory + 'Mutation_Scanner/Input/ReferenceGenome'
        samples_path = window_directory + 'Mutation_Scanner/Input/Samples'
        data_missing = []
        if not glob.glob(reference_genome_path + '/*.fa'):
            data_missing.append("Genome")
        if not glob.glob(samples_path + '/*.xlsx'):
            data_missing.append("Sample")
        if data_missing:
            messages = [f"The directory has no {', '.join(data_missing)} data present, please add them to continue.",
                        "This program will now restart."]
            print_with_delays(messages)
            return True
    return False

def directory_creation(window_directory):
    # Main function to control directory creation and data presence check
    directory_restart = create_directories(window_directory)
    restart_required = directory_restart or check_data_presence(window_directory, directory_restart)
    if restart_required:
        sleep(5)
        os.execl(sys.executable, sys.executable, *sys.argv)  # Restart the program

username = os.getlogin()
window_directory = "C:/Users/" + username + "/Documents/"
directory_creation(window_directory)

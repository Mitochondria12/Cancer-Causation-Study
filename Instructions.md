# User Guide: Setting up and Running the Program

## Prerequisites
1. **Ensure Python is Installed**
    - Make sure Python is installed on your computer.

2. **Download Necessary Modules**
 
## Initial Setup
1. **Open the Primary Python Module File**
    - Locate and open the primary Python script to run the program.

2. **Run the Program**
    - Execute the program for the first time. This will create necessary directories in your `Documents` folder.

## Adding Reference Data
1. **Navigate to 'Reference' Folder**
    - Go to the `Reference` folder within the generated directories.

2. **Add Chromosome Autosome FASTA Files**
    - Ensure they resemble the example 'chromosome 22' FASTA folder in the `Resources` file.

## Collecting Cancer Data
1. **Visit the COSMIC Database**
    - Go to the COSMIC database and select your cancer type.

2. **Download Data**
    - Download exome and intron data for your selected cancer type.

3. **Save Files**
    - Save these files in the `Raw` directory. Name them as `Cancer_Coding.csv` and `Cancer_Non_Coding.csv`.

    > **Note**: Double-check that the columns match those in the example CSV files in `Resources`.

## Adding Blood Protein List
1. **Download Blood Proteins**
    - Download the list from the `Resources` file.

2. **Move File to 'Raw' Folder**
    - Place the downloaded file into the `Raw` folder within the `Input` directory.

## Running the Analysis
1. **Run the Program Again**
    - Execute the program.

2. **Input Parameters**
    - You will be prompted for the fragment size and the cancer type. Ensure you match the capitalization.

3. **Wait for Completion**
    - Wait for the analysis to be complete.

## Viewing Results
1. **Navigate to 'Output' Folder**
    - Go to the `Output` folder to find your results.

2. **Interpret Results**
    - Files like `Mutation_Frequency_Distributions` and `Average_Mutation_Frequency_Table` will be available.

3. **Graphing**
    - You can use this data to generate Excel graphs for visualizing mutation distribution across each chromosome.

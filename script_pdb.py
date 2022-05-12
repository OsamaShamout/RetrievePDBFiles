# -*- coding: utf-8 -*-


#################################################################
#                                                               #
#               COURSE: INTRO. TO SCRIPTIN CSC498P              #
#               INSTRUCTOR: DR. JOSEPH REBEHMED                 #
#               TASK: PYTHON FINAL ASSIGNMENT                   #
#               STUDENT: OSAMA SHAMOUT                          #  
#               LAU ID: 201906271                               #
#                                                               #
#################################################################    


"""
Python Final assignment file. Used to demonstrate various scripting codes to web scrape,
extract, and manipulate data using Python.
"""

#urllib library needed to download html documents for later usage.

##BeautifulSoup library imported to enable web scraping of websites. 
##Read more about library usage and documentation @ https://www.crummy.com/software/BeautifulSoup/bs4/doc/). Version upon creating document: 4.8.1.
##Note some libraries might not be automatically installed. Check IDE Envrionment and libraries documentation to understand how to run the library.
from bs4 import BeautifulSoup  

##Add system libraries to assist in file creation and relevant operating system commands.
import os

##Import regex
import re

import urllib

from pathlib import Path

import glob

import pandas as pd


import plotly.graph_objects as go
from plotly.subplots import make_subplots

import numpy as np

def check_file_exists(user_input):
    #Create list of files donwloaded
    files_downloaded = [f for f in os.listdir('.') if os.path.isfile(f)]
    
    print("\nChecking if file has been downloaded...")
    for f in files_downloaded:
        ## print(f)
        if (f == user_input):
            #print("File exists locally.")
            return True
    print("File does not exist locally.\n")
    return False

def generate_pdb (query=''):
    """
    Function takes in the user's inputted query 
    and generates PDB file as a .txt
    

    Parameters
    ----------
    query : STRING
        THE USER'S INPUTTED ID. The default is ''.

    Returns
    -------
    page : 
        DESCRIPTION.

    """
           
    ## Build URL of PDB File
    url = 'http://files.rcsb.org/view/' + query + ".pdb"  
    ## Get Parameters of Query
    params = {
    'from':"PDB_ID",
    'to':'ACC',
    'format':'txt',
    'query':query
    }
    ## Parse into Readable Format
    data = urllib.parse.urlencode(params)
    data = data.encode('ascii')
    ## Get from URL
    request = urllib.request.Request(url, data)
    try:
        print("\nChecking if file in the online database...")
        with urllib.request.urlopen(request) as response:
            ## Read Lines
            res = response.read()
            page=BeautifulSoup(res, features = "lxml").get_text()
            
            ## Check if pdb has XRAY
            xray_pdb = re.search("EXPDTA\s+X-RAY\s+DIFFRACTION", page)
            if not xray_pdb:
                return (None, 1)
            
            ## Check if protein atoms > 0
            atoms = re.search("REMARK\s+\d+\s+PROTEIN\s+ATOMS\s+:\s+(\d+)", page)
            num_atoms = int(atoms.group(1))
            if num_atoms == 0:
                return (None, 2)
            
            ## Split Lines into a List
            page=page.splitlines()
            
    except Exception:
        ## print("File not found in the online database.")
        return (None, 3)
    
    
    ## Save in the form of `query_id`.txt
    file_name = query + ".txt"
    with open(file_name, 'w') as f:
        for line in page:
            f.write("%s\n" % line)
    print("File found in online database.")
    return (page, 0)
      
def download_pdb(query, ignore_invalid = "n"):
    """
    Function that gets user input for the PDB ID and
    validates by alphanumeric of size 4
    &&
    Then downloads the pdb if it exists 

    
    Loops user until a correct input is given.

    Returns
    -------
    query : STRING
        THE USER'S INPUTTED ID.

    """
    pdb_create = False
    error = False
    invalid_input = False
    not_available = False
    
    if ignore_invalid=="y":
        ignore_invalid = True
    else:
        ignore_invalid = False
    

    while not pdb_create:
        ## Print message if an error is obtained
        if error:
            print(f"{query} is already downloaded.")
            error = False
            ## Get Input
            if not ignore_invalid:
                query = input("Enter PDB ID: ")
            else:
                return
            
        ## Print message if the user entered an invalid Query
        if invalid_input:
            print(f"{query} is not Valid: Must be 4 Alphanumeric Characters.")
            invalid_input = False
            ## Get Input
            if not ignore_invalid:
                query = input("Enter PDB ID: ")
            else:
                return
        ## Print message if the Query is available in the pdb
        if not_available:
            print(f"{query} is not Available in Protein Bank Database.")
            not_available = False
            ## Get Input
            if not ignore_invalid:
                query = input("Enter PDB ID: ")
            else:
                return
        #> For Debugging
        ## print("\nChecking if alhphanumeric condition satisfied...")
        
        ## Restrict input to be alphanumeric of 4 characters.
        expression = re.search("^([A-Za-z0-9]){4}$", query)
        
        ## If Query has invalid format, re-enter loop
        if not expression:
            #> print("Invalid ID. Enter 4 alphanumeric valid ID.")
            invalid_input = True
            ## Continue to restart loop
            continue
            
        ## Query has valid format.
        ## print("Valid Input")
        file_name = query + ".txt"
        ## Check if File exists
        if check_file_exists(file_name):
            pdb_create = False
            error = True
            continue
        ## File does not exist
        else:
            pdb_create = True 
          
        ## Generate the pdb 
        pdb, return_value = generate_pdb(query)
        if return_value == 1:
            print(f"{query} does not contain X-Ray Diffraction and was not downloaded.")
            return
        
        if return_value == 2:
            print(f"{query} has 0 atoms and was not downloaded.")
            return
        
        ## If no pdb is returned => pdb query is not in database
        if return_value == 3:
            pdb_create = False
            not_available = True
            continue
        else:
            print(f"{query}.txt downloaded.")

def read_input():
    """
    Reads the input from the user
    1- Asks user for the number of pdb files to download
    2- Downloads the pdb files.

    Returns
    -------
    None.

    """
    num_files = int(input("Enter the number of pdb files you would like to download: "))
    if num_files == 0:
        print("No files added.\n-------------")
        return
    
    ## Build the list of pdb ids
    pdb_ids = []
    for i in range(num_files):
        pdb_id = input(f"Enter pdb #{i+1}: ")
        pdb_ids.append(pdb_id)
    
    ignore_invalid = input("Do you wish to ignore Invalid IDs (y or n): ")
    print("------------------")
    for pdb_id in pdb_ids:
        print(f"Trying to download {pdb_id}...")
        download_pdb(pdb_id, ignore_invalid)
        print("------------------")

def get_data():
    """
    This method will read the pdb files
        and add the aminacids information to dataframe    

    Returns
    -------
    Dataframe of (aminaocid_cod, frequency, percentage)
    """
    ## --------
    ## Strategy:
    ## --------
    ## Data:
    ##    dataframe of (residue, frequency)
    ## Plan:
    ## 1. Go through each file
    ## 2. Go through each line
    ## 3. Regex for Atoms Signature
    ## 4. Add Data about new Residue
    
    columns = ["aminoacid", "frequency"]
    df = pd.DataFrame(columns = columns)
    
    dictionary_of_aa = {}
    
    for file in glob.glob("*.txt"):
        txt = Path(file).read_text()
        regex = re.compile('ATOM\s+\d+\s+CA\s+(\w{3})\s+')
        result = regex.findall(txt)
        for amino_acid in result:
            ## check if inside dictionary
            if amino_acid in dictionary_of_aa:
                dictionary_of_aa[amino_acid] = dictionary_of_aa[amino_acid] + 1
            else:
                dictionary_of_aa[amino_acid] = 1
    return dictionary_of_aa
        

def all_combined(amino_acids_frq):
    
    three_to_one = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                    'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
                    'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
                    'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    ## Calculate Total Frequency
    total_frequency = 0
    for key,value in amino_acids_frq.items():
        total_frequency += value
    
    ## Calculate Each Percentage
    percentages = []
    for key,value in amino_acids_frq.items():
        percentages.append(round(value/total_frequency*100, ndigits = 2))
    
    aa_list = []
    frq_list = []
    for key,value in amino_acids_frq.items():
        ## Change from 3 letters to one letter
        key = three_to_one[key]
        ## Add Aminoacids
        aa_list.append(key)
        ## Add Frequencies
        frq_list.append(value)
    
    ## Create Dataframe
    df = pd.DataFrame(columns = ["amino acid", "frequency", "percentage"])
    df["amino acid"] = aa_list
    df["frequency"]  = frq_list
    df["percentage"] = percentages
    
    ## Sort by percentage and return
    df = df.sort_values(by="percentage", ascending = False, ignore_index = True)
    return(df)
 

def generate_plot(df, by):
    
    aa  = df[by]
    frq = df["frequency"]
    prc = df["percentage"]
    
    # Create figure with secondary y-axis
    fig = make_subplots(specs=[[{"secondary_y": True}]])
    
    # Add traces
    fig.add_trace(
        go.Scatter(x=aa, y=frq, name="Frequency"),
        secondary_y=True,
    )
    
    fig.add_trace(
        go.Scatter(x=aa, y=prc, name="Percentages"),
        secondary_y=False,
    )
    
    # Add figure title
    fig.update_layout(
        title_text="Frequency & Percentages of Amino Acids"
    )
    
    # Set x-axis title
    fig.update_xaxes(title_text=by)
    
    # Set y-axes titles
    fig.update_yaxes(title_text="<b>Frequency<b>", secondary_y=True)
    fig.update_yaxes(title_text="<b>Percentages<b>", secondary_y=False)
    
    max_range = max(prc)*2
    if max_range > 100:
        max_range = 100
    fig.update_yaxes(range = [0,max_range], secondary_y=False)


    
    return fig
  
def by_category(amino_acids_frq):
    
    categories = {'CYS': 'Polar', 'ASP': 'Negatively Charged', 'SER': 'Polar',
                 'GLN': 'Polar', 'LYS': 'Positively Charged',
                 'ILE': 'Non-Polar', 'PRO': 'Non-Polar', 'THR': 'Polar',
                 'PHE': 'Non-Polar', 'ASN': 'Polar', 'GLY': 'Non-Polar',
                 'HIS': 'Positively Charged', 'LEU': 'Non-Polar',
                 'ARG': 'Positively Charged', 'TRP': 'Non-Polar', 
                 'ALA': 'Non-Polar', 'VAL':'Non-Polar', 'GLU': 'Negatively Charged',
                 'TYR': 'Polar', 'MET': 'Non-Polar'}
    
    ## Calculate Total Frequency
    cat_frequencies = {
            'Polar'              : 0,
            'Non-Polar'          : 0,
            'Positively Charged' : 0,
            'Negatively Charged' : 0}
    
    total_frequency = 0
    for key,value in amino_acids_frq.items():
        total_frequency += value
    
        cat_key = categories[key]
        cat_frequencies[cat_key] = cat_frequencies[cat_key]+value
    
    ## Calculate Each Percentage
    percentages = []
    for key,value in cat_frequencies.items():
        percentages.append(round(value/total_frequency*100, ndigits = 2))
    
    cat_list = []
    frq_list = []
    for key,value in cat_frequencies.items():
        ## Add Aminoacids
        cat_list.append(key)
        ## Add Frequencies
        frq_list.append(value)
    
    ## Create Dataframe
    df = pd.DataFrame(columns = ["category", "frequency", "percentage"])
    df["category"] = cat_list
    df["frequency"]  = frq_list
    df["percentage"] = percentages
    
    ## Sort by percentage and return
    df = df.sort_values(by="percentage", ascending = False, ignore_index = True)
    return(df)

def all_in_category(amino_acids_frq, category):    
    categories = {'CYS': 'Polar', 'ASP': 'Negatively Charged', 'SER': 'Polar',
                 'GLN': 'Polar', 'LYS': 'Positively Charged',
                 'ILE': 'Non-Polar', 'PRO': 'Non-Polar', 'THR': 'Polar',
                 'PHE': 'Non-Polar', 'ASN': 'Polar', 'GLY': 'Non-Polar',
                 'HIS': 'Positively Charged', 'LEU': 'Non-Polar',
                 'ARG': 'Positively Charged', 'TRP': 'Non-Polar', 
                 'ALA': 'Non-Polar', 'VAL':'Non-Polar', 'GLU': 'Negatively Charged',
                 'TYR': 'Polar', 'MET': 'Non-Polar'}
    
    ## Get Category Argument
    category_names = ["Non-Polar", "Polar", "Negatively Charged", "Positively Charged"]
    category_name = category_names[category-1]
    three_to_one = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                    'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
                    'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
                    'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    ## Calculate Total Frequency
    total_frequency = 0
    for key,value in amino_acids_frq.items():
        if categories[key] == category_name:
            total_frequency += value
    
    ## Calculate Each Percentage
    percentages = []
    for key,value in amino_acids_frq.items():
        if categories[key] == category_name:
            percentages.append(round(value/total_frequency*100, ndigits = 2))
            
    aa_list = []
    frq_list = []
    for key,value in amino_acids_frq.items():
        if not categories[key] == category_name:
            continue
        ## Change from 3 letters to one letter
        key = three_to_one[key]
        ## Add Aminoacids
        aa_list.append(key)
        ## Add Frequencies
        frq_list.append(value)
    
    ## Create Dataframe
    df = pd.DataFrame(columns = ["amino acid", "frequency", "percentage"])
    df["amino acid"] = aa_list
    df["frequency"]  = frq_list
    df["percentage"] = percentages
    
    ## Sort by percentage and return
    df = df.sort_values(by="percentage", ascending = False, ignore_index = True)
    return(df)

def one_aa(amino_acids_frq):
    
    ## Take User Input
    valid_codes = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M",
                   "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
    not_valid = True
    error_print = False
    while not_valid:
        if error_print:
            print("You Entered an Invalid AA 1 Letter Code")
        not_valid = False
        print("The Valid AA Codes Are: ",valid_codes)
        aa_code = input("Enter the AA One Letter Code: ")
        if aa_code not in valid_codes:
            not_valid = True
            error_print = True
    
    ## Calculate Frq & Percentages
    df = all_combined(amino_acids_frq)
    df = df.loc[df["amino acid"] == aa_code]
    return(df)

def calculate_interatomic_distances():    
    ## iad: interatmic distance
    dictionary_of_iad = {}
    proteins    = []
    iadistances = []
    
    for file in glob.glob("*.txt"):
        protein = Path(file).read_text()
        regex = re.compile('ATOM\s+\d+\s+CA\s+\w{3}\s+\w\s+\d+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)')
        result = regex.findall(protein)
        
        first = result[0]
        last = result[-1]
        
        first = [float(x) for x in first]
        last  = [float(x) for x in last]
        import math
        distance = round(math.sqrt((first[0]-last[0])**2 +
                                   (first[1]-last[1])**2 +
                                   (first[2]-last[2])**2),ndigits = 3)
        iadistances.append(distance)
        ## Get protein name
        protein_name = file.split(".")[0]
        proteins.append(protein_name)
    
    stats = {
        "minimum": min(iadistances),
        "maximum": max(iadistances),
        "median" : np.median(iadistances),
        "mean"   : np.mean(iadistances),
        "std"    : np.std(iadistances)
    }
    
    for key,value in stats.items():
        print(key,":",round(value,ndigits=3))
        
    import seaborn as sns
    sns.distplot(iadistances, hist = True,
                 axlabel = iadistances)
    
    import plotly.express as px
    fig = px.histogram(x=iadistances, y=proteins, labels={'x':'Interatomic Distance', 'y':'Protein ID'})
    fig.show()
    


        

def __main__():
    
    ## Add Files
    read_input()
    
    ## Get Dictionary of AA Frequency
    amino_acids_frq = get_data()
    ## print(amino_acids_frq)
    not_exit = True
    while not_exit:
        print("Options:\n1) ALL\n2) Per category\n3) Within category\n4) Specific AA\
    \nOther) Exit")
        choice = int(input("Choice: "))
        if choice == 1:
            ## Get Results
            df = all_combined(amino_acids_frq)
            ## Show Results DF
            print(df)
            ## Generate Plot
            fig = generate_plot(df, "amino acid")
            fig.show()
            
        elif choice == 2:
            df = by_category(amino_acids_frq)
            ## Show results
            print(df)
            ## Generate Plot
            fig = generate_plot(df, "category")
            fig.show()
        
        elif choice == 3:
            category = int(input("Enter number of category\n1. Non-polar\n\
2. Polar\n3. Negatively Charged\n4. Positively Charged\nChoice: "))
            df = all_in_category(amino_acids_frq, category)
            print (df)    
            ## Generate Plot
            fig = generate_plot(df, "amino acid")
            fig.show()
        elif choice == 4:
            df = one_aa(amino_acids_frq)
            print(df.to_string(index = False))
        else:
            not_exit = False
            print("Exiting")
    
    ## Calculate Interatomic Distances
    choice = input("Do you wish to calculate Interatomic Distances (y/n): ")
    if choice == "y":
        calculate_interatomic_distances()
    else:
        print("Done.")
    
__main__()   
        
    
    

    

    





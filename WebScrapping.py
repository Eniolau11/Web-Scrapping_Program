## Mini project

# Importing necessary libraries
from bs4 import BeautifulSoup  # Library for web scraping
import requests  # Library for sending HTTP requests
import random  # Library for generating random values
import subprocess  # Library for running shell commands
import pandas as pd  # Library for data manipulation and analysis
from Bio import SearchIO  # Biopython library for handling bioinformatics data
import matplotlib.pyplot as plt  # Library for data visualization
import os  # Library for interacting with the operating system
import seaborn as sns # Library for visualisation

# URL of the webpage containing worm data
worm_url = "https://parasite.wormbase.org/ftp.html"

# Sending an HTTP request to the webpage and creating a BeautifulSoup object
page = requests.get(worm_url)
soup = BeautifulSoup(page.text, 'html.parser')

# Extracting links ending with 'protein.fa.gz' from the webpage
protein_links = [link['href'] for link in soup.find_all('a', href=True) if link['href'].endswith('protein.fa.gz')]

# Selecting three random protein links
selected_protein_links = random.sample(protein_links, 3)

# Downloading and extracting the selected protein files
for link in selected_protein_links:
    file_url = link
    file_name = link.split("/")[-1]
    
    # Downloading the file
    with requests.get(file_url) as r:
        with open(file_name, 'wb') as file:
            for chunk in r.iter_content(chunk_size=8192):  # This was done to control the amount of RAM being used.
                file.write(chunk)

    # Extracting the downloaded file using subprocess (gunzip)
    subprocess.run(['gunzip', file_name])





# Using Pandas to parse through a CSV file and extracting rows with Accession starting with 'PF'
df = pd.read_csv('SearchResults-succinatedehydrogenase.tsv', sep='\t')
pf_df = df[df["Accession"].str.startswith('PF')]

# Downloading Pfam files using wget for each Accession in pf_df
for pf_accession in pf_df["Accession"]:
    pf_file_name = f'{pf_accession}.gz'
    command = f'wget https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/{pf_accession}?annotation=hmm -O {pf_accession}.gz'
    subprocess.run(command, shell=True)

    # Extracting the downloaded Pfam file using subprocess (gunzip)
    gunzip_command = f'gunzip {pf_file_name}'
    subprocess.run(gunzip_command, shell=True)






# Creating a shell script for HMMER (hmmer_script.sh)
hmmer_script = """#!/bin/bash
#SBATCH --job-name=hmmer_test
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --mem=8gb
#SBATCH --time=00:30:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=eu15@student.le.ac.uk

for pfam_file in PF*; do
    for fasta_file in *.fa; do
        output_file="${pfam_file}-${fasta_file%.fa}.out"
        hmmsearch --tblout "$output_file" -E 0.1 --noali "$pfam_file" "$fasta_file"
    done
done
"""

# Writing the HMMER script to a file
with open("hmmer_script.sh", "w") as script_file:
    script_file.write(hmmer_script)

#closing the file after writing
script_file.close()

#Making the shell script executable
subprocess.call("chmod +x hmmer_script.sh", shell=True)

#This will automatically run the shell script within the terminal remove the (#) to allow it to run
subprocess.call("./hmmer_script.sh", shell=True)




# Processing HMMER search results and storing them in a list (results)
results = []

for pfam_file in pf_df["Accession"]:
    file_prefix = f'{pfam_file}-'
    
    for file_name in os.listdir():
        if file_name.startswith(file_prefix) and file_name.endswith('.out'):
            with open(file_name, 'r') as file:
                for record in SearchIO.parse(file, 'hmmer3-tab'):
                    for hit in record.hits:
                        results.append({
                            'Accession': pfam_file,
                            'Query_Name': record.id,  
                            'Target_Name': hit.id,
                            'E-value': hit.evalue,
                            'Score': hit.bitscore,
                        })



# Creating a dataframe with the required values for each heading
summary_df = pd.DataFrame(results)
# Extract the first 3 characters from Target_Name. This was required to group the target names.
# I previously tried to group the target names but they have different end characters which didnt work.
summary_df['Target_Name_Short'] = summary_df['Target_Name'].str[:3]


print(summary_df)

# Exporting the DataFrame to a TSV file
csv_filename = 'summary_table.tsv'
summary_df.to_csv(csv_filename, sep='\t', index=False)
print(f'Table exported to {csv_filename}')


#Creating a bar chart showing the frquency of hits against Pfam domain grouped with the respected species.
summary_df.groupby(['Accession', 'Target_Name_Short']).size().unstack().plot(kind='bar', stacked=True)
plt.xlabel('Pfam Domain')
plt.ylabel('Frequency')
plt.title('Frequency of Hits for Each Pfam Domain')
plt.legend(title='Target Names', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.show()


#Box plot to show to distribution of E-values with against Pfam domain grouped with the respected species.
sns.boxplot(x='Accession', y='E-value', hue='Target_Name_Short', data=summary_df)
plt.yscale('log')  # Set y-axis to logarithmic scale for better visualization of E-values
plt.xlabel('Pfam Domain')
plt.ylabel('E-value')
plt.title('Distribution of E-values for Each Pfam Domain')
plt.legend(title='Target Names', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.show()
















    



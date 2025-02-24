# BOT563 Cui Amanita Pipeline

## 1 - Project information

### 1.1 - Project objectives

I am using the dataset that Cui et al. 2018 used to create the modern accepted phylogeny of Amanita. This is an important aspect of my project as I will use the Cui et al. phylogeny as a reference with which to put my own data into perspective.
The objectives of this project are:

- To recreate the phylogeny used by Cui et al. 2018 and begin to understand how to use their data
- To learn about phylogenetic pipelines and gain experience in coding and trouble shooting issues
- To build multilocus phylogenies
- To gain foundational knowledge in phylogenetic inference that I can apply to my own data in the future

### 1.2 - Information on the data set I am using

Fungal phylogenetics use similar loci as targets.
For the Amanita genus, it is common to see sequences that derive the Internal Transcribed Spacer (ITS) and the Large Ribosomal Subunit (nrLSU or 28S)
Furthermore, there are 3 other loci that are used intermittenly: Beta-Tubulin (btub), Rna Polymerase Subunit 2 (rpb2), and Transelongation Factor 1 - alpha (ef1a).
Using all five loci will give us more confidence in the phylogenetic tree that is created.
Cui et al. 2018 created a word document with all accessions used to create the 5 loci phylogenies. The data was derived from 864 specimens however a majority did not have data for all five loci with 1-4 of the loci not sequenced
I decided to target the specimens that had all five loci to parse down the size of the data set and to make it easier to work with the multiple loci.

- This gave me a set with 207 total specimens
I then went through and deleted specimens until there were only three specimens per species
- This gave me a set with 158 total specimens
While this brings the total number of sequences to (158x5=) 790 however I only have to work with one loci at a time so I don't think my computing time will be much longer than anyone else.

### 1.3 - Working with git

This is necessary for backing up my work in git and for Claudia to be able to follow my work along.

- Cloning git repo

```shell
git clone https://github.com/mdhensley7/BOT563_MH.git #pulls repository from github online and makes a copy in your directory
```

- Committing files to online git hub repo

```shell
git add . #stages the push
git commit -m "update readme" #commits the changes as long as all will merge properly
git push #sets changes and next step is creating a pull request on github online
```

### 1.4 - Packages I downloaded along the way

These will also be found in the code where they were being used but I downloaded some packages I did not use and just want to keep track of the fact that I have them

```shell
#Install the package needed to translate accessions nuc into .sra
sudo apt install sra-toolkit #I didn't end up using 'sra-toolkit'

#Package that downloads data via accession numbers
sudo apt install ncbi-acc-download #This worked for a data set I am no longer using but didn't work for the Cui data

pip install #Specific python code to install things

#The code below installed the Entrez Direct package which is what ultimately worked in dowloading the .fasta files
sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
sh -c "$(wget -q https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"

pip install biopython #Not used 
```

### 1.5 - Potentially useful scripts that ended up not being used in this pipeline

Taking the .gbk (SRA accession) files and converting them to .fasta

```shell
pip install biopython #needed for transition
for file in /mnt/c/Users/Michael/Desktop/BOT563/BOT563_MH/SDR_sraAccessions/*.gbk
do
python -c "from Bio import SeqIO; SeqIO.convert('$file', 'genbank', '${file%.*}.fasta', 'fasta');"
done 
#I now have a .gbk and a .fasta file for each succession in my SDR_sraAccessions
```

## 2 - Downloading and organizing the data

Figuring out how to take accession numbers from the word document and using those to download the .fasta files

### 2.1 - Copying accession numbers into .txt files

From the word document, I copied each accession number into a text file
Each loci has its own .txt with the respective accession numbers from each specimen
These files are found in the "cui_am5_raw" folder housed in BOT563_MH/data/

### 2.2 - Translating "loci_acc.txt" accessions into .fasta files

Each accession number is on its own seperate line in the .txt files. Below is a walk through of downloading the .fasta files from these accession numbers.

```shell
#First I move into the directory with the .txt files
cd /mnt/c/Users/Michael/Desktop/BOT563/BOT563_MH/data/cui_am5_raw

#to read the text files in the terminal
cat loci_acc.txt # swap "loci" in "loci_acc.txt" with name of loci in question
```

Downloading proper tools for process

```shell
# Installing entrez direct as the means of downloading fasta files
sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
sh -c "$(wget -q https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"
```

Executing 2.2

```shell
#This script automates the process so I don't have to download each individual fasta manually
accs=$(cat "loci_acc.txt") #swap "loci" with correct loci name
for i in $accs
do
efetch -db nuccore -format fasta -id $i #command to retrieve fasta file from corresponding accession number
done
#This worked by displaying all fasta files in the terminal but I don't see them downloaded anywhere
#Command below actually downloaded the fastas
efetch -db nuccore -format fasta -id $i > $i.fasta #uses same other commands as above just adds the command to download them into fasta files
#It worked!
```

Now we have a bunch of .fasta files in with the accession .txt files so need to organize

```shell
cd cui_am5_raw
mkdir fasta #only do once
cd fasta
mkdir nrlsu_fasta #do for each loci
cd .. #Should put you back into cui_am5_raw
mv *.fasta fasta #Moves all fasta files into fasta folder
cd fasta
mv *.fasta nrlsu_fasta #Moves all fasta files in fasta folder into nrlsu folder. Make sure to put correct loci in corresponding folder
```

- When downloading the "nrlsu" fastas there were four errors "missing argument describing data source".
  - 2 errors for btub
  - 1 error for its
  - 0 errors for rpb2
  - 0 errors for tef1a
    - However, all accessions downloaded a .fasta file so hoping there is no issue moving forward
- Repeated process above for the remaining four loci
  - At this point the ~/BOT563_MH/data/cui_am5_working/fasta folder should contain 5 separate folders each containing the .fasta files of the sequences corresponding to each loci
- There was also an empty .fasta file that was made for each loci so I deleted that one
  - each folder should contain only 158 fasta files
- Finally, I transfer the 'fasta' folder from cui_am5_raw to cui_am5_working

```shell
cd /mnt/c/Users/Michael/Desktop/BOT563/BOT563_MH/data/cui_am5_raw
mv fasta /mnt/c/Users/Michael/Desktop/BOT563/BOT563_MH/data/cui_am5_working
```

### From this point, I stick with only using the nrlsu data for class homework

## 3 - Aligning the sequences

Now I am going to try to align the nrlsu sequences in ~/BOT563_MH/data/cui_am5_working/fasta using MAFFT and Clustalw

This is a trial run of sorts so have to go back to make adjustments on the software's assumptions

-First thing I need to do is compile the nrlsu fasta folder into a single fasta file

```shell
cd
cat *.fasta > /mnt/c/Users/Michael/Desktop/BOT563/BOT563_MH/data/cui_am5_working/fasta/nrlsu.fasta
```

### 3.1 - Aligning via MAFFT

Make sure it is downloaded

```shell
conda install -c bioconda mafft
```

Running MAFFT

```shell
mafft --auto nrlsu.fasta > nrlsu_aligned #-- auto flag in the command ensures MAFFT selects the most appropriate alignment algorithm for a given set of sequences

mafft nrlsu_aligned #shows alignments in the terminal
```

### 3.2 - Aligning via ClustalW

nrlsu.fasta is the file I use for both alignments. Found in the "cui_am5_working/data" folder

```shell
grep ">" nrlsu.fasta | wc -l #checks how many sequences you have (should be 158)
clustalw2 -ALIGN -INFILE=nrlsu.fasta -OUTFILE=nrlsu_aligned_clustal.fasta -OUTPUT=FASTA
```

- This initial run gave an alignment score of 17900898
- Also produced a .dnd file which I am unsure about
- The alignment score will change with the parameters you manipulate
- Want the highest score but confirm this by reading the alignment manual

GO THROUGH DOCUMENTATION FOR THE ALIGNMENTS AND MANIPULATE THE PARAMETERS ACCORDINGLY

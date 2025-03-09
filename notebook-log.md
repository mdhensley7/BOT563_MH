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

- Finally just chose 23 samples from the set and 1 outgroup

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
cd /mnt/c/Users/Michael/Desktop/BOT563/BOT563_MH/data/3loci_24samp_working

#to read the text files in the terminal
cat loci_acc.txt # swap "loci" in "loci_acc.txt" with name of loci in question
```

Downloading proper tools for process

```shell
# Installing entrez direct as the means of downloading fasta files
bash -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
```

Executing 2.2

```shell
#Enter correct directory
cd /Users/hensley/Desktop/BOT563/BOT563_MH/data/3loci_24samp_working/3loci_acc
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
cd 3loci_24samp_working
mkdir fasta #only do once
cd fasta
mkdir lsu_fastas
mkdir rpb2_fastas
mkdir tef1a_fastas
cd .. #Should put you back into 3loci_24samp
mv *.fasta fasta #Moves all fasta files into fasta folder
cd fasta
mv *.fasta lsu_fasta #Moves all fasta files in fasta folder into nrlsu folder. Make sure to put correct loci in corresponding folder
```

```shell
cd 3loci_24samp
cat *.fasta > "loci".fasta"
#takes all .fasta files and compiles them into single file called ""loci".fasta"
#Be sure only the .fastas of the correct loci are in the folder. All .fastas will be compiled
```

At this point, the .txt files are in 3loci_24samp_working
this is where .fasta files. Once done, the .txt were moved to 3loci_acc

```shell
cd /3loci_24samp_working
mkdir 3loci_acc
mv *.txt 3loci_acc
```


```shell
cd /mnt/c/Users/Michael/Desktop/BOT563/BOT563_MH/data/cui_am5_raw
mv fasta /mnt/c/Users/Michael/Desktop/BOT563/BOT563_MH/data/cui_am5_working
```

### From this point, I stick with only using the nrlsu data for class homework

## 3 - Aligning the sequences

Now I am going to try to align the nrlsu sequences in ~/BOT563_MH/data/cui_am5_working/fasta using MAFFT and Clustalw

This is a trial run of sorts so have to go back to make adjustments on the software's assumptions

### 3.1 - Aligning via MAFFT

Make sure it is downloaded

```shell
conda install -c bioconda mafft
```

Running MAFFT

```shell
mafft --auto lsu.fasta > lsu_aligned_mafft #-- auto flag in the command ensures MAFFT selects the most appropriate alignment algorithm for a given set of sequences

mafft nrlsu_aligned_mafft #shows alignments in the terminal
```

### 3.2 - Aligning via ClustalW

nrlsu.fasta is the file I use for both alignments. Found in the "cui_am5_working/data" folder

#### Clustal W summary

| Software | Description | Strengths | Weaknesses | Assumptions | My Choices |
| -------- | ----------- | --------- | ---------- | ----------- | ---------- |
| Clustal W | Progressive MSA method boosting sensitivity through dynamic use of weights and costs | Creates an alignment that more accurately represents mutation types and the probability of their occurence in different areas of the genome | Still deals with the local minimum issue that comes with PMSAs and more options add more difficulty in choosing appropriate alignment parameters | Homologous sequences are evolutionarily related and in protein alignments, gaps do not occur randomly | TYPE=DNA GAPOPEN=1 GAPEXT=1

```shell
grep ">" lsu.fasta | wc -l #checks how many sequences you have (should be 24)
clustalw2 -ALIGN -INFILE=lsu.fasta -OUTFILE=lsu_aligned_clustal.fasta -OUTPUT=FASTA
```

- This initial run gave an alignment score of 17900898
- Also produced a .dnd file which I am unsure about
- The alignment score will change with the parameters you manipulate
- Want the highest score but confirm this by reading the alignment manual

## 4 - Creating distance and parsimony trees in R

Distance and parsimony are not really used for nything besides creating initial trees that bayesian and likelihood methods use as an initial tree
That said, still good to know how they work and are good for quick inference because their trees will tell you if there is a major problem with your data

### 4.1 - Distance method using ape

Install necessary packages

```shell
install.packages("adegenet", dep=TRUE)
library(ape) #This basically activates in current project
library(adegenet)
```

Loading data and plotting tree

```shell
dna <- fasta2DNAbin(file="~/lsu_aligned_clustal.fasta")
D <- dist.dna(dna, model="TN93") #computing genetic distances with Tamura and Nei 1993 model
tre <- nj(D) #getting the actual tree using neighbor-joining
tre <- ladderize(tre) #reorganizes tree to get a ladderized effect
plot(tre, cex=.6,)
```

Ths tree given had the outgroup placed within the tree which is not a good start. 
Used code below to reroot tree and specify the outgroup

```shell
root(tre, outgroup="KT833807.1", resolve.root = TRUE) #KT833807 is the accession number for the lsu loci of my outgroup Limacella
```

Now my tree is rooted properly however, I still am not sure if my data is okay if I have to specify the outgroup
I feel like it should be able to distinguish the outgroup.

### 4.2 - Maximum Parsimony using phangorn

| Software | Description | Strengths | Weaknesses | Assumptions | My Choices |
| -------- | ----------- | --------- | ---------- | ----------- | ---------- |
| Phangorn | Phangorn allows for the reconstruction of phylogenies using distance, parsimony, or likelihood. It also offers the possibility to estimate mixture and partition models | Seems to have more functionality than ape allowing for multiple different construction methods and the ability to compare them | Simplistic model of evolution, statistically inconsistent | Assumes rate of evolution is slow.

Installing needed packages

```shell
install.packages("adegenet", dep=TRUE)
install.packages("phangorn", dep=TRUE)
```

Activating packages

```shell
library(ape)
library(adegenet)
library(phangorn)
```

Loading in data

```shell
dna <- fasta2DNAbin(file=/data/3loci_24samp_working/fasta/lsu_aligned_clustal.fasta)
dna2 <- as.phyDat(dna)
```

Designating a starting tree to search tree space and compute parsimony score

```shell
tre.ini <- nj(dist.dna(dna,model="raw"))
parsimony(tre.ini, dna2)
```

Parisomony score:673

Searching for tree with maximum parsimony

```shell
tre.pars <- optim.parsimony(tre.ini, dna2)
```

Final pscore of 668 after 2 nni operations

Plot tree

```shell
plot(tre.pars, cex=0.6)
```

Tree again looks rather horrible. Tried to root tree using code in ape but still doesn't look great and didn't designate the correct outgroup
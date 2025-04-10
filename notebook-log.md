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
dna <- fasta2DNAbin(file="fasta/lsu_aligned_clustal.fasta")
D <- dist.dna(dna, model="TN93") #computing genetic distances with Tamura and Nei 1993 model
tre <- nj(D) #getting the actual tree using neighbor-joining
tre <- ladderize(tre) #reorganizes tree to get a ladderized effect
plot(tre, cex=.6,)
```

Ths tree given had the outgroup placed within the tree which is not a good start. 
Used code below to reroot tree and specify the outgroup

```shell
tre2 <- root(tre, outgroup="KT833807.1", resolve.root = TRUE) #KT833807 is the accession number for the lsu loci of my outgroup Limacella
plot(tre2, cex=0.6)
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
dna <- fasta2DNAbin(file="data/3loci_24samp_working/fasta/lsu_aligned_clustal.fasta)
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

Tree again looks rather horrible. So rooted it like the distance tree

```shell
tre.pars2 <- root(tre.pars, outgroup="KT833807.1", resolve.root = TRUE) #KT833807 is the accession number for the lsu loci of my outgroup Limacella
plot(tre.pars2, cex=0.6)
```

## 5 - Creating maximum likelihood trees using RAxML and IQTREE

### 5.1 - Quick ML notes

Maximum likelihood is a method of finding the probability of a tree given the data. Basically, we search the tree space to find the tree that has the highest probability of fitting the data.

The 4 most common things affecting ML
1. The starting tree
    - Maximum likelihood requires a starting tree to begin processing
    - Can be a random tree but lots of models will start with a tree generated with distance or parsimony
    - Also best practise to do multiple searches using multiple topologically different starting trees
2. The chosen model
    - Affects the actual tree landscape
    - May lead to identifiability
    - Might optimize the wrong function
3. The data itself
    - "garbage in, garbage out"
    - May lead to lack of signal
    - Also may lead to identifiability
4. Convergence
    - When do we estimate convergence and stop the search?

### 5.2 - Assumptions of ML

All ML models make these assumptions
1. The mutation process is the same for every branch of the tree
2. Sites evolve independently
3. All sites evolve the same
    - Some models can break this assumption by applying Γ model of rate heterogeneity

### 5.3 - For major influencers on ML performance

1. Starting tree selection
    - Affects optimization as possible to get stuck in poor likelihood region
    - Bad tree can slow process or end in a suboptimal tree
2. Model Chosen
    - Affects shape of treescape we optimize
    - Could optimize wrong function
    - Identifiability
3. Data
    - Lack of signal due to small sample size or poorly chosen sequence region
    - Identifiability
4. Convergence
    - When to terminate exploration of treescape?
    - Affects optimization

### 5.4 - iqtree overview

| Software | Description | Strengths | Weaknesses | Assumptions | My Choices |
| -------- | ----------- | --------- | ---------- | ----------- | ---------- |
| IQtree2 | Fast, stochastic algorithm for estimating ML phylogenies using bootstrap resampling to randomly reasses sites, rebuild tree, and provie branch support via a bootstrap value. |  1. Evaluates different substitution models to determine best fit for data. 2. Sotachistic algorithm usage helps escape local optima. 3. Allows for reversible and non-reversible models  | 1. More computationally and memory expensive. | Assumptions in 5.2 | Used GTR + GAMMA with n=1000 replicates like stated in the Cui et al 2018 paper |

### 5.5 - Downloading software

Both packages had a necessary online download

```shell
conda install -c bioconda iqtree
```

### 5.6 - Running iqtree2

1. Generating a tree with just the preset parameters

```shell
cd ~/BOT563_MH/data/3loci_24samp_working/fasta
iqtree2 -s lsu_aligned_clustal.fasta
```

ML score of -4480.739

2. Generating a tree with my chosen parameters

```shell
iqtree2 -s lsu_aligned_clustal.fasta -m GTR+G -alrt 1000 -nt AUTO
```
ML score of -4483.217 so actually worse than just using preset parameters!

3. Generating a tree like above but increasing iterations to 10,000

```shell
iqtree2 -s lsu_aligned_clustal.fasta -m GTR+G -alrt 10000 -nt AUTO
```
ML score of -4483.222 somehow got even worse... Maybe it's not that big of a deal if they are this close together and best to look at the trees.

Or maybe this shows an identifiability plateau?

### 5.7 - Tree visualization for iqtree

1. Primary way to just open the .treefile produced by iqtree in Figtree
    - Usually will stick with FigTree as it is a good tree visualization tool 

2. Can also visualize tree in R

```shell
library(ape) #Important to activate these packages even just for tree visualization
library(adegenet)
tre <- read.tree(text="parenthetical form")
#Access the parenthetical form by right clicking on the .bionj file go to "open file with" > textedit
#Copy the parenthetical including ";" and paste in the above code 
tre2 <- root(tre, outgroup="KT833807.1", resolve.root=TRUE)
plot(tre2, cex=.6)
```

### 5.8 - IQtree file organization

- Files from initial tree with preset parameters are found in ~/3loci_24samp_working/iqtree/preset
- Files from the third tree made are found in ~/3loci_24samp_working/iqtree/gtr+r_10000

### 5.9 - RAxML information (did not use for HW)

Keep in mind RAxML can only be run in the actual folder with the download (found Desktop/phylogenetics_software/maximum_likelihood)

## 6 - Bayesian phylogenetics

### 6.1 - Why bayesian?

- Bayesian inference allows us to refine Maximum Likelihood by implenting any previous knowledge we have about our data to influence the result
    - The previous knowledge is called the prior
- If we do not have any prior, we can still use what is called an uniformative prior.
    - Why not just stick with Maximum Likelihood if the prior is uninformative you may ask. Great question!
        - Bayesian phylogenetics also gives us an idea of how reliable the results are based on the variance of the "global" maximum the tree search comes to.

### 6.2 - What tools do I need?

- Beast2
    - This is likely the main tool you will use to create high quality tree so learn it well

- Mr. Bayes

### 6.3 Mr. Bayes

#### Software summary

| Software | Description | Strengths | Weaknesses | Assumptions | My Choices |
| -------- | ----------- | --------- | ---------- | ----------- | ---------- |
| MrBayes3 | Performs Bayesian inference and gives the user felxibility in model selection, choice of prior, and a diversity of parameters to manipulate. |  1. Flexibilty in parameter and prior selection 2. Model selection 3. Uses both MCMC and Metropolis-coupled MCMC; the latter of which allows for heated chains allowing for a more even traversal of tree space (the algorithm won't be as likely to get caught in a local maximum) | Extremely computationally and memory intensive; often need to send to computing servers/clusters  | Assume the prior given is accurate and appropriate to specific data set | I did one run with general class-set parameters and one run following parameters given in B.E. Wolfe paper linked below |

#### Installation process

-I don't know the difference but the first shell is the installation using the directions given at the mrbayes source github. The second method is defined in the 563 class directions

```shell
git clone --depth=1 https://github.com/NBISweden/MrBayes.git
cd MrBayes ##PATH coding_tools/phylogenetics_software/MrBayes
./configure
make && sudo make install
```

#### Mr. Bayes performs the analysis on an aligned .nexus file

Below is the code needed to convert aligned .fasta file into .nexus file (biopython) needed
- seqmagick is a handy dandy tool for quick and reliable file conversion for sequences
    - Super cool. it will also convert DNA->RNA RNA->protein or even DNA->protein
    - Just type "seqmagick convert" into terminal to see usages

```shell
conda install seqmagick #make sure you have it
seqmagick convert --output-format nexus --alphabet dna lsu_aligned_clustal.fasta lsu_aligned_clustal.nex
```
Now that we have the correct file format, we can move it into our "mrbayes" folder and run Mr Bayes
- We need to create a block that specifies all the parameters, priors, etc. for our Mr Bayes run
- Block is created in a separate .txt file

```shell
nano mbblock_general.txt
touch mbblock_general.txt
```
Below is the general format of the text going in the .txt file
- prset = setting up priors
- lset = defining substitution
- mcmp = setting for the MCMC
- ngen = should be 1,000,000 or more
- mcmc; = the command that actually runs the MCMC
- sumt; = Tells mrbayes to obtain a summary tree 

begin mrbayes;  
 set autoclose=yes;  
 prset brlenspr=unconstrained:exp(10.0);  
 prset shapepr=exp(1.0);  
 prset tratiopr=beta(1.0,1.0);  
 prset statefreqpr=dirichlet(1.0,1.0,1.0,1.0);  
 lset nst=2 rates=gamma ngammacat=4;  
 mcmcp ngen=1000000 samplefreq=10 printfreq=100 nruns=1 nchains=3 savebrlens=yes;  
 outgroup KT833807.1; <mark>specify your correct outgroup</mark>  
 mcmc;  
 sumt;  
end;  

##### Once .txt is inhabited by wanted code, append the file to end of the .nex "nexus" file

```shell
cat lsu_aligned_clustal.nex mbblock_general.txt > lsu_aligned_clustal_mb.nex
```
Finally, run the program

```shell
mb lsu_aligned_clustal_mb.nex
```

##### First Mr Bayes run stored in ~/mrbayes/mbrun1
- Has pretty good support values except for a early diverging node at 53 which is not super ideal
- Need to discuss more with Claudia to gain an understanding of these numbers

##### Wolfe, Pringle paper gives parameters used for an Amanita phylogeny inferred by Mr Bayes:
- Parameters used:
    - GTR+I+gamma
    - ngen = 10,000,000
    - 4 chains
    - Tree sampled every 1,000 generations
    - burnin = 2.5 million generations (default setting for Mr Bayes)
- <mark>Interesting note:</mark> 
    - ran phylogeney multiple times with different outgroups to make sure the tree didn't change

link to paper : https://www.researchgate.net/publication/229427712_The_Irreversible_Loss_of_a_Decomposition_Pathway_Marks_the_Single_Origin_of_an_Ectomycorrhizal_Symbiosis

Below is the new mrbayes block to copy the method used for building an already published Amanita phylogeny:

begin mrbayes;  
 set autoclose=yes;  
 prset brlenspr=unconstrained:exp(10.0);  
 prset shapepr=exp(1.0);  
 prset tratiopr=beta(1.0,1.0);  
 prset statefreqpr=dirichlet(1.0,1.0,1.0,1.0);  
 lset nst=6 rates=invgamma ngammacat=4;  
 mcmcp ngen=10000000 samplefreq=1000 printfreq=10000 nruns=1 nchains=4 savebrlens=yes;  
 outgroup KT833807.1; <mark>specify your correct outgroup</mark>  
 mcmc;  
 sumt;  
end;

- Append to .nex file using

```shell
cat lsu_aligned_clustal.nex mbblock_wolfe.txt > lsu_wolfe_mb.nex
```

- Run mrbayes

```shell
mb lsu_wolfe_mb.nex
```

## 7 - Coalescence

Coalescence to concatenation is what BBQ potato chips are to regular potato chips. Sure regular potato chips are fine and have their moments to shine but BBQ chips are just overall better, suit most ocassions, and will usually end up giving us more robust species trees.

Concatenation assumes that all genes within a species genome have undegone the same rate of evolution which is usually quite false. Incomplete Lineage Sorting (ILS), gene/genome duplication, and horizontal gene transfer are rather common occurances that lead to differences in how genes evolved.

The Multispecies Coalescent Models allows us to take the trees generated by a group of genes and infer the species tree.

### 7.3 - Handy-dandy tools to use

#### ASTRAL


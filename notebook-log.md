# BOT563 SDR Amanita Pipeline
## Project objectives
I want to create a pipeline that will organize the taxonomy of midwestern Amanita. It will only be using one locus of comparison so I will keep in mind that I may not have the most robust phylogeny. That is not the point of this. The point is to gain familiarity with a phylogenetic pipeline and begin the journey of tinkering with Amanita genetic sequences.
- The folder for this project is named "BOT563_MH"
## Cloning git repository
```shell
git clone https://github.com/mdhensley7/BOT563_MH.git #pulls repository from github online and makes a copy in your directory
```

## Pulling data from NCBI:
- I pulled a total of 785 accession numbers from NCBI of different specimens that were sequenced for their Internal Transcribed Spacer (ITS) region by Steven D. Russell.
- These accession numbers are found in the .txt file labelled "SDR_Accessions" in the project folder
- I basically manually copy and pasted these accessions from the website so will need to figure out a better way to strip data but fine for now

## Translating "SDR_Accession.txt" accessions into fast.q files
```shell
cat SDR_Accessions #reads texts file in the terminal

sudo apt install sra-toolkit #Install the package needed to translate accessions into fast.q
## This is for SRA accessions and I have nucleotide accession

sudo apt install ncbi-acc-download #Will hopefully allow me to download SRA accessions from my nucleotide accessions

#Using ncbi-acc-download to download SRA accessions from my nucleotide accessions
accs=$(cat SDR_Accessions)
for i in $accs
do
ncbi-acc-download -m 'nucleotide' $i
done
#These come out as .gbk files and are found in SDR_sraAccessions
```

The next step is taking the .gbk (SRA accession) files and converting them to fast-q
```shell
pip install biopython #needed for transition
for file in /mnt/c/Users/Michael/Desktop/BOT563/BOT563_MH/SDR_sraAccessions/*.gbk
do
python -c "from Bio import SeqIO; SeqIO.convert('$file', 'genbank', '${file%.*}.fasta', 'fasta');"
done 
#I now have a .gbk and a .fasta file for each succession in my SDR_sraAccessions
```
## Committing files to online git hub repo
```shell
git add . #stages the push
git commit -m "update readme" #commits the changes as long as all will merge properly
git push #sets changes and next step is creating a pull request on github online
```
## Creating a smaller subset of the data
Since this is a large data set, it will be slow to work with and will make troubleshooting a lot more time consuming.

First, I need to pull the .gbk files out of the same folder that my .fasta files are in.
```shell
mv /mnt/c/Users/Michael/Desktop/BOT563/BOT563_MH/SDR_sraAccessions/*.gbk /mnt/c/Users/Michael/Desktop/BOT563/BOT563_MH/sra_accessions # the * indicates that the command will be applied to any file with .gbk
#this command moved the .gbk files from SDR_sraAccessions to sra_accessions. The naming is very confusing so I will change the name of the first folder to "sdr_fasta"
mv -v SDR_sraAccessions sdr_fasta #now my folders are named properly
```

The "sdr_fasta" directory contains 785 fasta files. I need to take an n=50 sample

Seqkit is a tool that can randomly choose this sample for me but it needs to be installed
```shell
sudo apt install seqkit
```
I then combined all of my 785 .fasta files into a single .fasta file
```shell
cat *.fasta > /mnt/c/Users/Michael/Desktop/BOT563/BOT563_MH/sdr_fasta/sdr_am.fasta
#now I can use different scripts to inspect the data together
```
Manipulating the .fasta using seqkit
```shell
#Will show you the sequence of the first accession in the file; NOTE: changing the 1 to a 2 will show both the first AND second. Not just the second. -n 35 will show the first 35 sequences etc.
seqkit head -n 1 sdr_am.fasta

#converting the .fasta to a .fasta.gz
seqkit seq sdr_am.fasta -o sdr_am.fasta.gz

#pulling a random size 50 samples from the 785 accessions I have
zcat sdr_am.fasta.gz | seqkit sample -n 50 -o sdr_am_50samp.fasta.gz
#It worked! "sdr_am_50samp.fasta.gz" is the data I will use for the class

#Just to make sure it pulled randomly we are going to check the first 5 sequences of each
seqkit head -n 5 sdr_am.fasta.gz
seqkit head -n 5 sdr_am_50samp.fasta.gz
#They are indeed different suggesting a random pull of sequences
```
## Aligning the sequences
Now I am going to try to align the sequences in the file sdr_am_50samp.fasta.gz using MAFFT and potentially clustalw

This is a trial run of sorts so have to go back to make adjustments on the software's assumptions

### Aligning via MAFFT
Make sure it is downloaded
```shell
conda install -c bioconda mafft
```
Have to convert the file back into a fasta form
```shell
seqkit seq sdr_am_50samp.fasta.gz -o sdr_am_50samp.fasta
```
Running MAFFT
```shell
mafft --auto sdr_am_50samp.fasta > 50samp_aligned #-- auto flag in the command ensure MAFFT selects the most appropriate alignment algorithm for a given set of sequences

mafft 50samp_aligned #shows alignments in the terminal
```

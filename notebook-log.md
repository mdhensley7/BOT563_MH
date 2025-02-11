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
done #I now have a .gbk and a .fasta file for each succession in my SDR_sraAccessions
```
## Creating a smaller subset of the data
Since this is a large data set, it will be slow to work with and will make troubleshooting a lot more time consuming.

Instead I will pull out a subset of the data using the tool "seqkit"
```shell
sudo apt install seqkit
```


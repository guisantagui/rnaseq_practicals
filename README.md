# RNAseq analysis of *Pseudomonas aeurginosa data*

## Objective
Do alignment, differential expression analysis and functional enrichment of Pseudomonas aeruginosa RNAseq data: 6 strains (3 rhamnolipid producers and 3 non-producers).
Data from ["Evolution and regulation of microbial secondary metabolism"](https://elifesciences.org/articles/76119#s4).

## Workflow
1. Download FASTQs from SRA
2. Download PA14 reference genome and annotation.
3. Do quality control of the FASTQs.
4. Remove adapters and low quality reads with Trim Galore
5. Do quality control after trimming
6. Align FASTQs to reference (PA14 strain) with STAR.
7. Obtain gene counts
8. Perform differential expression analysis
9. Perform functional enrichment.

## Software requirements
- SRA-Tools
- STAR
- FastQC and multiqc
- Trim Galore
- R, and several packages
Most of these tools run only in UNIX-based systems (Linux and MacOS). So, this workflow needs to work in one of these OS.
Can be set up in Windows, but a bit tricky (VM, WSL), so we won’t do it today as it takes time.

## Setup
### SRA-tools
For downloading FASTQs from NCBI's [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra). There is more information about installation [here](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit
).
```bash
# Go to home directory
cd $HOME
# Download program
curl --output sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-mac64.tar.gz

# Uncompress file
tar -vxzf sratoolkit.tar.gz

# Add to path
# First check what shell you're using
echo $SHELL
# If /bin/bash
echo 'export PATH=$HOME/sratoolkit.3.2.1-mac-x86_64/bin:$PATH' >> ~/.bashrc
# If /bin/zsh
echo 'export PATH=$HOME/sratoolkit.3.2.1-mac-x86_64/bin:$PATH' >> ~/.zshrc

# Check installation
fasterq-dump –-version
```
### STAR
For performing the alignment. More information about it [here](https://github.com/alexdobin/STAR) and [here](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf
).

- Install git, if not already
In Mac.
```bash
# For Mac: If Homebrew is not installed
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

brew install git
git --version
```
In Linux
```bash
sudo apt install git-all
git --version
```
Download STAR
```bash
wget https://github.com/alexdobin/STAR/archive/2.7.11b.tar.gz
tar -xzf 2.7.11b.tar.gz

# Add to path
# First check what shell you're using
echo $SHELL
# If /bin/bash
echo 'export PATH=$HOME/STAR-2.7.11b/bin/MacOSX_x86_64/:$PATH' >> ~/.bashrc
# If /bin/zsh
echo 'export PATH=$HOME/STAR-2.7.11b/bin/MacOSX_x86_64/:$PATH' >> ~/.zshrc

# Check installation
star --version
```

### FastQC and MultiQC
Quality control of sequence data. More info [here](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [here](https://seqera.io/multiqc/).
```bash
# Install Homebrew (if not installed).
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install FastQC
brew install fastqc

# Check installation
fastqc –-version

# Install MultiQC
pip install multiqc

# Check installation
multiqc -–version
```
### Trim Galore
Remove bad quality reads and filter out adapters. More info [here](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/).
```bash
# Install cutadapt
pip install cutadapt

# Download Trim Galore
wget https://github.com/FelixKrueger/TrimGalore/archive/refs/tags/0.6.10.zip

# Add to path
# First check what shell you're using
echo $SHELL
# If /bin/bash
echo 'export PATH=$HOME/TrimGalore-0.6.10/:$PATH' >> ~/.bashrc
# If /bin/zsh
echo 'export PATH=$HOME/TrimGalore-0.6.10/:$PATH' >> ~/.zshrc
```

### Picard
For performing quality control (among other things) of the actual alignment.
```bash
brew install picard-tools
```

## Run pipeline
```bash
cd scripts
bash RNAseq_pipe.sh
```
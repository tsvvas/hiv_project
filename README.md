[![DudeItsCool](https://static.tildacdn.com/tild6464-3064-4237-a433-383539613333/bi_logo.png)](https://bioinf.me/)

# HIV evolution

## Introduction:
Human immunodeficiency virus (**HIV**) is a retrovirus that infects human cells, especially CD4+ T-cells. HIV infection is a dreadful disease that requires lifelong treatment. Although there have been developed several therapeutics that inhibit HIV replication, they are not curative. That's why dealing with this virus is important task. HIV has several mechanisms for escaping from human's immunity. One of them is high mutation rate. In this project we are to investigate the molecular evolution of HIV in time, changes in HIV's proteins due to mutations and their connection with escape from immune system.

## Data:
* Longitudinal HIV sequencing data was taken from https://hiv.biozentrum.unibas.ch/  
* Represantative proteomes were taken from  https://www.uniprot.org/, https://www.ebi.ac.uk/ena. You can download all represantative proteomes used in our project via https://drive.google.com/open?id=1Ghbuc4Tenl-fcvYZ6pyGLe-GHJ96jFjD (Be careful, download will start immediately).

## Problem:
HIV's high mutation rate is aimed to accommodate virus for exact human's immunity in order to escape it. HIV population typically mutates in the first few months after infection. The mutations are thoroughly studied and documented, but little focus was placed on the protein level. It is known that these mutations are preferentially found in CTL epitopes. But an interesting question is the fact whether these mutations are directed toward self, i.e., if they are aimed at becoming more like a host proteome, or not?

## Main hypothesis:
Main hypothesis of our project is that HIV mutates in direction toward its host to escape recognition. Although it's known that HIV's mutations occures in order to escape immune answer but we assume that this mutations is also providing similarity to human proteins.  

## Solution:
In order to validate the hypothesis we are going to measure the probability of HIV's peptide sequences to be similar to host (i.e. human). We derive the probability from the machine learning models trained to discriminate between human and pathogen peptide sequences.

## Methods:
First of all for each patient we've made evolution paths of HIV. For that reason single-patient phylogenetic trees for chosen region were created with pairwise alignments. Paths in those trees (from latest date to the reference) were taken as HIV's evolution paths. Further analyses include translating DNA into peptide sequences and applying quantitative and qualitative methods to them.  
In the quantitative analysis, frequencies of all possible amino acid k-mers (2-mers in our work) were calculated. To investigate the changes of similarity of HIV populations to self in time we trained a Random Forest to discriminate between human and pathogen protein sequences. Then HIV’s evolution paths were used to validate our hypothesis.  
In qualitative analysis translated sequences were used to calculate hydrophobicity metrics according to Kyle-Doolittle hydrophobicity score or kidera factors. Obtained data was also combined with virus’ paths to test our hypothesis.  

## Prerequisites
All preequisites can be found in **prerequisites.txt** at **help_dir**. Prequisites are written for anaconda's enviroment.

It's recommended to use computer with at least **2 GB RAM** and **2,00 GHz** CPU, especially while training machine learning models using our scripts.

## Results

**Quantitative analysis:**
Using quantitative method we had proved main hypothesis to be valid only for a part of the patients studied: patient p1, patient p4 (approximately), patient p11 (approximately). For another patients our model shows inconsistent results.
![2_mer_illu](https://github.com/tsvvas/hiv_project/raw/dev_andrew_classification/results/frequency_2_mer_plot/2_mer_plot.png)
We need more data to extensively validate the hypothesis or to eventually dismiss it. Also models in our pipeline should be improved to be more precise and sensitive.  

**Qualitative analysis:**
_TBD_

## Usage

Examples of how to use functions can be found in jupyter notebooks.  
With our code you can make:  

* single-patient phylogenetic trees with the structure like on the picture:  
![tree_illu](https://github.com/tsvvas/hiv_project/raw/dev_andrew_classification/results/examples/tree_example.png)
* download data from [this database](https://hiv.biozentrum.unibas.ch/)  
* train Random Forest on represantative proteomes (but proteomes you should download by yourself with the use of get_bact_fasta.sh from scripts)  
* use our explained pipeline  

## Roadmap

- Make library
- Improve accuracy and consistency of our methods
- Include analyzes for sequencing data from other papers (i.e. patients analyzed in this papers) about HIV

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## Acknowledgments
We want to show our apprecition to:
- Fabio Zanini and Co. for making helpful interesting [paper](https://elifesciences.org/articles/11282) and giving free access to their [database](https://hiv.biozentrum.unibas.ch)
- [Bioinformatics Institute](https://bioinf.me/en) for additional help in creating our team and other stuff

## Authors

- Vasily Tsvetkov, Institute of Bioorganic Chemistry
- Alexandra Ovsyannikova, Institute of Cyber Intelligence systems - National Research Nuclear University MEPhI
- Andrey Kravets, Moscow Institute of Physics and Technology

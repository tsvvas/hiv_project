<a href="http://google.com"><img src="https://static.tildacdn.com/tild6464-3064-4237-a433-383539613333/bi_logo.png" title="HIV" alt="HIV"></a>
<!-- [![DudeItsCool](https://static.tildacdn.com/tild6464-3064-4237-a433-383539613333/bi_logo.png)](http://google.com) -->

# HIV evolution

## Introduction:
HIV is still one of the remaining hard nuts to crack in XXI's century and it's epidemy is one of largest in the world. That's why dealing with this virus is important task. HIV has several mechanisms for escaping from human's immunity. One of them is high mutation rate. In this project we are to investigate the molecular evolution of HIV in time. For one group of patients we have found significant dependencies between their probability to escape and time. However, that result wasn't shown for all patient we have taken. All in all our approach shows unstable correlations, but work is in progress.

## Data:
(here i think we should show our data, with all links to it; also representative proteomes should be here i guess)

## Problem:
HIV's high mutation rate is aimed to accommodate virus for exact human's immunity in order to escape it. But connection between escaping and HIV's genome, which is changing in time due to mutations, is still unknown. Also why HIV chose exact evolution paths in exact human remains unrevealed. Solving this problem will help to predict how HIV can change in particular body, also it can help in creating stable vaccine for HIV.

## Solution:
In our work we've used proteins-based methods to analyze HIV's evolution.
First of all for each patient we've made evolution paths of HIV. For that reason single-patient phylogenetic trees for chosen region using pairwise alignments were created. Paths in those trees (from latest date to reference) were taken as HIV's evolution paths. Main analyzes include translating sequences from paths and then applying on them quantitative and qualitative analyzing methods.
In quantitative analyzes we've found frequencies of all possible aminoacid k-mers (2-mers in our work). To investigate HIV's accommodation for patient we've used representative proteomes for human, some bacteria and viruses to train classifier (human gene and anti-human gene) using machine learning approaches.
In qualitative method translated sequences were used to calculate hydrophobicity metrics according to Kyle-Doolittle table and then another machine learning algorithms were used. (there is part about hydrophobicity and Alexandra's methods; need some help here)

## Prerequisites
> jupyter notebook   
> Biopython  
> joblib  
> itertools  
> pandas  
> sklearn  
> urllib3  
> json  
> re  
> matplotlib.pyplot  
> shutil  
> os

## Getting started

To start one should clone this repository.


## Usage

Examples of how to use functions can be found in jupyter notebooks.

## Roadmap

- make library
- improve accuracy and consistency of our methods
- include analyzes for sequencing data from other papers about HIV

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## Acknowledgments
We want to show our apprecition to:
- Fabio Zanini and Co. for making helpful interesting [paper](https://elifesciences.org/articles/11282) and giving free access to their [database](https://hiv.biozentrum.unibas.ch)
- [Bioinformatics Institute](https://bioinf.me/en) for additional help in creating our team and other stuff

## Support

For the sake of saving time or in case of emergency, you can connect with us via [telegram](https://telegram.org/) @belsawan

## Authors

- Vasily Tsvetkov
- Alexandra Ovsyannikova
- Andrey Kravets

## Donations (Optional)

[![Support via Gratipay](https://cdn.rawgit.com/gratipay/gratipay-badge/2.3.0/dist/gratipay.png)](https://gratipay.com/spasibo_vasiliy_ne_obmanul/)

## License
[MIT](https://choosealicense.com/licenses/mit/)

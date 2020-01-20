# HIV evolution

In this project we investigate the molecular evolution of HIV in time. The longitudinal sequence data from HIV patients are obtained from [this database](https://hiv.biozentrum.unibas.ch). One should also consider to read [this paper](https://elifesciences.org/articles/11282) for additional info about HIV patients and how data was collected.

In our project we will analyze evolution paths of HIV for each patient. For finding out paths phylogenetic trees were created. Main analyzes include using quantitative and qualitative (hydrophobicity) composition of aminoacids from translated HIV sequences. 

For quantitative analyzes we have used along with HIV data represantative proteomes for human, bacteries and  viruses to train k-mer (2-mer) classifier (human gene and anti-human gene) using machine learning method Random Forest.  

For qualitive... (there is part about hydrophobicity and Alexandra's methods; need some help here)

All gathered results are stored in jupyter notebooks. Folder scripts contain code for creating phylogenetic trees, classificators, data preparations and etc.

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


## Usage

TODO

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

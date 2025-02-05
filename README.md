# Spistillata_nutrients
Repository for "Changes in host gene expression patterns underpin responses of the coral Stylophora pistillata to nutrient stress". Currently in review.

Experimental appraoch

# RNA-Seq mapping

Raw reads were de-multiplexed, trimmed, and adapters removed prior to data being released from Novogene (Cambridge, UK). Reads were aligned to the _[Stylophora pistillata](http://spis.reefgenomics.org/)_ host genome assembly using HISAT2 (_v_ 2.2.1) ([Kim et al., 2019](https://www.nature.com/articles/s41587-019-0201-4)). Assembly and count matrix generation were done using StringTie (_v_ 2.2.1) ([Pertea et al., 2015](https://www.nature.com/articles/nbt.3122)). Mapping precision of the generated, merged GFF files to the _S. pistillata_ reference assembly was assessed using GFFcompare (v 0.12.6) ([Pertea & Pertea, 2020](https://f1000research.com/articles/9-304)).

# Functional annotation
Functional annotation of protein sequences was done using Diamond (_v_ 2.1.8) ([Bunchfink et al., 2015](https://www.nature.com/articles/nmeth.3176)), Uniprot ([Apweiler et al., 2004](https://academic.oup.com/nar/article/51/D1/D523/6835362)), and eggNOG ([Cantalapiedra et al., 2021](https://www.biorxiv.org/content/10.1101/2021.06.03.446934v2)).

# Differential expression analysis
Differential expression analysis of the coral host was conducted using Bioconductor edgeR (_v_ 3.2) ([Robinson et al., 2010](https://academic.oup.com/nar/article/53/2/gkaf018/7973897?login=false)) in R considering phosphate and nitrate as factors. 


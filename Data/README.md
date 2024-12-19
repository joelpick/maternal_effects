# Data README

**Project title:** Simple maternal effect animal models provide biased estimates of additive genetic and maternal variation

**Data License:** CC-BY 4.0


## Raw/

Directory containing raw data

### Bonnet/
Directory containing data files from Bonnet et al 2022 Science.
Original data available at: https://doi.org/10.1126/science.abk0853

Directory contains a file for each of the 19 populations used in the study: bsR = bighorn sheep on Ram Mountain, btM = blue tits at Muro, btP = blue tits at Pirio, btR = blue tits at la Rouvi`ere, cfG = collared flycatchers on Gotland, gtH = great tits in Hoge Veluwe, gtW = great tits in Wytham Woods, hhT = hihi on Tiritiri Matangi Island, hhK = hihi at Karori, mkK = meerkats in the Kalahari, rdR = red deer on the Isle of Rum, rmC = rhesus macaques at Cayo Santiago, rsK=red squirrels in Kluane, sfC = superb fairy-wrens in Canberra, shN = spotted hyenas in the Ngorongoro Crater, spM= song sparrows on Mandarte Island, ssS = Soay sheep on St Kilda, ,svG = snow voles in Graubunden, ybA = yellow baboons at Amboseli. 

Each file contains the variable:

dam: mother ID


### BT_Pick.csv
Pedigree from Pick et al 2022 Nature Ecology and Evolution
Original data repository: [https://doi.org/10.5281/zenodo.6323689](https://doi.org/10.5281/zenodo.6323689)

**Variables:**

animal: individual ID

dam: mother ID

sire: father ID


### covariances.csv
Relatedness of different individuals to each other that makes up the covariance between different relatives

**Variables:**

relative: The relationship to a focal individual

relative_mother: The relationship to a focal individuals mother

r1: relatedness between a focal individual and the relative

r2: relatedness between a focal individuals mother and the relatives mother

r3 f-rm: relatedness between a focal individual and the relatives mother

r4 r-m: relatedness between a focal individuals mother and the relative


### MGE - all_est.csv
Data extracted from literature on maternal genetic variances

**Variables:**

Study: Study ID (see Supplementary Table S1 to match IDs to studies)

Species: Common name for species

Population: Whether or no the population is wild (levels 'wild', 'breeding design' and 'managed')

Trait: Phenotype measured

juv_trait: Whether or not a trait is classed as a juvenile trait (1-yes, 0-no)

h2_1: heritability from simple maternal effect model

c2_1: proportion of variance due to maternal id effects from simple maternal effect model

h2_2: heritability from full maternal effect model

m2_2: proportion of variance due to maternal genetic effects from full maternal effect model

c2_2: proportion of variance due to maternal environmental effects from full maternal effect model

Source: where in the paper the effects sizes were extracted from


### moore_2019_datawithSE.csv
Data from Moore et al. 2019 Ecology Letters. 
Original data repository: [https://doi.org/10.5061/dryad.360v97q](https://doi.org/10.5061/dryad.360v97q)

**Variables:**

StudyCode: ID of the study

h2: Heritability estimate

h2.samp.err: Standard error of heritability estimate


### postma2014_SM.csv
Data from POstma 2014 Quantitative Genetics in the Wild
Original data repository: [https://global.oup.com/booksites/fdscontent/booksites/uk/booksites/content/9780199674244/supplementary/Supplementary_material_Chapter2.xls](https://global.oup.com/booksites/fdscontent/booksites/uk/booksites/content/9780199674244/supplementary/Supplementary_material_Chapter2.xls)

**Variables:**

study ID: ID of the study

method: method of estimating heritability

h2: Heritability estimate

SE of h2: Standard error of heritability estimate


### young_postma_2023_data.csv
Data from Young and Postma 2023
Original data repository: [https://doi.org/10.34894/RIVFHW](https://doi.org/10.34894/RIVFHW)

**Variables:**
estimate.ID: ID of the effect size

study.ID: ID of the study

method: method of estimating heritability

h2: Heritability estimate

SE.h2: Standard error of heritability estimate


### RD_Gauzere.txt
Pedigree from Gauzere et al. 2020
Original data repository: [https://doi.org/10.5061/dryad.gtht76hhs](https://doi.org/10.5061/dryad.gtht76hhs)

**Variables:**

ID: individual ID

Sire: mother ID

Dam: father ID




## Intermediate/

Contains processed data files

### mge_sims_small_ped.Rdata

### mge_sims3.Rdata

### parameters.Rdata

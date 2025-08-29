
# All code and generated data for manuscript *Simple maternal effect animal models may provide biased estimates of additive genetic and maternal variation*
**Authors:** Joel Pick, Craig Walling and Loeske Kruuk

**Abstract:** 
Maternal effects (the consistent effect of a mother on her offspring) can inflate estimates of additive genetic variation ($V_A$) if not properly accounted for. As they are typically assumed to cause similarities only among maternal siblings, they are often accounted for by modelling maternal identity effects. However, if maternal effects have a genetic basis, they create additional similarities among relatives with related mothers that are not captured by maternal identity effects. Unmodelled maternal genetic variance ($V_{Mg}$) may therefore still inflate $V_A$ estimates in common quantitative genetic models, which is underappreciated in the literature. Using published data and simulations, we explore the extent of this problem. Published estimates from eight species suggest that a large proportion of total maternal variation ($V_M$) is genetic ($\sim$65\%). Both these data and simulations confirmed that unmodelled $V_{Mg}$ can cause overestimation of $V_A$ and underestimation of $V_M$, the bias increasing with the proportion of non-sibling maternal relatives in a pedigree. Simulations show these biases are further influenced by the size and direction of any direct-maternal genetic covariance. The estimation of total additive genetic variation ($V_{A_t}$; the weighted sum of $V_A$ and $V_{Mg}$) is additionally affected, limiting inferences about evolutionary potential from simple maternal effect models. Unbiased estimates require modelling $V_{Mg}$ explicitly, but these models are often avoided due to perceived data limitations. We demonstrate that estimating $V_{Mg}$ is possible even with small pedigrees, reducing bias in $V_A$ estimates and maintaining accuracy in estimates of $V_A$, $V_M$, and $V_{A_t}$. We therefore advocate for the broader use of these models.

**Pre-print DOI:** [https://ecoevorxiv.org/repository/view/8226/](https://doi.org/10.32942/X2V33J)

**Funders:** Royal Society (RSRP-R1-211017) and ERC (101020503 — EVOWILD — ERC-2020-ADG)

**Data License:** CC-BY 4.0

**Code License:** MIT


## Repository organisation Structure

### Data

- Contains two folders
    - Raw with all raw data files
    - Intermediate with all processed files

### extract

- Contains files needed to re-produce data extraction from figures

### Figures

- Contains all generated figures


### R

- Contains all R code

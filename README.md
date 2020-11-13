# treemut

### Overview 
Code implementing maximum likelihood assignment of mutations to trees.  The method uses an EM method to soft assign mutations to branches and simultaneously estimate branch lengths.

### Algorithm
Each mutation is assigned to branch with probability:

P(Mutation has genotype ) ~ P(Read Count | genotype)P(genotype)

Where a given branch implies a genotype and P(genotype) is proportional to branch length.  The read counts are modelled using the binomial distribution with up to 4 distinct per sample specific error rates and an assumed VAF=0.5 for variant sites and VAF=0 for wild type sites.

## Installation
```r
library(devtools)
install_git("https://github.com/NickWilliamsSanger/treemut")
```
## View vignette for usage 
```r
browseVignettes(“treemut”)
```



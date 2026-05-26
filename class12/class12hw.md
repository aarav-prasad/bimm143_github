# Class 12 HW: Population Scale Analysis
Aarav Prasad (PID: A17440940)

## Population Scale Analysis

One sample is obviously not enough to know what is happening in a
population. You are interested in assessing genetic differences on a
population scale. So, you processed about ~230 samples and did the
normalization on a genome level. Now, you want to find whether there is
any association of the 4 asthma-associated SNPs (rs8067378…) on ORMDL3
expression.

> Q13: Read this file into R and determine the sample size for each
> genotype and their corresponding median expression levels for each of
> these genotypes

``` r
expr <- read.table("rs8067378_ENSG00000172057.6.txt")
head(expr)
```

       sample geno      exp
    1 HG00367  A/G 28.96038
    2 NA20768  A/G 20.24449
    3 HG00361  A/A 31.32628
    4 HG00135  A/A 34.11169
    5 NA18870  G/G 18.25141
    6 NA11993  A/A 32.89721

``` r
table(expr$geno)
```


    A/A A/G G/G 
    108 233 121 

``` r
b <- boxplot(exp ~ geno, data = expr)
```

![](class12hw_files/figure-commonmark/unnamed-chunk-3-1.png)

``` r
b$stats[3, ]
```

    [1] 31.24847 25.06486 20.07363

> Q14: Generate a boxplot with a box per genotype, what could you infer
> from the relative expression value between A/A and G/G displayed in
> this plot? Does the SNP effect the expression of ORMDL3?

``` r
library(ggplot2)
```

``` r
ggplot(expr) + aes(geno, exp, fill=geno) +
  geom_boxplot(notch=TRUE)
```

![](class12hw_files/figure-commonmark/unnamed-chunk-5-1.png)

Based on the boxplot, ORMDL3 expression appears to differ by genotype.
A/A individuals show the highest expression, G/G individuals show the
lowest expression, and A/G individuals fall in between. This suggests
rs8067378 may affect ORMDL3 expression, with the G allele associated
with lower expression.

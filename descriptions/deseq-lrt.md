Runs DESeq2 using LRT (Likelihood Ratio Test)
=============================================

The LRT examines two models for the counts, a full model with a certain number of terms and a reduced model, in which some of the terms of the full model are removed. The test determines if the increased likelihood of the data using the extra terms in the full model is more than expected if those extra terms are truly zero.

The LRT is therefore useful for testing multiple terms at once, for example testing 3 or more levels of a factor at once, or all interactions between two variables. The LRT for count data is conceptually similar to an analysis of variance (ANOVA) calculation in linear regression, except that in the case of the Negative Binomial GLM, we use an analysis of deviance (ANODEV), where the deviance captures the difference in likelihood between a full and a reduced model.

When one performs a likelihood ratio test, the p values and the test statistic (the stat column) are values for the test that removes all of the variables which are present in the full design and not in the reduced design. This tests the null hypothesis that all the coefficients from these variables and levels of these factors are equal to zero.

The likelihood ratio test p values therefore represent a test of all the variables and all the levels of factors which are among these variables. However, the results table only has space for one column of log fold change, so a single variable and a single comparison is shown (among the potentially multiple log fold changes which were tested in the likelihood ratio test). This indicates that the p value is for the likelihood ratio test of all the variables and all the levels, while the log fold change is a single comparison from among those variables and levels.

**Technical notes**

1. At least two biological replicates are required for every compared category
2. Metadata file describes relations between compared experiments, for example
    
    ```
    ,time,condition
    DH1,day5,WT
    DH2,day5,KO
    DH3,day7,WT
    DH4,day7,KO
    DH5,day7,KO
    ```
    where `time, condition, day5, day7, WT, KO` should be a single words (without spaces) and `DH1, DH2, DH3, DH4, DH5` correspond to the experiment aliases set in **RNA-Seq experiments** input.
3. Design and reduced formulas should start with **~** and include categories or, optionally, their interactions from the metadata file header. See details in DESeq2 manual [here](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#interactions) and [here](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#likelihood-ratio-test)
4. Contrast should be set based on your metadata file header and available categories in a form of `Factor Numerator Denominator`, where `Factor` - column name from metadata file, `Numerator`  - category from metadata file to be used as numerator in fold change calculation, `Denominator` - category from metadata file to be used as denominator in fold change calculation. For example `condition WT KO`.
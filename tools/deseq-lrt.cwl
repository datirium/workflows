cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/scidap-deseq:v0.0.22


inputs:

  expression_files:
    type: File[]
    inputBinding:
      position: 5
      prefix: "--input"
    doc: "Grouped by gene / TSS/ isoform expression files, formatted as CSV/TSV"

  expression_file_names:
    type: string[]
    inputBinding:
      position: 6
      prefix: "--name"
    doc: "Unique names for input files, no special characters, spaces are allowed. Number and order corresponds to --input"

  metadata_file:
    type: File
    inputBinding:
      position: 7
      prefix: "--meta"
    doc: "Metadata file to describe relation between samples, where first column corresponds to --name, formatted as CSV/TSV"

  design_formula:
    type: string
    inputBinding:
      position: 8
      prefix: "--design"
    doc: "Design formula. Should start with ~. See DeSeq2 manual for details"

  reduced_formula:
    type: string
    inputBinding:
      position: 9
      prefix: "--reduced"
    doc: "Reduced formula to compare against with the term(s) of interest removed. Should start with ~. See DeSeq2 manual for details"

  contrast:
    type:
      - string
      - string[]
    inputBinding:
      position: 10
      prefix: "--contrast"
    doc: "Contrast to be be applied for output, formatted as Factor Numerator Denominator or 'Factor Numerator Denominator'"

  output_prefix:
    type: string?
    inputBinding:
      position: 9
      prefix: "--output"
    doc: "Output prefix for generated files"

  threads:
    type: int?
    inputBinding:
      position: 10
      prefix: "--threads"
    doc: "Threads number"


outputs:

  diff_expr_file:
    type: File
    outputBinding:
      glob: "*_table.tsv"

  ma_plot:
    type: File
    outputBinding:
      glob: "*_ma_plot.png"

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: [run_deseq_lrt.R]
stdout: deseq_stdout.log
stderr: deseq_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


s:name: "DESeq2 (LRT) - differential gene expression analysis using likelihood ratio test"
label:  "DESeq2 (LRT) - differential gene expression analysis using likelihood ratio test"
s:alternateName: "Differential gene expression analysis based on the LRT (likelihood ratio test)"


s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/deseq-lrt.cwl
s:codeRepository: https://github.com/Barski-lab/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Cincinnati Children's Hospital Medical Center"
  s:location:
  - class: s:PostalAddress
    s:addressCountry: "USA"
    s:addressLocality: "Cincinnati"
    s:addressRegion: "OH"
    s:postalCode: "45229"
    s:streetAddress: "3333 Burnet Ave"
    s:telephone: "+1(513)636-4200"
  s:logo: "https://www.cincinnatichildrens.org/-/media/cincinnati%20childrens/global%20shared/childrens-logo-new.png"
  s:department:
  - class: s:Organization
    s:legalName: "Allergy and Immunology"
    s:department:
    - class: s:Organization
      s:legalName: "Barski Research Lab"
      s:member:
      - class: s:Person
        s:name: Michael Kotliar
        s:email: mailto:misha.kotliar@gmail.com
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898

doc: |
  Runs DESeq2 using LRT (Likelihood Ratio Test)

  The LRT examines two models for the counts, a full model with a certain number of terms and a reduced model,
  in which some of the terms of the full model are removed. The test determines if the increased likelihood of
  the data using the extra terms in the full model is more than expected if those extra terms are truly zero.

  The LRT is therefore useful for testing multiple terms at once, for example testing 3 or more levels of a factor
  at once, or all interactions between two variables. The LRT for count data is conceptually similar to an analysis
  of variance (ANOVA) calculation in linear regression, except that in the case of the Negative Binomial GLM, we use
  an analysis of deviance (ANODEV), where the deviance captures the difference in likelihood between a full and a
  reduced model.

  When one performs a likelihood ratio test, the p values and the test statistic (the stat column) are values for
  the test that removes all of the variables which are present in the full design and not in the reduced design.
  This tests the null hypothesis that all the coefficients from these variables and levels of these factors are
  equal to zero.

  The likelihood ratio test p values therefore represent a test of all the variables and all the levels of factors
  which are among these variables. However, the results table only has space for one column of log fold change, so
  a single variable and a single comparison is shown (among the potentially multiple log fold changes which were
  tested in the likelihood ratio test). This indicates that the p value is for the likelihood ratio test of all the
  variables and all the levels, while the log fold change is a single comparison from among those variables and levels.

  Note: at least two biological replicates are required for every compared category.

  All input CSV/TSV files should have the following header (case-sensitive)
  <RefseqId,GeneId,Chrom,TxStart,TxEnd,Strand,TotalReads,Rpkm>         - CSV
  <RefseqId\tGeneId\tChrom\tTxStart\tTxEnd\tStrand\tTotalReads\tRpkm>  - TSV

  Format of the input files is identified based on file's extension
  *.csv - CSV
  *.tsv - TSV
  Otherwise used CSV by default

  The output file's rows order corresponds to the rows order of the first CSV/TSV file.
  Output file is always saved in TSV format

  Output file includes only intersected rows from all input files. Intersected by
  RefseqId, GeneId, Chrom, TxStart, TxEnd, Strand

  Additionally we calculate -LOG10(pval) and -LOG10(padj)

  Example of CSV metadata file set with --meta

  ,time,condition
  DH1,day5,WT
  DH2,day5,KO
  DH3,day7,WT
  DH4,day7,KO
  DH5,day7,KO

  where time, condition, day5, day7, WT, KO should be a single words (without spaces)
  and DH1, DH2, DH3, DH4, DH5 correspond to the --names (spaces are allowed)

  --contrast should be set based on your metadata file in a form of Factor Numerator Denominator
  where Factor      - columns name from metadata file
        Numerator   - category from metadata file to be used as numerator in fold change calculation
        Denominator - category from metadata file to be used as denominator in fold change calculation
  for example condition WT KO
  if --contrast is set as a single string "condition WT KO" then is will be splitted by space


s:about: |
  usage: run_deseq_lrt.R
        [-h] -i INPUT [INPUT ...] -n NAME [NAME ...] -m META -d DESIGN -r
        REDUCED -c CONTRAST [CONTRAST ...] [-o OUTPUT] [-p THREADS]

  Run DeSeq2 for multi-factor analysis using LRT (likelihood ratio or chi-
  squared test)

  optional arguments:
    -h, --help            show this help message and exit
    -i INPUT [INPUT ...], --input INPUT [INPUT ...]
                          Grouped by gene / TSS/ isoform expression files,
                          formatted as CSV/TSV
    -n NAME [NAME ...], --name NAME [NAME ...]
                          Unique names for input files, no special characters,
                          spaces are allowed. Number and order corresponds to
                          --input
    -m META, --meta META  Metadata file to describe relation between samples,
                          where first column corresponds to --name, formatted as
                          CSV/TSV
    -d DESIGN, --design DESIGN
                          Design formula. Should start with ~. See DeSeq2 manual
                          for details
    -r REDUCED, --reduced REDUCED
                          Reduced formula to compare against with the term(s) of
                          interest removed. Should start with ~. See DeSeq2
                          manual for details
    -c CONTRAST [CONTRAST ...], --contrast CONTRAST [CONTRAST ...]
                          Contrast to be be applied for output, formatted as
                          Factor Numerator Denominator or "Factor Numerator
                          Denominator"
    -o OUTPUT, --output OUTPUT
                          Output prefix for generated files
    -p THREADS, --threads THREADS
                          Threads number
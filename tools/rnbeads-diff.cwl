cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement


hints:
- class: DockerRequirement
  dockerPull: robertplayer/scidap-rnbeads:v1.0.0


inputs:

  threads:
    type: int
    inputBinding:
      prefix: "-t"
    doc: "Number of threads for parallel processing"

  genome:
    type: string
    inputBinding:
      prefix: "-g"
    doc: "Sample genome, available options: hg19, hg38, mm9, mm10, rn5"

  condition1_name:
    type: string
    inputBinding:
      prefix: "-a"
    doc: "Name for group 1, default: condition1"

  condition2_name:
    type: string
    inputBinding:
      prefix: "-b"
    doc: "Name for group 2, default: condition2"

  condition1_filepaths:
    type: File[]
    inputBinding:
      prefix: "-c"
      itemSeparator: ","
    doc: "List of absolute filepaths to BismarkCov formatted bed files for group 1"

  condition2_filepaths:
    type: File[]
    inputBinding:
      prefix: "-d"
      itemSeparator: ","
    doc: "List of absolute filepaths to BismarkCov formatted bed files for group 2"

  condition1_aliases:
    type: string[]
    inputBinding:
      prefix: "-j"
      itemSeparator: ","
    doc: "List of alias sample names for group 1"

  condition2_aliases:
    type: string[]
    inputBinding:
      prefix: "-k"
      itemSeparator: ","
    doc: "List of alias sample names for group 2"

  refgene_annotations:
    type: File
    inputBinding:
      prefix: "-m"
    doc: "refgene annotations text file"

  condition1_methperc_bigwigs:
    type: File[]
    inputBinding:
      prefix: "-p"
      itemSeparator: ","
    doc: "List of absolute filepaths to methylation percent bigWig files for group 1"

  condition2_methperc_bigwigs:
    type: File[]
    inputBinding:
      prefix: "-q"
      itemSeparator: ","
    doc: "List of absolute filepaths to methylation percent bigWig files for group 2"


outputs:

  c1_bigwigs:
    type: File[]
    outputBinding:
      glob: "c1_*.bigWig"
    doc: |
      Comma separated list of absolute filepaths to all condition1 bigWig files (from Bismark upstream).

  c2_bigwigs:
    type: File[]
    outputBinding:
      glob: "c2_*.bigWig"
    doc: |
      Comma separated list of absolute filepaths to all condition2 bigWig files (from Bismark upstream).

  samplesheet:
    type: File
    outputBinding:
      glob: sample_annotation.csv
    doc: |
      Prepared samplesheet from user inputs, used for condition labeling by RnBeads.

  samplesheet_overview:
    type: File
    outputBinding:
      glob: sample_annotation.md
    doc: |
      Prepared samplesheet from user inputs, used for display on Overview tab.

  report_tar:
    type: File
    outputBinding:
      glob: reports.tar.gz
    doc: |
      Compressed reports dir tarball containing html summary files and additional files.

  report_directory:
    type: Directory
    outputBinding:
      glob: reports
    doc: |
      output directory for all rnbeads results, preserves structure for html references

  report_data_import_html:
    type: File
    outputBinding:
      glob: reports/data_import.html
    doc: |
      Import data report

  report_qc_html:
    type: File
    outputBinding:
      glob: reports/quality_control.html
    doc: |
      Data QC report

  report_preprocessing_html:
    type: File
    outputBinding:
      glob: reports/preprocessing.html
    doc: |
      preprocessing report

  report_differential_methylation_html:
    type: File
    outputBinding:
      glob: reports/differential_methylation.html
    doc: |
      Differential methylation HTML report for visualization.

  meth_stats_genes:
    type: File
    outputBinding:
      glob: dm_genes.tsv
    doc: |
      Stats for all differentially methylated Genes: quot.log2,p.adj.fdr,num.sites,mean.covg.condition1,mean.covg.condition2

  meth_stats_promoters:
    type: File
    outputBinding:
      glob: dm_promoters.tsv
    doc: |
      Stats for all differentially methylated Promoters: quot.log2,p.adj.fdr,num.sites,mean.covg.condition1,mean.covg.condition2

  meth_stats_cpg:
    type: File
    outputBinding:
      glob: dm_cpg.tsv
    doc: |
      Stats for all differentially methylated CpG islands: quot.log2,comb.p.adj.fdr,num.sites,mean.covg.condition1,mean.covg.condition2

  meth_stats_tiling:
    type: File
    outputBinding:
      glob: dm_tiling.tsv
    doc: |
      Stats for all differentially methylated Tilings of 5kbp: quot.log2,comb.p.adj.fdr,num.sites,mean.covg.condition1,mean.covg.condition2

  meth_stats_sites:
    type: File
    outputBinding:
      glob: dm_sites.tsv
    doc: |
      Stats for all differentially methylated Sites: quot.log2,diffmeth.p.adj.fdr,mean.covg.condition1,mean.covg.condition2

  sig_dm_sites_igvtrack:
    type: File
    outputBinding:
      glob: sig_dm_sites.bed
    doc: |
      Bed file with locations for significantly (FDR<0.10) differentially methylated Sites with log2 fold change and FDR values

  sig_dm_sites_annotated:
    type: File
    outputBinding:
      glob: sig_dm_sites_annotated.tsv
    doc: |
      Stats for significantly (FDR<0.10) differentially methylated Sites with single closest gene annotation

  sig_dm_cpg_igvtrack:
    type: File
    outputBinding:
      glob: sig_dm_cpg.bed
    doc: |
      Bed file with locationss for significantly (FDR<0.10) differentially methylated CpG Islands with log2 fold change and FDR values

  sig_dm_cpg_annotated:
    type: File
    outputBinding:
      glob: sig_dm_cpg_annotated.tsv
    doc: |
      Stats for significantly (FDR<0.10) differentially methylated CpG Islands with single closest gene annotation

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["run_rnbeads_diff.sh"]
stdout: rnbeads_diff_stdout.log
stderr: rnbeads_diff_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "rnbeads-diff"
s:downloadUrl: https://github.com/datirium/workflows/blob/master/tools/rnbeads-diff.cwl
s:codeRepository: https://github.com/datirium/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Datirium LLC"
  s:location:
  - class: s:PostalAddress
    s:addressCountry: "USA"
    s:addressLocality: "Cincinnati"
    s:addressRegion: "OH"
    s:postalCode: ""
    s:streetAddress: ""
    s:telephone: ""
  s:logo: "https://avatars.githubusercontent.com/u/33202955?s=200&v=4"
  s:department:
  - class: s:Organization
    s:legalName: "Datirium LLC"
    s:department:
    - class: s:Organization
      s:legalName: "Bioinformatics"
      s:member:
      - class: s:Person
        s:name: Robert Player
        s:email: mailto:support@datirium.com
        s:sameAs:
        - id: https://orcid.org/0000-0001-5872-259X


doc: |
  Wrapper for RnBeads differential methylation pipeline script, with downstream processing for tables and IGV

  Primary Output files:
    - sig_dm_sites.bed (bed for IGV; sig diff meth sites)
    - sig_dm_sites_annotated.tsv (tsv for TABLE; for each site above, closest single gene annotation)

  Output rnbeads reports directory in container at '/tmp/reports/', includes:
    reports/
  ├── configuration
  ├── data_import.html
  ├── data_import_data
  ├── data_import_images
  ├── data_import_pdfs
  ├── differential_methylation.html
  ├── differential_methylation_data
  ├── differential_methylation_images
  ├── differential_methylation_pdfs
  ├── preprocessing.html
  ├── preprocessing_data
  ├── preprocessing_images
  ├── preprocessing_pdfs
  ├── quality_control.html
  ├── quality_control_data
  ├── quality_control_images
  ├── quality_control_pdfs
  ├── tracks_and_tables.html
  ├── tracks_and_tables_data
  ├── tracks_and_tables_images
  └── tracks_and_tables_pdfs

  Reported methylation is in the form of regions (genes, promoters, cpg, tiling) and specific sites:
   - genes - Ensembl gene definitions are downloaded using the biomaRt package.
   - promoters - A promoter is defined as the region spanning 1,500 bases upstream and 500 bases downstream of the transcription start site of the corresponding gene
   - cpg - the CpG islands from the UCSC Genome Browser
   - tiling - a window size of 5 kilobases are defined over the whole genome
   - sites - all cytosines in the context of CpGs in the respective genome

  SCRIPT PARAMS:
  -h  help	show this message
  -g  STRING   Sample genome, available options: hg19, hg38, mm9, mm10, rn5
  -t  INT	number of threads
  -a  STRING     name of condition1
  -b  STRING     name of condition2
  -c  LIST	comma separated list of absolute filepaths to all condition1 bed files (BismarkCov format)
  -d  LIST	comma separated list of absolute filepaths to all condition2 bed files (BismarkCov format)
  -j  LIST   comma separated list of sample names in condition1
  -k  LIST   comma separated list of sample names in condition2
  -m  FILE   refGene.txt file for annotating DM sites with gene information

  BismarkCov formatted bed:
      https://www.bioinformatics.babraham.ac.uk/projects/bismark/Bismark_User_Guide.pdf
      The genome-wide cytosine report (optional) is tab-delimited in the following format (1-based coords):
      <chromosome> <position> <strand> <count methylated> <count unmethylated> <C-context> <trinucleotide context>

  ____________________________________________________________________________________________________
  References:
      https://rnbeads.org/materials/example_3/differential_methylation.html
          Makambi, K. (2003) Weighted inverse chi-square method for correlated significance tests.
          Journal of Applied Statistics, 30(2), 225234
      https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4216143/
          Assenov Y, Müller F, Lutsik P, Walter J, Lengauer T, Bock C. Comprehensive analysis of DNA
          methylation data with RnBeads. Nat Methods. 2014 Nov;11(11):1138-1140. doi: 10.1038/nmeth.3115.
          Epub 2014 Sep 28. PMID: 25262207; PMCID: PMC4216143.

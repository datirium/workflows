cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/scidap-deseq:v0.0.15


inputs:

  untreated_files:
    type:
      - File
      - File[]
    inputBinding:
      position: 5
      prefix: "-u"
    doc: |
      Untreated input CSV/TSV files

  treated_files:
    type:
      - File
      - File[]
    inputBinding:
      position: 6
      prefix: "-t"
    doc: |
      Treated input CSV/TSV files

  untreated_name:
    type: string?
    inputBinding:
      position: 7
      prefix: "-un"
    doc: |
      Name for untreated condition, use only letters and numbers

  treated_name:
    type: string?
    inputBinding:
      position: 8
      prefix: "-tn"
    doc: |
      Name for treated condition, use only letters and numbers

  untreated_sample_names:
    type:
      - "null"
      - string
      - string[]
    inputBinding:
      prefix: "-ua"
    doc: |
      Unique aliases for untreated expression files. Default: basenames of -u without extensions

  treated_sample_names:
    type:
      - "null"
      - string
      - string[]
    inputBinding:
      prefix: "-ta"
    doc: |
      Unique aliases for treated expression files. Default: basenames of -t without extensions

  output_prefix:
    type: string?
    inputBinding:
      position: 9
      prefix: "-o"
    doc: |
      Output prefix. Default: deseq

  threads:
    type: int?
    inputBinding:
      position: 10
      prefix: '-p'
    doc: |
      Run script using multiple threads


outputs:

  diff_expr_file:
    type: File
    outputBinding:
      glob: "*report.tsv"

  read_counts_file:
    type: File
    outputBinding:
      glob: "*counts.gct"

  phenotypes_file:
    type: File
    outputBinding:
      glob: "*phenotypes.cls"

  plot_lfc_vs_mean:
    type: File?
    outputBinding:
      glob: "*_ma_plot.png"

  gene_expr_heatmap:
    type: File?
    outputBinding:
      glob: "*_expression_heatmap.png"

  plot_pca:
    type: File?
    outputBinding:
      glob: "*_pca_plot.png"

  plot_lfc_vs_mean_pdf:
    type: File?
    outputBinding:
      glob: "*_ma_plot.pdf"

  gene_expr_heatmap_pdf:
    type: File?
    outputBinding:
      glob: "*_expression_heatmap.pdf"

  plot_pca_pdf:
    type: File?
    outputBinding:
      glob: "*_pca_plot.pdf"

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: [run_deseq.R]
stdout: deseq_stdout.log
stderr: deseq_stderr.log

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: ./metadata/deseq-metadata.yaml

s:name: "deseq-advanced"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/deseq-advanced.cwl
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
  Tool runs DESeq/DESeq2 script similar to the original one from BioWArdrobe.
  untreated_files and treated_files input files should have the following header (case-sensitive)
  <RefseqId,GeneId,Chrom,TxStart,TxEnd,Strand,TotalReads,Rpkm>         - CSV
  <RefseqId\tGeneId\tChrom\tTxStart\tTxEnd\tStrand\tTotalReads\tRpkm>  - TSV

  Format of the input files is identified based on file's extension
  *.csv - CSV
  *.tsv - TSV
  Otherwise used CSV by default

  The output file's rows order corresponds to the rows order of the first CSV/TSV file in
  the untreated group. Output is always saved in TSV format

  Output file includes only intersected rows from all input files. Intersected by
  RefseqId, GeneId, Chrom, TxStart, TxEnd, Strand

  DESeq/DESeq2 always compares untreated_vs_treated groups.
  Normalized read counts and phenotype table are exported as GCT and CLS files for GSEA downstream analysis.


s:about: |
  usage: run_deseq.R
        [-h] -u UNTREATED [UNTREATED ...] -t TREATED [TREATED ...]
        [-ua [UALIAS [UALIAS ...]]] [-ta [TALIAS [TALIAS ...]]] [-un UNAME]
        [-tn TNAME] [-o OUTPUT] [-p THREADS]

  Run BioWardrobe DESeq/DESeq2 for untreated-vs-treated groups

  optional arguments:
    -h, --help            show this help message and exit
    -u UNTREATED [UNTREATED ...], --untreated UNTREATED [UNTREATED ...]
                          Untreated CSV/TSV isoforms expression files
    -t TREATED [TREATED ...], --treated TREATED [TREATED ...]
                          Treated CSV/TSV isoforms expression files
    -ua [UALIAS [UALIAS ...]], --ualias [UALIAS [UALIAS ...]]
                          Unique aliases for untreated expression files.
                          Default: basenames of -u without extensions
    -ta [TALIAS [TALIAS ...]], --talias [TALIAS [TALIAS ...]]
                          Unique aliases for treated expression files. Default:
                          basenames of -t without extensions
    -un UNAME, --uname UNAME
                          Name for untreated condition, use only letters and
                          numbers
    -tn TNAME, --tname TNAME
                          Name for treated condition, use only letters and
                          numbers
    -o OUTPUT, --output OUTPUT
                          Output prefix. Default: deseq
    -p THREADS, --threads THREADS
                          Threads

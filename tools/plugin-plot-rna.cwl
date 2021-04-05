cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/plugin-plot-rna:v0.0.3


inputs:

  annotation_file:
    type: File
    inputBinding:
      position: 5
      prefix: "--annotation"
    doc: |
      Path to the annotation TSV/CSV file

  bambai_pair:
    type: File
    inputBinding:
      position: 6
      prefix: "--bam"
    secondaryFiles:
    - .bai
    doc: |
      Path to the indexed BAM file

  isoforms_file:
    type: File
    inputBinding:
      position: 7
      prefix: "--isoforms"
    doc: |
      Path to the isoforms TSV/CSV file

  mapped_reads_number:
    type: int
    inputBinding:
      position: 8
      prefix: "--mapped"
    doc: |
      Mapped reads number

  output_prefix:
    type: string?
    inputBinding:
      position: 9
      prefix: "--output"
    doc: |
      Output prefix. Default: ./coverage

  pair:
    type: boolean?
    inputBinding:
      position: 10
      prefix: "--pair"
    doc: |
      Run as paired end. Default: false

  dutp:
    type: boolean?
    inputBinding:
      position: 11
      prefix: "--dutp"
    doc: |
      Run as dUTP. Default: false

  minimum_rpkm:
    type: float?
    inputBinding:
      position: 12
      prefix: "--minrpkm"
    doc: |
      Ignore isoforms with RPKM smaller than --minrpkm.
      Default: 10

  minimum_isoform_length:
    type: float?
    inputBinding:
      position: 13
      prefix: "--minlength"
    doc: |
      Ignore isoforms shorter than --minlength.
      Default: 1000

  threads:
    type: int?
    inputBinding:
      position: 14
      prefix: "--threads"
    doc: |
      Threads. Default: 1


outputs:

  gene_body_report_file:
    type: File
    outputBinding:
      glob: "*gene_body_report.tsv"

  gene_body_plot_png:
    type: File
    outputBinding:
      glob: "*gene_body_plot.png"

  gene_body_plot_pdf:
    type: File
    outputBinding:
      glob: "*gene_body_plot.pdf"

  rpkm_distribution_plot_png:
    type: File
    outputBinding:
      glob: "*rpkm_distribution_plot.png"

  rpkm_distribution_plot_pdf:
    type: File
    outputBinding:
      glob: "*rpkm_distribution_plot.pdf"

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["plot_rna.R"]


stderr: gene_body_stderr.log
stdout: gene_body_stdout.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Gene body average tag density plot and RPKM distribution histogram"
s:name: "Gene body average tag density plot and RPKM distribution histogram"
s:alternateName: "Gene body average tag density plot and RPKM distribution histogram"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/plugin-plot-rna.cwl
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
  Runs R script to produce gene body average tag density plot and RPKM distribution histogram

s:about: |
  usage: plugin_plot_rna.R
        [-h] --annotation ANNOTATION --bam BAM --isoforms ISOFORMS
        [--minrpkm MINRPKM] [--minlength MINLENGTH] --mapped MAPPED [--pair]
        [--dutp] [--output OUTPUT] [--threads THREADS]

  Gene body average tag density plot and RPKM distribution histogram for
  isoforms

  optional arguments:
    -h, --help            show this help message and exit
    --annotation ANNOTATION
                          Path to the annotation TSV/CSV file
    --bam BAM             Path to the indexed BAM file
    --isoforms ISOFORMS   Path to the isoforms TSV/CSV file
    --minrpkm MINRPKM     Ignore isoforms with RPKM smaller than --minrpkm.
                          Default: 10
    --minlength MINLENGTH
                          Ignore isoforms shorter than --minlength. Default:
                          1000
    --mapped MAPPED       Mapped reads number
    --pair                Run as paired end. Default: false
    --dutp                Run as dUTP. Default: false
    --output OUTPUT       Output prefix. Default: ./coverage
    --threads THREADS     Threads. Default: 1
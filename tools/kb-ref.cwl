cwlVersion: v1.0
class: CommandLineTool


requirements:
  - class: InlineJavascriptRequirement
    expressionLib:
    - var get_output_prefix = function(ext) {
            if (inputs.output_prefix == ""){
              var root = inputs.genome_fasta_file.basename.split('.').slice(0,-1).join('.');
              return (root == "")?inputs.genome_fasta_file.basename+ext:root+ext;
            } else {
              return inputs.output_prefix + ext;
            }
          };


hints:
  - class: DockerRequirement
    dockerPull: biowardrobe2/kb-python:v0.0.2


inputs:

  genome_fasta_file:
    type: File
    inputBinding:
      position: 100
    doc: |
      Genome FASTA file

  annotation_gtf_file:
    type: File
    inputBinding:
      position: 101
    doc: |
      GTF annotation file

  feature_tsv_file:
    type: File?
    inputBinding:
      position: 102
    doc: |
      TSV file containing barcodes and feature names.
      Should be provided only when workflow_type is kite

  workflow_type:
    type:
    - "null"
    - type: enum
      name: "workflow_type"
      symbols:
      - standard
      - lamanno
      - nucleus
      - kite
    inputBinding:
      position: 20
      prefix: "--workflow"
    doc: |
      Type of workflow. Use lamanno to calculate RNA velocity based
      on La Manno et al. 2018 logic. Use nucleus to calculate RNA
      velocity on single-nucleus RNA-seq reads.
      Default: standard

  output_prefix:
    type: string?
    default: ""
    doc: |
      Output prefix for generated files


outputs:

  kallisto_index_file:
    type: File
    outputBinding:
      glob: $(get_output_prefix("_index.idx"))
    doc: |
      Kallisto index file

  tx_to_gene_mapping_file:
    type: File
    outputBinding:
      glob: $(get_output_prefix("_tx_to_gene.tsv"))      
    doc: |
      Transcript-to-gene mapping TSV file

  tx_fasta_file:
    type: File
    outputBinding:
      glob: $(get_output_prefix("_tx.fasta"))
    doc: |
      Transcrip FASTA (for lamanno or nucleus workflow_type)
      or mismatch FASTA (for kite workflow_type) file

  intron_fasta_file:
    type: File?
    outputBinding:
      glob: $(get_output_prefix("_intron.fasta"))
    doc: |
      Intron FASTA file

  tx_to_capture_mapping_file:
    type: File?
    outputBinding:
      glob: $(get_output_prefix("_tx_to_capture.tsv"))
    doc: |
      Transcripts-to-capture mapping TSV file

  intron_tx_to_capture_mapping_file:
    type: File?
    outputBinding:
      glob: $(get_output_prefix("_intron_tx_to_capture.tsv"))
    doc: |
      Intron transcripts-to-capture mapping TSV file

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["kb", "ref", "--verbose"]


arguments:
- valueFrom: $(get_output_prefix("_index.idx"))
  position: 5
  prefix: -i
- valueFrom: $(get_output_prefix("_tx_to_gene.tsv"))
  position: 6
  prefix: -g
- valueFrom: $(get_output_prefix("_tx.fasta"))
  position: 7
  prefix: -f1

# arguments used only when workflow_type is lamanno or nucleus but 
# we still append them to basecommand to make cwl less complicate

- valueFrom: $(get_output_prefix("_intron.fasta"))
  position: 8
  prefix: -f2
- valueFrom: $(get_output_prefix("_tx_to_capture.tsv"))
  position: 9
  prefix: -c1
- valueFrom: $(get_output_prefix("_intron_tx_to_capture.tsv"))
  position: 10
  prefix: -c2


stdout: kb_ref_stdout.log
stderr: kb_ref_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "kb-ref"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/kb-ref.cwl
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
  Builds a kallisto index, transcript-to-gene mapping and cDNA/mismatch FASTA  files.
  If workflow_type is lamanno or nucleus tool will produce additional three
  outputs:
  - intron FASTA file
  - cDNA transcripts-to-capture TSV file
  - intron transcripts-to-capture TSV file
  Otherwise the correspondent outputs will be null.

  Notes:
  --verbose was hardcoded
  --keep-tmp, -d, --overwrite doesn't make sense when running from container
  -f2, -c1 and -c2 are always appended to the basecommand regardless of workflow_type (makes cwl less complicate)

  `annotation_gtf_file` input should have correct "gene_id" field.
  
  To generate correct GTF from refgene annotations use:
  `docker run --rm -ti -v `pwd`:/tmp/ biowardrobe2/ucscuserapps:v358 /bin/bash -c "cut -f 2- refGene.txt | genePredToGtf file stdin refgene.gtf"`
  to generate a proper gtf file from `refGene.txt` downloaded from http://hgdownload.cse.ucsc.edu/goldenPath/${GEN}/database/refGene.txt.gz

s:about: |
  usage: kb ref [-h] [--workflow {standard,lamanno,nucleus,kite}] [--keep-tmp] [--verbose] -i INDEX -g T2G -f1 FASTA [-f2 FASTA] [-c1 T2C] [-c2 T2C] [-d {human,mouse,linnarsson}]
                [--lamanno] [--overwrite]
                fasta gtf [feature]

  Build a kallisto index and transcript-to-gene mapping

  positional arguments:
    fasta                 Genomic FASTA file
    gtf                   Reference GTF file
    feature               [`kite` workflow only] Path to TSV containing barcodes and feature names.

  optional arguments:
    -h, --help            Show this help message and exit
    --workflow {standard,lamanno,nucleus,kite}
                          Type of workflow. Use `lamanno` to calculate RNA velocity based on La Manno et al. 2018 logic. Use `nucleus` to calculate RNA velocity on single-nucleus
                          RNA-seq reads (default: standard)
    --keep-tmp            Do not delete the tmp directory
    --verbose             Print debugging information
    -d {human,mouse,linnarsson}
                          Download a pre-built kallisto index (along with all necessary files) instead of building it locally
    --lamanno             Deprecated. Use `--workflow lamanno` instead.
    --overwrite           Overwrite existing kallisto index

  required arguments:
    -i INDEX              Path to the kallisto index to be constructed
    -g T2G                Path to transcript-to-gene mapping to be generated
    -f1 FASTA             [Optional with -d] Path to the cDNA FASTA (lamanno, nucleus) or mismatch FASTA (kite) to be generated

  required arguments for `lamanno` and `nucleus` workflows:
    -f2 FASTA             Path to the intron FASTA to be generated
    -c1 T2C               Path to generate cDNA transcripts-to-capture
    -c2 T2C               Path to generate intron transcripts-to-capture
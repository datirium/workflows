cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/samtools:v1.4
  dockerFile: >
    $import: ./dockerfiles/samtools-Dockerfile

inputs:

  fastcompression:
    type: boolean?
    inputBinding:
      position: 6
      prefix: '-1'
    doc: |
      Compress level 1

  compression_level:
    type: int?
    inputBinding:
      position: 7
      prefix: -l
    doc: |
      Compression level, from 0 to 9 [-1]

  uncompressed:
    type: boolean?
    inputBinding:
      position: 8
      prefix: -u
    doc: |
      uncompressed BAM output

  samheader:
    type: File?
    inputBinding:
      position: 9
      prefix: -h
    doc: |
      Copy the header in FILE to <out.bam> [in1.bam]

  reference_fasta:
    type: File?
    inputBinding:
      position: 10
      prefix: --reference
    doc: |
      Reference sequence FASTA FILE [null]

  attach_rg_tag:
    type: string?
    inputBinding:
      position: 11
      prefix: -r
    doc: |
      Attach RG tag (inferred from file names)

  threads:
    type: int?
    inputBinding:
      position: 12
      prefix: -@
    doc: |
      number of BAM compression threads [0]

  output_name:
    type: string
    inputBinding:
      position: 99
    default: "merged.bam"
    doc: |
      Output to file

  input:
    type:
      - File
      - type: array
        items: File
    inputBinding:
      position: 100
    doc: |
      Input SAM, BAM, or CRAM file.

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.output_name)

baseCommand: [samtools, merge]

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: ./metadata/samtools-metadata.yaml

s:downloadUrl: https://github.com/SciDAP/workflows/blob/master/tools/samtools-view.cwl
s:codeRepository: https://github.com/SciDAP/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:author:
  class: s:Person
  s:name: Andrey Kartashov
  s:email: mailto:Andrey.Kartashov@cchmc.org
  s:sameAs:
  - id: http://orcid.org/0000-0001-9102-5681
  s:worksFor:
  - class: s:Organization
    s:name: Cincinnati Children's Hospital Medical Center
    s:location: 3333 Burnet Ave, Cincinnati, OH 45229-3026
    s:department:
    - class: s:Organization
      s:name: Barski Lab
doc: |
  samtools-merge.cwl is developed for CWL consortium
  Usage: samtools merge [-nurlf] [-h inh.sam] [-b <bamlist.fofn>] <out.bam> <in1.bam> [<in2.bam> ... <inN.bam>]

  Options:
    -n         Input files are sorted by read name
    -r         Attach RG tag (inferred from file names)
    -u         Uncompressed BAM output
    -f         Overwrite the output BAM if exist
    -1         Compress level 1
    -l INT     Compression level, from 0 to 9 [-1]
    -R STR     Merge file in the specified region STR [all]
    -h FILE    Copy the header in FILE to <out.bam> [in1.bam]
    -c         Combine @RG headers with colliding IDs [alter IDs to be distinct]
    -p         Combine @PG headers with colliding IDs [alter IDs to be distinct]
    -s VALUE   Override random seed
    -b FILE    List of input BAM filenames, one per line [null]
    -@, --threads INT
               Number of BAM/CRAM compression threads [0]
        --input-fmt-option OPT[=VAL]
                 Specify a single input file format option in the form
                 of OPTION or OPTION=VALUE
    -O, --output-fmt FORMAT[,OPT[=VAL]]...
                 Specify output format (SAM, BAM, CRAM)
        --output-fmt-option OPT[=VAL]
                 Specify a single output file format option in the form
                 of OPTION or OPTION=VALUE
        --reference FILE
                 Reference sequence FASTA FILE [null]


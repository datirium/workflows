cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function(ext) {
        var ext = inputs.alignment_files[0].basename.split('.').slice(-1)[0];
        var root = inputs.alignment_files[0].basename.split('.').slice(0,-1).join('.');
        return inputs.output_filename?inputs.output_filename:root+"_merged."+ext;
    };


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/samtools:v1.4


inputs:

  fastcompression:
    type: boolean?
    inputBinding:
      position: 6
      prefix: "-1"
    doc: "Compress level 1"

  compression_level:
    type: int?
    inputBinding:
      position: 7
      prefix: "-l"
    doc: "Compression level, from 0 to 9 [-1]"

  uncompressed:
    type: boolean?
    inputBinding:
      position: 8
      prefix: "-u"
    doc: "uncompressed BAM output"

  samheader:
    type: File?
    inputBinding:
      position: 9
      prefix: "-h"
    doc: "Copy the header in FILE to <out.bam> [in1.bam]"

  reference_fasta:
    type: File?
    inputBinding:
      position: 10
      prefix: "--reference"
    doc: "Reference sequence FASTA FILE [null]"

  attach_rg_tag:
    type: string?
    inputBinding:
      position: 11
      prefix: "-r"
    doc: "Attach RG tag (inferred from file names)"

  threads:
    type: int?
    inputBinding:
      position: 12
      prefix: "-@"
    doc: "Number of BAM compression threads [0]"

  output_filename:
    type: string?
    inputBinding:
      position: 99
      valueFrom: $(default_output_filename())
    default: ""
    doc: "Output filename"

  alignment_files:
    type: File[]
    inputBinding:
      position: 100
    doc: "Input SAM, BAM, or CRAM files"


outputs:

  merged_alignment_file:
    type: File
    outputBinding:
      glob: $(default_output_filename())

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: [samtools, merge]
stdout: samtools_merge_stdout.log
stderr: samtools_merge_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: ./metadata/samtools-metadata.yaml

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
  If output_filename is not provided, the default output name for merged file is generated as follows:
  name without extestion of alignment_files[0] + suffix "_merged" + extension of alignment_files[0]

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


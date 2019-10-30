cwlVersion: v1.0
class: CommandLineTool


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/ucscuserapps:v358


requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() {
          if (inputs.output_filename == ""){
            var root = inputs.twobit_file.basename.split('.').slice(0,-1).join('.');
            return (root == "")?inputs.twobit_file.basename+".fa":root+".fa";
          } else {
            return inputs.output_filename;
          }
        };


inputs:

  script:
    type: string?
    default: |
      #!/bin/bash
      if [ "$#" -gt 2 ]; then
        echo ${@:2} | tr " " "\n" > selected_chr.tsv
        twoBitToFa -seqList=selected_chr.tsv $0 $1 
      else
        twoBitToFa $0 $1
      fi
    inputBinding:
      position: 5
    doc: |
      Bash function to run twoBitToFa with optional chromosome filtering

  twobit_file:
    type: File
    inputBinding:
      position: 6
    doc: "Reference genome 2bit file"

  output_filename:
    type: string?
    default: ""
    inputBinding:
      position: 7
      valueFrom: $(default_output_filename())
    doc: "Output file name"

  chr_list:
    type:
      - "null"
      - string[]
    inputBinding:
      position: 8
    doc: "List of the chromosomes to be included into the output file"

outputs:

  fasta_file:
    type: File
    outputBinding:
      glob: $(default_output_filename())
    doc: "Reference genome FASTA file"

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["bash", "-c"]
stdout: twobit_to_fa_stdout.log
stderr: twobit_to_fa_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: ./metadata/ucsc-metadata.yaml

s:name: "ucsc-twobit-to-fa"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/ucsc-twobit-to-fa.cwl
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
  twoBitToFa - Convert all or part of .2bit file to fasta.
  Outputs only those chromosomes that are set in chr_list intput.
  Tool will fail if you include in chr_list those chromosomes that are absent in 2bit file.

s:about: |
  usage:
    twoBitToFa input.2bit output.fa
  options:
    -seq=name       Restrict this to just one sequence.
    -start=X        Start at given position in sequence (zero-based).
    -end=X          End at given position in sequence (non-inclusive).
    -seqList=file   File containing list of the desired sequence names 
                    in the format seqSpec[:start-end], e.g. chr1 or chr1:0-189
                    where coordinates are half-open zero-based, i.e. [start,end).
    -noMask         Convert sequence to all upper case.
    -bpt=index.bpt  Use bpt index instead of built-in one.
    -bed=input.bed  Grab sequences specified by input.bed. Will exclude introns.
    -bedPos         With -bed, use chrom:start-end as the fasta ID in output.fa.
    -udcDir=/dir/to/cache  Place to put cache for remote bigBed/bigWigs.

  Sequence and range may also be specified as part of the input
  file name using the syntax:
        /path/input.2bit:name
    or
        /path/input.2bit:name
    or
        /path/input.2bit:name:start-end
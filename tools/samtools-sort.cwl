#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: ./metadata/envvar-global.yml
- class: ShellCommandRequirement
- class: InitialWorkDirRequirement
  listing: |
    ${
      return  [
                {
                  "entry": inputs.sort_input,
                  "entryname": inputs.sort_input.basename,
                  "writable": true
                }
              ]
    }
- class: InlineJavascriptRequirement
  expressionLib:
  - var ext = function() {
      if (inputs.out_format && inputs.out_format == 'SAM'){
        return '.sam';
      } else if (inputs.out_format && inputs.out_format == 'CRAM') {
        return '.cram';
      } else {
        return '.bam';
      }
    };
  - var default_output_name = function() {
      if (inputs.trigger == true){
        return inputs.sort_input.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.') + ext();
      } else {
        return inputs.sort_input.location.split('/').slice(-1)[0];
      }
    };


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/samtools:v1.4
  dockerFile: >
    $import: ./dockerfiles/samtools-Dockerfile

inputs:

  bash_script_sort:
    type: string?
    default: |
      #!/bin/bash
      if [ "$0" = True ]
      then
        samtools sort "${@:1}"
      else
        echo "Skip samtools sort " ${@:1}
      fi
    inputBinding:
      position: 5
    doc: |
      Bash function to run samtools sort with all input parameters or skip it if trigger is false

  trigger:
    type: boolean?
    default: true
    inputBinding:
      position: 6
    doc: |
      If true - run samtools, if false - return sort_input, previously staged into output directory

  sort_compression_level:
    type: int?
    inputBinding:
      position: 7
      prefix: -l
    doc: |
      SORT: desired compression level for the final output file, ranging from 0 (uncompressed)
      or 1 (fastest but minimal compression) to 9 (best compression but slowest to write),
      similarly to gzip(1)'s compression level setting.
      If -l is not used, the default compression level will apply.

  sort_by_name:
    type: boolean?
    inputBinding:
      position: 8
      prefix: -n
    doc: |
      Sort by read names (i.e., the QNAME field) rather than by chromosomal coordinates

  sort_output_filename:
    type: string?
    inputBinding:
      position: 9
      prefix: -o
      valueFrom: |
        ${
            if (self == null || inputs.trigger == false){
              return default_output_name();
            } else {
              return self;
            }
        }
    default: null
    doc: |
      Write the final sorted output to FILE, rather than to standard output

  out_format:
    type: string?
    inputBinding:
      position: 10
      prefix: -O
    doc: |
      Specify output format (SAM, BAM, CRAM)

  threads:
    type: int?
    inputBinding:
      position: 11
      prefix: -@
    doc: |
      Set number of sorting and compression threads

  sort_input:
    type: File
    inputBinding:
      position: 12
    doc: |
      Input only in.sam|in.bam

outputs:
  sorted_file:
    type: File
    outputBinding:
      glob: |
        ${
            if (inputs.sort_output_filename == null || inputs.trigger == false){
              return default_output_name();
            } else {
              return inputs.sort_output_filename;
            }
        }

baseCommand: [bash, '-c']


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: ./metadata/samtools-metadata.yaml

s:name: "samtools-sort"
s:downloadUrl: https://raw.githubusercontent.com/SciDAP/workflows/master/tools/samtools-sort.cwl
s:codeRepository: https://github.com/SciDAP/workflows
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
        s:email: mailto:michael.kotliar@cchmc.org
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898

doc: |
  This tool is used to sort input BAM/SAM file by means of samtools sort
  Input Trigger (default: true) allows to skip all calculation and return
  all input files unchanged. To set files to be returned in case of Trigger == false,
  use the following inputs:
    sort_input

s:about: |
  Usage: samtools sort [options...] [in.bam]
  Options:
    -l INT     Set compression level, from 0 (uncompressed) to 9 (best)
    -m INT     Set maximum memory per thread; suffix K/M/G recognized [768M]
    -n         Sort by read name
    -o FILE    Write final output to FILE rather than standard output
    -T PREFIX  Write temporary files to PREFIX.nnnn.bam
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
    -@, --threads INT
                 Number of additional threads to use [0]

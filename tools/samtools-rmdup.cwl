#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: ./metadata/envvar-global.yml
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() {
      if (inputs.trigger == true){
        return inputs.input_file.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+".r.bam";
      } else {
        return inputs.input_file.location.split('/').slice(-1)[0];
      }
    };
- class: InitialWorkDirRequirement
  listing:
    - $(inputs.input_file)



hints:
- class: DockerRequirement
  dockerPull: scidap/samtools:v1.4
  dockerFile: >
    $import: ./dockerfiles/samtools-Dockerfile

inputs:

  bash_script:
    type: string?
    default: |
      #!/bin/bash
      if [ "$0" = True ]
      then
        echo "Run: samtools rmdup " ${@:1}
        samtools rmdup "${@:1}"
      else
        echo "Skip samtools rmdup " ${@:1}
      fi
    inputBinding:
      position: 1
    doc: |
      Bash function to run samtools rmdup with all input parameters or skip it if trigger is false

  trigger:
    type: boolean?
    default: true
    inputBinding:
      position: 2
    doc: |
      If true - run samtools rmdup, if false - return input_file, previously staged into output directory

  input_file:
    type: File
    inputBinding:
      position: 10
    doc: |
      Input sorted bam file.

  output_filename:
    type:
      - "null"
      - string
    inputBinding:
      position: 11
      valueFrom: |
        ${
            if (self == null || inputs.trigger == false){
              return default_output_filename();
            } else {
              return self;
            }
        }
    default: null
    doc: |
      Writes the output bam file to output_filename if set,
      otherwise generates output_filename on the base of input_file

  single_end:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 3
      prefix: '-s'
    doc: |
      rmdup for SE reads

  force_single_end:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 4
      prefix: '-S'
    doc: |
      treat PE reads as SE in rmdup (force -s)

outputs:
  rmdup_output:
    type: File
    outputBinding:
      glob: |
        ${
            if (inputs.output_filename == null || inputs.trigger == false){
              return default_output_filename();
            } else {
              return inputs.output_filename;
            }
        }
    secondaryFiles: |
      ${
          if (inputs.input_file.secondaryFiles && inputs.trigger == false){
            return inputs.input_file.secondaryFiles;
          } else {
            return "null";
          }
        }
    doc: File with removed duplicates or input_file with optional secondaryFiles

  rmdup_log:
    type: File
    outputBinding:
      glob: |
        ${
          if (inputs.output_filename == null || inputs.trigger == false){
            return default_output_filename() + '.rmdup';
          } else {
            return inputs.output_filename + '.rmdup';
          }
        }

baseCommand: [bash, '-c']

arguments:
  - valueFrom: |
      ${
        if (inputs.output_filename == null || inputs.trigger == false){
          return " > " + default_output_filename() + ".rmdup 2>&1";
        } else {
          return " > " + inputs.output_filename + ".rmdup 2>&1";
        }
      }
    position: 100000
    shellQuote: false


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: ./metadata/samtools-metadata.yaml

s:downloadUrl: https://github.com/SciDAP/workflows/blob/master/tools/samtools-rmdup.cwl
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
      - class: s:Person
        s:name: Andrey Kartashov
        s:email: mailto:Andrey.Kartashov@cchmc.org
        s:sameAs:
        - id: http://orcid.org/0000-0001-9102-5681

doc: |
  Usage:  samtools rmdup [-sS] <input.srt.bam> <output.bam>

  Option: -s    rmdup for SE reads
          -S    treat PE reads as SE in rmdup (force -s)
        --input-fmt-option OPT[=VAL]
                 Specify a single input file format option in the form
                 of OPTION or OPTION=VALUE
        --output-fmt FORMAT[,OPT[=VAL]]...
                 Specify output format (SAM, BAM, CRAM)
        --output-fmt-option OPT[=VAL]
                 Specify a single output file format option in the form
                 of OPTION or OPTION=VALUE
        --reference FILE
                 Reference sequence FASTA FILE [null]

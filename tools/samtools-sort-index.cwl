#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: ./metadata/envvar-global.yml
- class: ShellCommandRequirement
- class: InitialWorkDirRequirement
  listing:
    - $(inputs.sortInput)
- class: InlineJavascriptRequirement
  expressionLib:
  - var ext = function() {
      if (!inputs.sortOutputFileName){
        return '.bai';
      }
      if (inputs.sortOutputFileName.split('.').slice(-1)[0] == 'cram'){
        return '.crai';
      } else if (inputs.indexCsi && !inputs.indexBai){
        return '.csi';
      } else {
        return '.bai';
      }
    };
  - var default_bam = function() {
      return inputs.sortInput.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+".sorted.bam";
    };


hints:
- class: DockerRequirement
  dockerPull: scidap/samtools:v1.2-242-4d56437
  dockerFile: >
    $import: ./dockerfiles/samtools-Dockerfile

inputs:

  sortCompressionLevel:
    type: int?
    inputBinding:
      position: 2
      prefix: -l
    doc: |
      SORT: desired compression level for the final output file, ranging from 0 (uncompressed)
      or 1 (fastest but minimal compression) to 9 (best compression but slowest to write),
      similarly to gzip(1)'s compression level setting.
      If -l is not used, the default compression level will apply.

  sortByName:
    type: boolean?
    inputBinding:
      position: 4
      prefix: -n
    doc: |
      Sort by read names (i.e., the QNAME field) rather than by chromosomal coordinates

  sortOutputFileName:
    type: string?
    inputBinding:
      position: 3
      prefix: -o
      valueFrom: |
        ${
            if (self == null){
              return default_bam();
            } else {
              return self;
            }
        }
    default: null
    doc: |
      Write the final sorted output to FILE, rather than to standard output.
      Only out.bam|out.cram

  threads:
    type: int?
    inputBinding:
      position: 5
      prefix: -@
    doc: |
      Set number of sorting and compression threads [1] (Only for sorting)

  sortInput:
    type: File
    inputBinding:
      position: 6
      valueFrom: $(self.basename)
    doc: |
      Input only in.sam|in.bam|in.cram

  interval:
    type: int?
    inputBinding:
      position: 12
      prefix: -m
    doc: |
      Set minimum interval size for CSI indices to 2^INT [14]

  indexCsi:
    type: boolean?
    inputBinding:
      position: 11
      prefix: -c
      valueFrom: |
        ${
          if (ext() == '.crai' || ext() == '.bai'){
            return false;
          }
          return true;
        }
    doc: |
      Generate CSI-format index for BAM files. If input isn't cram.

  indexBai:
    type: boolean?
    inputBinding:
      position: 10
      prefix: -b
      valueFrom: |
        ${
          if (ext() == '.crai' || ext() == '.csi'){
            return false;
          }
          return true;
        }
    doc: |
      Generate BAI-format index for BAM files [default]. If input isn't cram.

outputs:
  bamBaiPair:
    type: File
    outputBinding:
      glob: |
        ${
            if (!inputs.sortOutputFileName){
              return default_bam();
            } else {
              return inputs.sortOutputFileName;
            }
        }
    secondaryFiles: ${return self.location + ext()}


baseCommand: [samtools]
arguments:
  - valueFrom: sort
    position: 1
    # -l - position 2
    # -o sortOutputFileName - position 3
    # -n - position 4
    # --threads - position 5
    # sortInput - position 6
  - valueFrom: ";"
    position: 7
    shellQuote: false
  - valueFrom: samtools
    position: 8
    shellQuote: false
  - valueFrom: index
    position: 9
    # -b - position 10
    # -c - position 11
    # -m - position 12
  - valueFrom: |
      ${
          if (!inputs.sortOutputFileName){
            return default_bam();
          } else {
            return inputs.sortOutputFileName;
          }
      }
    position: 13
  - valueFrom: |
      ${
          if (!inputs.sortOutputFileName){
            return default_bam() + ext();
          } else {
            return inputs.sortOutputFileName + ext();
          }
      }
    position: 14
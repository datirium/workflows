#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: ./metadata/envvar-global.yml
- class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/fastqc:v0.11.5
  dockerFile: >
    $import: ./dockerfiles/fastqc-Dockerfile

inputs:
  fastq_file:
    type:
      - File
      - type: array
        items: File
    inputBinding:
      position: 50

  extract:
    type:
      - "null"
      - boolean
    default: True
    doc: |
      Do not uncompress the output file after creating it.  You
      should set this option if you do not wish to uncompress
      the output when running in non-interactive mode.
      Default is false.

  format:
    type:
      - "null"
      - type: enum
        name: "format"
        symbols: ['bam','sam','bam_mapped','sam_mapped','fastq']
    inputBinding:
      position: 6
      prefix: '--format'
    doc: |
      Bypasses the normal sequence file format detection and
      forces the program to use the specified format.  Valid
      formats are bam,sam,bam_mapped,sam_mapped and fastq

  threads:
    type:
      - "null"
      - int
    inputBinding:
      position: 7
      prefix: '--threads'
    doc: |
      Specifies the number of files which can be processed
      simultaneously.  Each thread will be allocated 250MB of
      memory so you shouldn't run more threads than your
      available memory will cope with, and not more than
      6 threads on a 32 bit machine

  contaminants:
    type:
      - "null"
      - File
    inputBinding:
      position: 8
      prefix: '--contaminants'
    doc: |
      Specifies a non-default file which contains the list of
      contaminants to screen overrepresented sequences against.
      The file must contain sets of named contaminants in the
      form name[tab]sequence.  Lines prefixed with a hash will
      be ignored.

  adapters:
    type:
      - "null"
      - File
    inputBinding:
      position: 9
      prefix: '--adapters'
    doc: |
      Specifies a non-default file which contains the list of
      adapter sequences which will be explicity searched against
      the library. The file must contain sets of named adapters
      in the form name[tab]sequence.  Lines prefixed with a hash
      will be ignored.

  limits:
    type:
      - "null"
      - File
    inputBinding:
      position: 10
      prefix: '--limits'
    doc: |
      Specifies a non-default file which contains a set of criteria
      which will be used to determine the warn/error limits for the
      various modules.  This file can also be used to selectively
      remove some modules from the output all together.  The format
      needs to mirror the default limits.txt file found in the
      Configuration folder.

  kmers:
    type:
      - "null"
      - int
    inputBinding:
      position: 11
      prefix: '--kmers'
    doc: |
      Specifies the length of Kmer to look for in the Kmer content
      module. Specified Kmer length must be between 2 and 10. Default
      length is 7 if not specified.

  casava:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 13
      prefix: '--casava'
    doc: |
      Files come from raw casava output. Files in the same sample
      group (differing only by the group number) will be analysed
      as a set rather than individually. Sequences with the filter
      flag set in the header will be excluded from the analysis.
      Files must have the same names given to them by casava
      (including being gzipped and ending with .gz) otherwise they
      won't be grouped together correctly.

  nofilter:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 14
      prefix: '--nofilter'
    doc: |
      If running with --casava then don't remove read flagged by
      casava as poor quality when performing the QC analysis.

  hide_group:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 15
      prefix: '--nogroup'
    doc: |
      Disable grouping of bases for reads >50bp. All reports will
      show data for every base in the read.  WARNING: Using this
      option will cause fastqc to crash and burn if you use it on
      really long reads, and your plots may end up a ridiculous size.
      You have been warned!

baseCommand: [fastqc, --outdir, .]
arguments:
  - valueFrom:
      ${
        if (inputs.extract){
          return "--extract"
        } else {
          return "--noextract"
        }
      }
    position: 5

outputs:
  zippedFile:
    type:
      - File[]
    outputBinding:
      glob: '*.zip'
  summary_file:
    type:
      - "null"
      - File[]
    outputBinding:
      glob: |
        ${
          return "*/summary.txt";
        }

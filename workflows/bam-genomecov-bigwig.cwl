#!/usr/bin/env cwl-runner
#
# Author: Andrey.Kartashov@cchmc.org (http://orcid.org/0000-0001-9102-5681) / Dr. Barski Lab / Cincinnati Children’s Hospital Medical Center
# Developed for CWL consortium http://commonwl.org/

cwlVersion: v1.0
class: Workflow

requirements:
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement

doc:
  creates genome coverage bigWig file from .bam file

inputs:
  input:
    type: File
  genomeFile:
    type: File
  scale:
    type: float?
  mappedreads:
    type: double?
  pairchip:
    type: boolean?
  fragmentsize:
    type: int?
  strand:
    type: string?
  bigWig:
    type: string?
  split:
    type: boolean?

outputs:
  outfile:
    type: File
    outputSource: bigwig/bigWigOut
  bed_file:
    type: File
    outputSource: genomecov/genomecoverage

steps:
  genomecov:
    run: ../../tools/bedtools-genomecov.cwl
    in:
      input: input
      genomeFile: genomeFile
      dept:
        default: "-bg"
      split:
        source: split
        valueFrom: |
          ${
            if (self == null){
              return true;
            } else {
              return self;
            }
          }
      pairchip: pairchip
      fragmentsize: fragmentsize
      scale: scale
      mappedreads: mappedreads
      strand: strand
    out: [genomecoverage]

  sort:
    run: ../../tools/linux-sort.cwl
    in:
      input: genomecov/genomecoverage
      key:
        default: ["1,1","2,2n"]
    out: [sorted]

  bigwig:
    run: ../../tools/ucsc-bedgraphtobigwig.cwl
    in:
      input: sort/sorted
      genomeFile: genomeFile
      bigWig: bigWig
    out: [bigWigOut]

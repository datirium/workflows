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
    type: int?
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
    outputSource: sort/sorted

steps:
  genomecov:
    run: ../tools/bedtools-genomecov.cwl
    in:
      input_file: input
      chrom_length_file: genomeFile # maybe we don't need it at all, because we always use BAM file as input
      depth:
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
      fragment_size: fragmentsize
      scale: scale
      mapped_reads_number: mappedreads
      strand: strand
    out: [genome_coverage_file]

  sort:
    run: ../tools/linux-sort.cwl
    in:
      input: genomecov/genome_coverage_file
      key:
        default: ["1,1","2,2n"]
    out: [sorted]

  bigwig:
    run: ../tools/ucsc-bedgraphtobigwig.cwl
    in:
      input: sort/sorted
      genomeFile: genomeFile
      bigWig: bigWig
    out: [bigWigOut]

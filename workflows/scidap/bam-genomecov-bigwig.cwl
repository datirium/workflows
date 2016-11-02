#!/usr/bin/env cwl-runner
#
# Author: Andrey.Kartashov@cchmc.org (http://orcid.org/0000-0001-9102-5681) / Dr. Barski Lab / Cincinnati Children’s Hospital Medical Center
# Developed for CWL consortium http://commonwl.org/

cwlVersion: v1.0
class: Workflow


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
    type: string

outputs:
  outfile:
    type: File
    outputSource: bigwig/bigWigOut

steps:
  genomecov:tools/bedtools-genomecov.cwl
    in:
      input: input
      genomeFile: genomeFile
      genomecoverageout:
        default: "./genomecov.bedGraph"
      dept:
        default: "-bg"
      split:
        default: true
      pairchip: pairchip
      fragmentsize: fragmentsize
      scale: scale
      mappedreads: mappedreads
      strand: strand
    out: [genomecoverage]

  sort:
    run: ../tools/linux-sort.cwl
    in:
      input: genomecov/genomecoverage
      key:
        default: ["1,1","2,2n"]
      output:
        default: tmp_sorted
    out: [sorted]

  bigwig:
    run: ../tools/ucsc-bedgraphtobigwig.cwl
    in:
      input: sort/sorted
      genomeFile: genomeFile
      bigWig: bigWig
    out: [bigWigOut]

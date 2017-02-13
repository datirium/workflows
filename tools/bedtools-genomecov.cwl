#!/usr/bin/env cwl-runner

$namespaces:
  dct: http://purl.org/dc/terms/
  foaf: http://xmlns.com/foaf/0.1/
  doap: http://usefulinc.com/ns/doap#
  adms: http://www.w3.org/ns/adms#
  dcat: http://www.w3.org/ns/dcat#

$schemas:
- http://dublincore.org/2012/06/14/dcterms.rdf
- http://xmlns.com/foaf/spec/20140114.rdf
- http://usefulinc.com/ns/doap#
- http://www.w3.org/ns/adms#
- http://www.w3.org/ns/dcat.rdf

cwlVersion: v1.0
class: CommandLineTool

adms:includedAsset:
  doap:name: "bedtools"
  doap:doc: |
    A software suite for the comparison, manipulation and annotation of genomic features in browser extensible data (BED) and general feature format (GFF) format.
    BEDTools also supports the comparison of sequence alignments in BAM format to both BED and GFF features.
    The tools are extremely efficient and allow the user to compare large datasets (e.g. next-generation sequencing data) with both public and custom genome annotation tracks.
    BEDTools can be combined with one another as well as with standard UNIX commands, thus facilitating routine genomics tasks as well as pipelines that can quickly answer intricate questions of large genomic datasets.
  doap:homepage: "http://bedtools.readthedocs.org"
  doap:repository:
  - class: doap:GitRepository
    doap:location: "https://github.com/arq5x/bedtools2"
  doap:release:
  - class: doap:Version
    doap:revision: "2.25.0"
  doap:license: "GPLv2"
  doap:category: "commandline tool"
  doap:programming-language: "C++"
  foaf:publications:
  - id: urn:pmid:20110278
    foaf:title: "Aaron R. Quinlan, Ira M. Hall (2010) BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics, 26(6) 841-842, http://dx.doi.org/10.1093/bioinformatics/btq033"
    foaf:homepage: "http://bioinformatics.oxfordjournals.org/content/26/6/841"
  doap:maintainer:
  - class: foaf:Person
    foaf:name: "Aaron R. Quinlan"
    foaf:mbox: "aaronquinlan at gmail.com"
    dct:isPartOf:
    - class: foaf:Organization
      foaf:name: "Department of Biochemistry and Molecular Genetics, University of Virginia School of Medicine"
    - class: foaf:Organization
      foaf:name: "Center for Public Health Genomics, University of Virginia, Charlottesville, VA 22908, USA"

doc: |
  bedtools-genomecov.cwl is developed for CWL consortium

  Original tool usage:
      Tool:    bedtools genomecov (aka genomeCoverageBed)
      Sources: https://github.com/arq5x/bedtools2
      Summary: Compute the coverage of a feature file among a genome.
      Usage: bedtools genomecov [OPTIONS] -i <bed/gff/vcf> -g <genome>

doap:name: "bedtools-genomecov.cwl"
dcat:downloadURL: "https://github.com/common-workflow-language/workflows/blob/master/tools/bedtools-genomecov.cwl"

doap:maintainer:
- class: foaf:Organization
  foaf:name: "Barski Lab, Cincinnati Children's Hospital Medical Center"
  foaf:member:
  - class: foaf:Person
    id: "http://orcid.org/0000-0001-9102-5681"
    foaf:openid: "http://orcid.org/0000-0001-9102-5681"
    foaf:name: "Andrey Kartashov"
    foaf:mbox: "mailto:Andrey.Kartashov@cchmc.org"

requirements:
  - $import: ./metadata/envvar-global.yml
  - class: InlineJavascriptRequirement

hints:
- class: DockerRequirement
    #dockerImageId: scidap/star:v2.5.0b #not yet ready
  dockerPull: scidap/bedtools2:v2.25.0-pefs2
  dockerFile: >
    $import: ./dockerfiles/bedtools-Dockerfile

inputs:
  input:
    type: File
    doc: |
      The input file can be in BAM format
          (Note: BAM _must_ be sorted by position)
      or <bed/gff/vcf>
    inputBinding:
      position: 1
      valueFrom: |
          ${
            var prefix = ((/.*\.bam$/i).test(inputs.input.path))?'-ibam':'-i';
            return [prefix,inputs.input.path];
          }
# TODO Need this only to make sure that bam file is sorted.
    secondaryFiles: |
           ${
            if ((/.*\.bam$/i).test(self.location))
               return {"location": self.location+".bai", "class": "File"};
            return [];
           }
#    secondaryFiles: |
#      ${
#      if ((/.*\.bam$/i).test(self.location)){
#        if (self.secondaryFiles.length > 0){
#          var secFilesSet = [];
#          for (var i = 0; i < self.secondaryFiles.length; i++){
#            secFilesSet.push ({"location": self.secondaryFiles[i].location, "class": "File"});
#          }
#          return secFilesSet;
#        }
#        return {"location": self.location+".bai", "class": "File"};
#      }
#      return [];
#      }

# Needed only when -i flag
  genomeFile:
    type: File?
    doc:
      Input genome file.
    inputBinding:
      position: 2
      prefix: "-g"

  dept:
    type:
      name: "JustDepts"
      type: enum
      symbols: ["-bg","-bga","-d"]
    inputBinding:
      position: 4

  scale:
    type: float?
    doc: |
      Scale the coverage by a constant factor.
      Each coverage value is multiplied by this factor before being reported.
      Useful for normalizing coverage by, e.g., reads per million (RPM).
      - Default is 1.0; i.e., unscaled.
      - (FLOAT)
    inputBinding:
      position: 4
      prefix: -scale

  mappedreads:
    type: double?
    doc: |
      Optional parameter to calculate scale as 1000000/mappedreads
    inputBinding:
      position: 4
      prefix: -scale
      valueFrom: |
        ${
          if (inputs.scale){
            return null;
          } else {
            return 1000000/inputs.mappedreads;
          }
        }

  dz:
    type: ["null",boolean]
    doc: |
      Report the depth at each genome position (with zero-based coordinates).
      Reports only non-zero positions.
      Default behavior is to report a histogram.
    inputBinding:
      position: 4
      prefix: "-dz"

  split:
    type: ["null",boolean]
    doc: |
      treat "split" BAM or BED12 entries as distinct BED intervals.
      when computing coverage.
      For BAM files, this uses the CIGAR "N" and "D" operations
      to infer the blocks for computing coverage.
      For BED12 files, this uses the BlockCount, BlockStarts, and BlockEnds
      fields (i.e., columns 10,11,12).
    inputBinding:
      position: 4
      prefix: "-split"

  strand:
    type: ["null", string]
    doc: |
      Calculate coverage of intervals from a specific strand.
      With BED files, requires at least 6 columns (strand is column 6).
      - (STRING): can be + or -
    inputBinding:
      position: 4
      prefix: "-strand"

  pairchip:
    type: ["null", boolean]
    doc: "pair-end chip seq experiment"
    inputBinding:
      position: 4
      prefix: "-pc"

  fragmentsize:
    type: ["null", int]
    doc: "fixed fragment size"
    inputBinding:
      position: 4
      prefix: "-fs"

  max:
    type: ["null",int]
    doc: |
      Combine all positions with a depth >= max into
      a single bin in the histogram. Irrelevant
      for -d and -bedGraph
      - (INTEGER)
    inputBinding:
      position: 4
      prefix: "-max"

  m5:
    type: ["null",boolean]
    doc: |
      Calculate coverage of 5" positions (instead of entire interval).
    inputBinding:
      position: 4
      prefix: "-5"

  m3:
    type: ["null",boolean]
    doc: |
      Calculate coverage of 3" positions (instead of entire interval).
    inputBinding:
      position: 4
      prefix: "-3"

  genomecoverageout:
    type: string

outputs:
  genomecoverage:
    type: File
    doc: "The file containing the genome coverage"
    outputBinding:
      glob: $(inputs.genomecoverageout)

stdout: $(inputs.genomecoverageout)

baseCommand: ["bedtools", "genomecov"]

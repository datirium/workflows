cwlVersion: v1.0
class: Workflow

requirements:
  - class: SubworkflowFeatureRequirement
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement

inputs:
#STAR.alighReads.cwl inputs
  readFileInLeft:
    type: File
    format: http://edamontology.org/format_1931 #FASTQ-illumina
  readFileInRight:
    type: File
    format: http://edamontology.org/format_1931 #FASTQ-illumina
  starIndex/genomeDir:
    type: Directory
  clip3pNbases:
    type: int?
  clip5pNbases:
    type: int?
#bam-genomecov-bigwig.cwl inputs
  chrLengthFile:
    type: File
#bowtie.cwl inputs
  bowtieIndex/bowtieEbwtFile:
    type: File
    format:
      - http://edamontology.org/format_3484 # ebwt
      - http://edamontology.org/format_3491 # ebwtl
#reads-counting.cwl inputs
  annotation:
    type: File

outputs:
  outfileBigWig:
    type: File
    outputSource: bamToBigwig/outfile
  alignedLog:
    type: File
    outputSource: alignReads/alignedLog
  fastxStats:
    type: File[]
    outputSource: fastxQualitystats/statistics
  bamBaiPair:
    type: File
    outputSource: samtoolsSortIndex/bamBaiPair
  output_bowtie_log:
    type: File
    outputSource: bowtie/output_bowtie_log
  rpkmFile:
    type: File
    outputSource: rpkm/rpkmFile


steps:
  alignReads:
    run: ../../tools/star-alignreads.cwl
    in:
      readFilesIn: [readFileInLeft, readFileInRight]
      genomeDir: starIndex/genomeDir
      outFilterMultimapNmax:
        default: 1
      outFilterMismatchNmax:
        default: 5
      alignSJDBoverhangMin:
        default: 1
      seedSearchStartLmax:
        default: 15
      clip3pNbases: clip3pNbases
      clip5pNbases: clip5pNbases
    out: [aligned, alignedLog]

  fastxQualitystats:
    run: ../../tools/fastx-quality-stats.cwl
    in:
      inputFile: [readFileInLeft, readFileInRight]
      outputFileDir:
        default: "."
    out: [statistics]
    scatter: inputFile

  samtoolsSortIndex:
    run: ../../tools/samtools-sort-index.cwl
    in:
      sortInput: alignReads/aligned
      sortOutputFileName:
        default: "sorted.bam"
    out: [bamBaiPair]

  bamtoolsStats:
    run: ../../tools/bamtools-stats.cwl
    in:
      inputFiles: samtoolsSortIndex/bamBaiPair
    out: [mappedreads]

  bamToBigwig:
    run: bam-genomecov-bigwig.cwl
    in:
      input: samtoolsSortIndex/bamBaiPair
      genomeFile: chrLengthFile # chamnge the name to the normal one
      mappedreads: bamtoolsStats/mappedreads
      bigWig:
        default: "output.bigwig"
    out: [outfile]

  bowtie:
    run: ../../tools/bowtie.cwl
    in:
      ebwtFile: bowtieIndex/bowtieEbwtFile
      filelist: readFileInLeft
      filelist_mates: readFileInRight
      filename:
        default: "bowtie"
      3p: clip3pNbases
      5p: clip5pNbases
      q:
        default: true
      v:
        default: 3
      m:
        default: 1
      best:
        default: false    # should be true only when not-paired mode
      strata:
        default: false    # should be true only when not-paired mode
      sam:
        default: true
    out: [output_bowtie_log]

  rpkm:
    run: ../../tools/reads-counting.cwl
    in:
      aligned: samtoolsSortIndex/bamBaiPair
      annotation: annotation
      rpkm-cutoff:
        default: 0.001
      rpkm-cutoff-val:
        default: 0
      rnaSeqType:
        default: RNA
    out: [rpkmFile]
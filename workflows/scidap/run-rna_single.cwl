cwlVersion: v1.0
class: Workflow

requirements:
  - class: SubworkflowFeatureRequirement
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement

inputs:
#STAR.alighReads.cwl inputs
  readFilesIn:
    type: File[]
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



steps:
  alignReads:
    run: ../star/starAlignReads/STAR.alignReads.cwl
    in:
      readFilesIn: readFilesIn
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
    run: ../fastx-toolkit/fastxQualitystats/fastxQualitystats.cwl
    in:
      inputFile: readFilesIn
      outputFileDir:
        default: "."
    out: [statistics]
    scatter: inputFile

  samtoolsSortIndex:
    run: ../samtools/samtoolsSortIndex/samtools-sort-index.cwl
    in:
      sortInput: alignReads/aligned
      sortOutputFileName:
        default: "sorted.bam"
    out: [bamBaiPair]

  bamtoolsStats:
    run: ../bamtools/bamtoolsStats/bamtoolsStats.cwl
    in:
      inputFiles: samtoolsSortIndex/bamBaiPair
    out: [mappedreads]

  bamToBigwig:
    run: ../bigwigWorkflow/bam-genomecov-bigwig.cwl
    in:
      input: samtoolsSortIndex/bamBaiPair
      genomeFile: chrLengthFile
      mappedreads: bamtoolsStats/mappedreads
      bigWig:
        default: "output.bigwig"
    out: [outfile]

  bowtie:
    run: ../bowtie/bowtie.cwl
    in:
      ebwtFile: bowtieIndex/bowtieEbwtFile
      filelist: readFilesIn
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
        default: true    # should be true only when not-paired mode
      strata:
        default: true    # should be true only when not-paired mode
      sam:
        default: true
    out: [output_bowtie_log]














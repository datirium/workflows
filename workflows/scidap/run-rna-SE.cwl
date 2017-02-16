cwlVersion: v1.0
class: Workflow

### TODO
# Add option to set threads number


requirements:
  - class: SubworkflowFeatureRequirement
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement

inputs:

  fastq_input_file:
    type: File
    label: "FASTQ input file"
    format: "http://edamontology.org/format_1930"
    doc: "SDLoadFile"

  __genome__organism__nucleus__star:
    type: Directory
    label: "STAR indices folder"
    doc: "SDSelect"

  clip_3p_nbases:
    type: int
    label: "Clip from 3p end"
    doc: "SDInput"

  clip_5p_nbases:
    type: int
    label: "Clip from 5p end"
    doc: "SDInput"

  __genome__organism__star__chrlength:
    type: File
    label: "Chromosome length file"
    format: "http://edamontology.org/format_2330"
    doc: "SDSelect"

  __genome__organism__ribosome__bowtie:
    type: Directory
    label: "BOWTIE indices folder"
    format:
      - http://edamontology.org/format_3484 # ebwt
      - http://edamontology.org/format_3491 # ebwtl
    doc: "SDSelect"

  __genome__organism__annotation__tsv:
    type: File
    label: "TAB-separated annoation file"
    format: "http://edamontology.org/format_3475"
    doc: "SDSelect"

  dutp:
    type: boolean
    label: "dUTP"
    doc: "SDCheckbox"

  spike:
    type: boolean
    label: "Use RNA spike-in sequences"
    doc: "SDCheckbox"

outputs:
  outfile_big_wig:
    type: File
    format: "http://edamontology.org/format_3006"
    label: "BigWig file"
    doc: "SDVisBigWig"
    outputSource: bam_to_bigwig/outfile

  aligned_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "STAR alignment statistics"
    doc: "SDVisStarStat"
    outputSource: star_align_reads/alignedLog

  fastx_stats:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "fastx_quality_stats statistics"
    doc: "SDVisFastxStat"
    outputSource: fastx_quality_stats/statistics

  bambai_pair:
    type: File
    format: "http://edamontology.org/format_2572"
    label: "Coordinate sorted BAM alignment file"
    doc: "SDVisBam"
    outputSource: samtools_sort_index/bamBaiPair

  output_bowtie_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "Bowtie alignment log"
    doc: "SDVisText"  # To display plain text
    outputSource: bowtie/output_bowtie_log

  rpkm_file:
    type: File
    format:
    - "http://edamontology.org/format_3752" # csv
    - "http://edamontology.org/format_3475" # tsv
    label: "RPKM list values"
    doc: "SDVisTable"  # To display csv/tsv file as a table
    outputSource: rpkm/rpkmFile

steps:
  star_align_reads:
    run: ../../tools/star-alignreads.cwl
    in:
      readFilesIn: fastq_input_file
      genomeDir: __genome__organism__nucleus__star
      outFilterMultimapNmax:
        default: 1
      outFilterMismatchNmax:
        default: 5
      alignSJDBoverhangMin:
        default: 1
      seedSearchStartLmax:
        default: 15
      clip3pNbases: clip_3p_nbases
      clip5pNbases: clip_5p_nbases
    out: [aligned, alignedLog]

  fastx_quality_stats:
    run: ../../tools/fastx-quality-stats.cwl
    in:
      inputFile: fastq_input_file
      outputFileDir:
        default: "."
    out: [statistics]

  samtools_sort_index:
    run: ../../tools/samtools-sort-index.cwl
    in:
      sortInput: star_align_reads/aligned
      sortOutputFileName:
        default: "sorted.bam"
    out: [bamBaiPair]

  bamtools_stats:
    run: ../../tools/bamtools-stats.cwl
    in:
      inputFiles: samtools_sort_index/bamBaiPair
    out: [mappedreads]

  bam_to_bigwig:
    run: bam-genomecov-bigwig.cwl
    in:
      input: samtools_sort_index/bamBaiPair
      genomeFile: __genome__organism__star__chrlength
      mappedreads: bamtools_stats/mappedreads
      bigWig:
        default: "output.bigwig"
    out: [outfile]

  bowtie:
    run: ../../tools/bowtie.cwl
    in:
      indices_folder: __genome__organism__ribosome__bowtie
      filelist: fastq_input_file
      filename:
        default: "bowtie"
      3p: clip_3p_nbases
      5p: clip_5p_nbases
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

  rpkm:
    run: ../../tools/reads-counting.cwl
    in:
      aligned: samtools_sort_index/bamBaiPair
      annotation: __genome__organism__annotation__tsv
      rpkm-cutoff:
        default: 0.001
      rpkm-cutoff-val:
        default: 0
      rnaSeqType:
        source: dutp
        valueFrom: |
          ${
            if (self){
              return 'dUTP';
            } else {
              return 'RNA';
            }
          }
      ignore_chrom:
        source: spike
        valueFrom: |
          ${
            if (self){
              return null;
            } else {
              return 'control';
            }
          }
    out: [rpkmFile]


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "run-rna-SE"
s:downloadUrl: https://raw.githubusercontent.com/SciDAP/workflows/v0.0.1b/workflows/scidap/run-rna-SE.cwl
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

s:about: >
  Current workflow should be used only with the single-end RNA-Seq data. It performs the following steps:
  1. Uses STAR to align reads from input FASTQ file according to the predefined reference indices; generates unsorted BAM file and alignment statistics file
  2. Usess fastx_quality_stats to analyze input FASTQ file and generates quality statistics file
  3. Uses samtools sort to generate coordinate sorted BAM(+BAI) file pair from the unsorted BAM file obtained on the step 1 (after running STAR)
  4. Calculates basic statistics for sorted BAM file using bamtools stats. Returns mapped reads number.
  5. Generates BigWig file on the base of sorted BAM file
  6. Maps input FASTQ file to predefined rRNA reference indices using Bowtie to define the level of rRNA contamination; exports resulted statistics to file
  7. Calculates isoform and gene expression levels for the sorted BAM file and tab-separated annotation file using SciDAP reads-counting utility; exports results to file
  Tasks #2 and #6 are started independently on other workflow steps.
  Workflow doesn't depend on the type of the organism











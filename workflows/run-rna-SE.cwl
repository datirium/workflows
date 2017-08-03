cwlVersion: v1.0
class: Workflow

requirements:
  - class: SubworkflowFeatureRequirement
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement

'sd:metadata':
  - "https://raw.githubusercontent.com/SciDAP/workflows/master/metadata/rnaseq-header.cwl"

inputs:

  fastq_input_file:
    type: File
    label: "FASTQ input file"
    format: "http://edamontology.org/format_1930"
    doc: "Reads data in a FASTQ format"

  star_indices:
    type: Directory
    label: "STAR indices folder"
    'sd:parent': "https://raw.githubusercontent.com/SciDAP/workflows/master/workflows/star-index.cwl"
    doc: "STAR generated indices"

  clip_3p_end:
    type: int?
    'sd:layout':
      advanced: true
    label: "Clip from 3p end"
    doc: "Number of bases to clip from the 3p end"

  clip_5p_end:
    type: int?
    'sd:layout':
      advanced: true
    label: "Clip from 5p end"
    doc: "Number of bases to clip from the 5p end"

  chrom_length:
    type: File
    label: "Chromosome length file"
    format: "http://edamontology.org/format_2330"
    'sd:parent': "https://raw.githubusercontent.com/SciDAP/workflows/master/workflows/star-index.cwl"
    doc: "Chromosome length file"

  bowtie_indices:
    type: Directory
    label: "BowTie Ribosomal Indices"
    format:
      - http://edamontology.org/format_3484 # ebwt
      - http://edamontology.org/format_3491 # ebwtl
    'sd:parent': "https://raw.githubusercontent.com/SciDAP/workflows/master/workflows/bowtie-index.cwl"
    doc: "Bowtie generated indices"

  annotation:
    type: File
    label: "GTF annotation file"
    format: "http://edamontology.org/format_2306"
    'sd:parent': "https://raw.githubusercontent.com/SciDAP/workflows/master/workflows/star-index.cwl"
    doc: "GTF annotation file"

  dutp:
    type: boolean
    label: "dUTP"
    doc: "Enable strand specific dUTP calculation"

  threads:
    type: int?
    'sd:layout':
      advanced: true
    label: "Number of threads to run tools"
    doc: "Number of threads for those steps that support multithreading"

outputs:
  bigwig:
    type: File
    format: "http://edamontology.org/format_3006"
    label: "BigWig file"
    doc: "Generated BigWig file"
    outputSource: bam_to_bigwig/outfile

  star_reads_alignment_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "STAR alignment log"
    doc: "STAR alignment log file"
    outputSource: star_reads_alignment/alignedLog

  fastx_quality_stats_statistics:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "fastx_quality_stats statistics"
    doc: "Fastx statistics file"
    outputSource: fastx_quality_stats/statistics

  bambai_pair:
    type: File
    format: "http://edamontology.org/format_2572"
    label: "Coordinate sorted BAM alignment file (+index BAI)"
    doc: "Coordinate sorted BAM file and BAI index file"
    outputSource: samtools_sort_index/bam_bai_pair

  bowtie_reads_alignment_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "Bowtie alignment log"
    doc: "Bowtie alignment log file"  # To display plain text
    outputSource: bowtie_reads_alignment/output_bowtie_log

  rpkm_calculation_table:
    type: File
    format:
    - "http://edamontology.org/format_3752" # csv
    - "http://edamontology.org/format_3475" # tsv
    label: "RPKM table file"
    doc: "Calculated rpkm values"  # To display csv/tsv file as a table
    outputSource: rpkm_calculation/rpkmFile

steps:
  star_reads_alignment:
    run: ../tools/star-alignreads.cwl
    in:
      readFilesIn: fastq_input_file
      genomeDir: star_indices
      outFilterMultimapNmax:
        default: 1
      outFilterMismatchNmax:
        default: 5
      alignSJDBoverhangMin:
        default: 1
      seedSearchStartLmax:
        default: 15
      clip3pNbases: clip_3p_end
      clip5pNbases: clip_5p_end
      threads: threads
    out: [aligned, alignedLog]

  fastx_quality_stats:
    run: ../tools/fastx-quality-stats.cwl
    in:
      input_file: fastq_input_file
    out: [statistics]

  samtools_sort_index:
    run: ../tools/samtools-sort-index.cwl
    in:
      sort_input: star_reads_alignment/aligned
      threads: threads
    out: [bam_bai_pair]

  bamtools_stats:
    run: ../tools/bamtools-stats.cwl
    in:
      input_files: samtools_sort_index/bam_bai_pair
    out: [mappedreads]

  bam_to_bigwig:
    run: bam-genomecov-bigwig.cwl
    in:
      input: samtools_sort_index/bam_bai_pair
      genomeFile: chrom_length
      mappedreads: bamtools_stats/mappedreads
    out: [outfile]

  bowtie_reads_alignment:
    run: ../tools/bowtie.cwl
    in:
      indices_folder: bowtie_indices
      filelist: fastq_input_file
      filename:
        source: fastq_input_file
        valueFrom: |
          ${
            return self.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+'_ribosomal'+'.sam';
          }
      clip_3p_end: clip_3p_end
      clip_5p_end: clip_5p_end
      q:
        default: true
      v:
        default: 3
      m:
        default: 1
      best:
        default: true
      strata:
        default: true
      sam:
        default: true
      threads: threads
    out: [output_bowtie_log]

  rpkm_calculation:
    run: ../tools/reads-counting.cwl
    in:
      aligned: samtools_sort_index/bam_bai_pair
      annotation: annotation
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
        default: "control"
      threads: threads
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

doc: >
  Runs RNA-Seq BioWardrobe basic analysis with single-end data file.

s:about: |
  The original [BioWardrobe's](https://biowardrobe.com) [PubMed ID:26248465](https://www.ncbi.nlm.nih.gov/pubmed/26248465)
  **RNA-Seq** basic analysis for a **single-end** experiment.
  A corresponded input [FASTQ](http://maq.sourceforge.net/fastq.shtml) file has to be provided.

  Current workflow should be used only with the single-end RNA-Seq data. It performs the following steps:
  1. Use STAR to align reads from input FASTQ file according to the predefined reference indices; generate unsorted BAM file and alignment statistics file
  2. Use fastx_quality_stats to analyze input FASTQ file and generate quality statistics file
  3. Use samtools sort to generate coordinate sorted BAM(+BAI) file pair from the unsorted BAM file obtained on the step 1 (after running STAR)
  4. Calculate basic statistics for sorted BAM file using bamtools stats. Return mapped reads number.
  5. Generate BigWig file on the base of sorted BAM file
  6. Map input FASTQ file to predefined rRNA reference indices using Bowtie to define the level of rRNA contamination; export resulted statistics to file
  7. Calculate isoform and gene expression levels for the sorted BAM file and GTF annotation file using SciDAP reads-counting utility; export results to file

  Tasks #2 and #6 are started independently on other workflow steps.
  Workflow doesn't depend on the type of the organism











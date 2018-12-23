cwlVersion: v1.0
class: Workflow

requirements:
  - class: SubworkflowFeatureRequirement
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement

'sd:metadata':
  - "../metadata/rnaseq-header.cwl"

'sd:upstream':
  genome_indices: "genome-indices.cwl"

inputs:

  star_indices_folder:
    type: Directory
    label: "STAR indices folder"
    'sd:upstreamSource': "genome_indices/star_indices"
    doc: "Path to STAR generated indices"

  star_indices_folder_mitochondrial:
    type: Directory
    label: "STAR indices mitochondrial folder"
    'sd:upstreamSource': "genome_indices/mitochondrial_indices"
    doc: "Path to STAR generated indices for mitochondrial dna"

  bowtie_indices_folder:
    type: Directory
    label: "BowTie Ribosomal Indices"
    'sd:upstreamSource': "genome_indices/ribosomal_indices"
    doc: "Path to Bowtie generated indices"

  indices_folder:
    type: Directory
    'sd:upstreamSource': "genome_indices/bowtie_indices"
    label: "Indexed genome folder (bowtie)"
    doc: "Path to indexed genome folder by **bowtie**"

  annotation_file:
    type: File
    'sd:upstreamSource': "genome_indices/annotation"
    label: "Annotation file"
    format: "http://edamontology.org/format_3475"
    doc: "Tab-separated annotation file"

  genome_size:
    type: string
    'sd:upstreamSource': "genome_indices/genome_size"
    label: "Effective genome size"
    doc: "MACS2 effective genome size: hs, mm, ce, dm or number, for example 2.7e9"

  chrom_length:
    type: File
    'sd:upstreamSource': "genome_indices/chrom_length"
    label: "Chromosomes length file"
    format: "http://edamontology.org/format_2330"
    doc: "Chromosomes length file"

  fastq_file:
    type: File
    label: "FASTQ input file"
    format: "http://edamontology.org/format_1930"
    doc: "Reads data in a FASTQ format, received after single end sequencing"

#  output_file:
#    type: string

# ADVANCED

  extract_method:
    type: enum
    default: "regex"
    'sd:layout':
      advanced: true
    label: "UMI extract method 'string' or 'regex'"
    symbols:
      - "string"
      - "regex"
    doc: |
      How to extract the umi +/- cell barcodes, Choose from
      'string' or 'regex'

  bc_pattern:
    type: string
    default: "(?P<umi_1>.{4})G{1}.*"
    'sd:layout':
      advanced: true
    label: "Barcode pattern"

  adapter:
    type: string
    default: "GTGTCAGTCACTTCCAGCGGG"
    'sd:layout':
      advanced: true
    label: "Adapter sequence to be trimmed"
    doc: |
      Adapter sequence to be trimmed. If not specified explicitly, Trim Galore will
      try to auto-detect whether the Illumina universal, Nextera transposase or Illumina
      small RNA adapter sequence was used. Also see '--illumina', '--nextera' and
      '--small_rna'. If no adapter can be detected within the first 1 million sequences
      of the first file specified Trim Galore defaults to '--illumina'.

  remove_duplicates:
    type: boolean?
    default: false
    'sd:layout':
      advanced: true
    label: "Remove duplicates"
    doc: "Calls samtools rmdup to remove duplicates from sortesd BAM file"

  exclude_chr:
    type: string?
    'sd:layout':
      advanced: true
    label: "Chromosome to be excluded in rpkm calculation"
    doc: "Chromosome to be excluded in rpkm calculation"

  clip_3p_end:
    type: int?
    default: 0
    'sd:layout':
      advanced: true
    label: "Clip from 3p end"
    doc: "Number of bases to clip from the 3p end"

  clip_5p_end:
    type: int?
    default: 0
    'sd:layout':
      advanced: true
    label: "Clip from 5p end"
    doc: "Number of bases to clip from the 5p end"


# SYSTEM DEPENDENT

  threads:
    type: int?
    default: 2
    'sd:layout':
      advanced: true
    doc: "Number of threads for those steps that support multithreading"
    label: "Number of threads"

outputs:

  output:
    type: File
    label: "clipped file"
    doc: "clipped fastq file"
    outputSource: extract_umi/output

  error_log:
    type: File
    label: "clipped error log file"
    doc: "clipped error log file"
    outputSource: extract_umi/error_log

  extract_log:
    type: File
    label: "clipped extract log file"
    doc: "clipped extract log file"
    outputSource: extract_umi/log

  star_final_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "STAR final log"
    doc: "STAR Log.final.out"
    outputSource: star_aligner/log_final

  star_out_log:
    type: File?
    format: "http://edamontology.org/format_2330"
    label: "STAR log out"
    doc: "STAR Log.out"
    outputSource: star_aligner/log_out

  star_progress_log:
    type: File?
    format: "http://edamontology.org/format_2330"
    label: "STAR progress log"
    doc: "STAR Log.progress.out"
    outputSource: star_aligner/log_progress

  star_stdout_log:
    type: File?
    format: "http://edamontology.org/format_2330"
    label: "STAR stdout log"
    doc: "STAR Log.std.out"
    outputSource: star_aligner/log_std

  star_sj_log:
    type: File?
    format: "http://edamontology.org/format_2330"
    label: "STAR sj log"
    doc: "STAR SJ.out.tab"
    outputSource: star_aligner/log_sj

  fastx_statistics_original:
    type: File
    label: "FASTQ statistics"
    format: "http://edamontology.org/format_2330"
    doc: "fastx_quality_stats generated FASTQ file quality statistics file"
    outputSource: fastx_quality_stats_original/statistics_file
    'sd:visualPlugins':
    - line:
      Title: 'Base frequency plot'
      xAxisTitle: 'Nucleotide position'
      yAxisTitle: 'Frequency'
      colors: ["#b3de69", "#99c0db", "#fb8072", "#fdc381", "#888888"]
      data: [$12, $13, $14, $15, $16]

  fastx_statistics_after:
    type: File
    label: "FASTQ statistics"
    format: "http://edamontology.org/format_2330"
    doc: "fastx_quality_stats generated FASTQ file quality statistics file"
    outputSource: fastx_quality_stats_after/statistics_file
    'sd:visualPlugins':
    - line:
      Title: 'Base frequency plot'
      xAxisTitle: 'Nucleotide position'
      yAxisTitle: 'Frequency'
      colors: ["#b3de69", "#99c0db", "#fb8072", "#fdc381", "#888888"]
      data: [$12, $13, $14, $15, $16]

  trim_report:
    type: File
    label: "TrimGalore report"
    doc: "TrimGalore generated log"
    outputSource: trim_fastq/report_file


steps:

  extract_fastq:
    run: ../tools/extract-fastq.cwl
    in:
      compressed_file: fastq_file
    out: [fastq_file]

  fastx_quality_stats_original:
    run: ../tools/fastx-quality-stats.cwl
    in:
      input_file: extract_fastq/fastq_file
    out: [statistics_file]

  extract_umi:
    run: ../tools/umi_tools-extract.cwl
    in:
      input_file: extract_fastq/fastq_file
      extract_method: extract_method
      bc_pattern: bc_pattern
    out: [output, log, error_log]

  trim_fastq:
    run: ../tools/trimgalore.cwl
    in:
      input_file: extract_umi/output
      dont_gzip:
        default: true
      length:
        default: 30
    out: [trimmed_file, report_file]

#  rename:
#    run: ../tools/rename.cwl
#    in:
#      source_file: trim_fastq/trimmed_file
#      target_filename:
#        source: extract_fastq/fastq_file
#        valueFrom: $(self.basename)
#    out:
#      - target_file

  fastx_quality_stats_after:
    run: ../tools/fastx-quality-stats.cwl
    in:
      input_file: trim_fastq/trimmed_file
    out: [statistics_file]

  star_aligner:
    run: ../tools/star-alignreads.cwl
    in:
      readFilesIn: trim_fastq/trimmed_file
      genomeDir: star_indices_folder
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
    out:
      - aligned_file
      - log_final
      - uniquely_mapped_reads_number
      - log_out
      - log_progress
      - log_std
      - log_sj



$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "CLIP-Seq pipeline for single-read experiment NNNNG"
label: "CLIP-Seq pipeline for single-read experiment NNNNG"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/clipseq-se.cwl
s:codeRepository: https://github.com/datirium/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Datirium, LLC"
  s:member:
  - class: s:Person
    s:name: Artem BArski
    s:email: mailto:Artem.Barski@datirum.com
  - class: s:Person
    s:name: Andrey Kartashov
    s:email: mailto:Andrey.Kartashov@datirium.com
    s:sameAs:
    - id: http://orcid.org/0000-0001-9102-5681

doc: |
  CLIP-Seq workflow for single-read experiment.

s:about: |
  '''CLIP''' ('''cross-linking immunoprecipitation''') is a method used in molecular biology that combines UV cross-linking with
  immunoprecipitation in order to analyse protein interactions with RNA or to precisely locate RNA modifications (e.g. m6A).
  (Uhl|Houwaart|Corrado|Wright|Backofen|2017)(Ule|Jensen|Ruggiu|Mele|2003)(Sugimoto|König|Hussain|Zupan|2012)(Zhang|Darnell|2011)
  (Ke| Alemu| Mertens| Gantman|2015) CLIP-based techniques can be used to map RNA binding protein binding sites or RNA modification
  sites (Ke| Alemu| Mertens| Gantman|2015)(Ke| Pandya-Jones| Saito| Fak|2017) of interest on a genome-wide scale,
  thereby increasing the understanding of post-transcriptional regulatory networks.

  The identification of sites where RNA-binding proteins (RNABPs) interact with target RNAs opens the door to understanding
  the vast complexity of RNA regulation. UV cross-linking and immunoprecipitation (CLIP) is a transformative technology in which RNAs
  purified from ~in vivo~ cross-linked RNA-protein complexes are sequenced to reveal footprints of RNABP:RNA contacts.
  CLIP combined with high-throughput sequencing (HITS-CLIP) is a generalizable strategy to produce transcriptome-wide maps of RNA
  binding with higher accuracy and resolution than standard RNA immunoprecipitation (RIP) profiling or purely computational approaches.

  The application of CLIP to Argonaute proteins has expanded the utility of this approach to mapping binding sites for microRNAs
  and other small regulatory RNAs. Finally, recent advances in data analysis take advantage of cross-link–induced mutation sites
  (CIMS) to refine RNA-binding maps to single-nucleotide resolution. Once IP conditions are established, HITS-CLIP takes ~8 d to prepare
  RNA for sequencing. Established pipelines for data analysis, including those for CIMS, take 3–4 d.

  ==Workflow==
  CLIP begins with the in-vivo cross-linking of RNA-protein complexes using ultraviolet light (UV).
  Upon UV exposure, covalent bonds are formed between proteins and nucleic acids that are in close proximity.
  (Darnell|2012) The cross-linked cells are then lysed, and the protein of interest is isolated via immunoprecipitation.
  In order to allow for sequence specific priming of reverse transcription, RNA adapters are ligated to the 3' ends,
  while radiolabeled phosphates are transferred to the 5' ends of the RNA fragments.
  The RNA-protein complexes are then separated from free RNA using gel electrophoresis and membrane transfer.
  Proteinase K digestion is then performed in order to remove protein from the RNA-protein complexes.
  This step leaves a peptide at the cross-link site, allowing for the identification of the cross-linked nucleotide.
  (König| McGlincy| Ule|2012) After ligating RNA linkers to the RNA 5' ends, cDNA is synthesized via RT-PCR.
  High-throughput sequencing is then used to generate reads containing distinct barcodes that identify the last cDNA nucleotide.
  Interaction sites can be identified by mapping the reads back to the transcriptome.










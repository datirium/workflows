cwlVersion: v1.0
class: Workflow


requirements:
  - class: StepInputExpressionRequirement
  - class: SubworkflowFeatureRequirement


'sd:upstream':
  genome_indices: ["genome-indices.cwl", "https://github.com/datirium/workflows/workflows/genome-indices.cwl"]


inputs:

  alias:
    type: string
    label: "Sample short name/Alias:"
    sd:preview:
      position: 1

  condition:
    type: string?
    label: "Experimental condition:"
    sd:preview:
      position: 2

  cells:
    type: string
    label: "Cells:"
    sd:preview:
      position: 3

  catalog:
    type: string?
    label: "Catalog No.:"
    sd:preview:
      position: 4

  chrom_length_file:
    type: File
    label: "Chromosome length file"
    format: "http://edamontology.org/format_2330"
    'sd:upstreamSource': "genome_indices/chrom_length"
    doc: "Chromosome length file"

  index_directory:
    type: Directory
    'sd:upstreamSource': "genome_indices/bowtie_indices"
    label: "Bowtie2 index:"
    'sd:localLabel': true
    doc: |
      Bowtie2 index directory of the reference genome for alignment and IGV visualization.
    sd:preview:
      position: 5

  refgenome:
    type: File
    'sd:upstreamSource': "genome_indices/fasta_output"
    label: "Reference Genome FASTA:"
    'sd:localLabel': true
    doc: |
      Reference genome FASTA file to be used for alignment.
    sd:preview:
      position: 6

  genome:
    type:
    - "null"
    - type: enum
      name: "genome short name"
      symbols:
      - hg19
      - hg38
      - mm10
      - rn7
      - dm3
    label: "Genome short name:"
    'sd:localLabel': true
    doc: |
      Genome short name used for setting organism name, genus, species, and tax ID.
    sd:preview:
      position: 7

  fastq_file:
    type:
      - File
      - type: array
        items: File
    label: "Input FASTQ file:"
    'sd:localLabel': true
    format: "http://edamontology.org/format_1930"
    doc: |
      FASTQ file from a single-end miRNA sequencing run.
    sd:preview:
      position: 11

  adapter:
    type: string?
    default: TCGTAT
    label: "Adapter:"
    'sd:localLabel': true
    doc: |
      Adapter sequence to be trimmed from miRNA sequence reads (Default: TCGTAT).
    'sd:layout':
      advanced: true

  threads:
    type: int?
    default: 4
    label: "Threads:"
    'sd:localLabel': true
    doc: |
      Number of threads to use for steps that support multithreading (Default: 4).
    'sd:layout':
      advanced: true


outputs:

  fastx_statistics:
    type: File
    label: "FASTQ statistics"
    format: "http://edamontology.org/format_2330"
    doc: "fastx_quality_stats generated quality statistics file for input FASTQ"
    outputSource: fastx_quality_stats/statistics_file
    'sd:visualPlugins':
    - line:
        tab: 'QC Plots'
        Title: 'FASTQ Base frequency plot'
        xAxisTitle: 'Nucleotide position'
        yAxisTitle: 'Frequency'
        colors: ["#b3de69", "#888888", "#fb8072", "#fdc381", "#99c0db"]
        data: [$13, $14, $15, $16, $17]
    - boxplot:
        tab: 'QC Plots'
        Title: 'FASTQ Quality Control'
        xAxisTitle: 'Nucleotide position'
        yAxisTitle: 'Quality score'
        colors: ["#b3de69", "#888888", "#fb8072", "#fdc381", "#99c0db"]
        data: [$11, $7, $8, $9, $12]

  bambai_pair:
    type: File
    format: "http://edamontology.org/format_2572"
    label: "Coordinate sorted BAM alignment file (+index BAI)"
    doc: "Coordinate sorted BAM file and BAI index file"
    outputSource: samtools_sort_index/bam_bai_pair
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        optional: true
        type: 'alignment'
        format: 'bam'
        name: "BAM Track"
        height: 120
        displayMode: "EXPANDED"

  bowtie_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "Bowtie alignment log"
    doc: "Bowtie alignment log file"
    outputSource: bowtie_aligner/log_file

  preseq_estimates:
    type: File?
    label: "Preseq estimates"
    format: "http://edamontology.org/format_3475"
    doc: "Preseq estimated results"
    outputSource: preseq/estimates_file
    'sd:visualPlugins':
    - scatter:
        tab: 'QC Plots'
        Title: 'Preseq Estimates'
        xAxisTitle: 'Total reads count'
        yAxisTitle: 'Expected distinct reads count'
        colors: ["#4b78a3"]
        height: 500
        data: [$1, $2]
        comparable: "preseq"

  mirdeep2_result:
    type: File
    format: "http://edamontology.org/format_2331"
    label: "Report: miRDeep2 results"
    doc: "standard mirdeep2 results html report of novel and known mirs in your data"
    outputSource: mirdeep2/mirdeep2_result
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  mirs_known:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "known mature miRNAs"
    doc: "known mature miRNAs detected by mirdeep2"
    outputSource: mirdeep2/mirs_known
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Known miRNAs'

  rpkm_isoforms:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "Input file for DESeq workflow"
    doc: "This is the input file for DESeq workflow, the output name is a misnomer: isoforms, genes, and common tss are all the same."
    outputSource: mirdeep2/deseq_input_isoforms

  rpkm_genes:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "Input file for DESeq workflow"
    doc: "This is the input file for DESeq workflow, the output name is a misnomer: isoforms, genes, and common tss are all the same."
    outputSource: mirdeep2/deseq_input_genes

  rpkm_common_tss:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "Input file for DESeq workflow"
    doc: "This is the input file for DESeq workflow, the output name is a misnomer: isoforms, genes, and common tss are all the same."
    outputSource: mirdeep2/deseq_input_common_tss

  mirs_novel:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "known novel miRNAs"
    doc: "known novel miRNAs detected by mirdeep2"
    outputSource: mirdeep2/mirs_novel
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Novel miRNAs'

  mirs_known_exocarta_deepmirs:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "known mature miRNA overlapping with ExoCarta exosome miRNA"
    doc: "known mature miRNA overlapping with ExoCarta exosome miRNA"
    outputSource: mirdeep2/mirs_known_exocarta_deepmirs
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Detected Exosome miRNAs'

  mirs_known_gene_targets:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "gene targets of detected known mature miRNA"
    doc: "gene targets of detected known mature miRNA"
    outputSource: mirdeep2/mirs_known_gene_targets

  known_mirs_mature:
    type: File
    format: "http://edamontology.org/format_1929"
    label: "FASTA of known mature miRNA sequences"
    doc: "FASTA of known mature miRNA sequences"
    outputSource: mirdeep2/known_mirs_mature

  known_mirs_precursor:
    type: File
    format: "http://edamontology.org/format_1929"
    label: "FASTA of known precursor miRNA sequences"
    doc: "FASTA of known precursor miRNA sequences"
    outputSource: mirdeep2/known_mirs_precursor

  novel_mirs_mature:
    type: File
    format: "http://edamontology.org/format_1929"
    label: "FASTA of novel mature miRNA sequences"
    doc: "FASTA of novel mature miRNA sequences"
    outputSource: mirdeep2/novel_mirs_mature

  novel_mirs_precursor:
    type: File
    format: "http://edamontology.org/format_1929"
    label: "FASTA of novel precursor miRNA sequences"
    doc: "FASTA of novel precursor miRNA sequences"
    outputSource: mirdeep2/novel_mirs_precursor

  overview:
    type: File
    format: "http://edamontology.org/format_3835"
    label: "mirdeep2 script metrics"
    doc: "markdown parsed metrics from run_mirdeep2.sh script run by mirna-mirdeep2-se.cwl"
    outputSource: mirdeep2/overview
    'sd:visualPlugins':
    - markdownView:
        tab: 'Overview'

  mirna_mirdeep2_stdout:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "stdout logfile"
    doc: "captures standard output from vc-germline-pe.cwl"
    outputSource: mirdeep2/log_file_stdout

  mirna_mirdeep2_stderr:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "stderr logfile"
    doc: "captures standard error from vc-germline-pe.cwl"
    outputSource: mirdeep2/log_file_stderr


steps:

  extract_fastq:
    label: "Loading unmapped sequence data for read 1"
    doc: |
      Most DNA cores and commercial NGS companies return unmapped sequence data in FASTQ format.
      The data can be uploaded from users computer, downloaded directly from an ftp server of
      the core facility by providing a URL or from GEO by providing SRA accession number.
    run: ../tools/extract-fastq.cwl
    in:
      compressed_file: fastq_file
      output_prefix:
        default: "merged_fastq_file"
    out: [fastq_file]

  fastx_quality_stats:
    label: "Quality control of unmapped sequence data for read 1"
    doc: |
      Evaluates the quality of your sequence data. Provides per base quality scores as well as
      base frequencies along the reads. These metrics can be used to identify whether your data
      has any problems that should be taken into account in the subsequent analysis steps.
    run: ../tools/fastx-quality-stats.cwl
    in:
      input_file: extract_fastq/fastq_file
    out: [statistics_file]

  trim_fastq:
    run: ../tools/trimgalore.cwl
    in:
      input_file: extract_fastq/fastq_file
      adapter: adapter
      dont_gzip:
        default: true
      length:
        default: 30
    out:
      - trimmed_file
      - report_file

# The m param is set to mirdeep mapper.pl default
  bowtie_aligner:
    run: ../tools/bowtie-alignreads.cwl
    in:
      upstream_filelist: trim_fastq/trimmed_file
      indices_folder: index_directory
      v:
        default: 3
      m:
        default: 5
      unaligned_prefix:
        default: "unaligned_reads"
      best:
        default: true
      strata:
        default: true
      sam:
        default: true
      threads: threads
    out: [sam_file, log_file, unaligned_fastq, multimapped_fastq]

  samtools_sort_index:
    run: ../tools/samtools-sort-index.cwl
    in:
      sort_input: bowtie_aligner/sam_file
      threads: threads
    out: [bam_bai_pair]

  clean_sam_headers_for_preseq:
    run: ../tools/samtools-clean-headers.cwl
    in:
      bam_file: samtools_sort_index/bam_bai_pair
    out: [preseq_bam]

  preseq:
    label: "Sequencing depth estimation"
    doc: |
      Estimates the complexity of the sequencing library, evaluates how many reads can
      be expected from the additional sequencing of the same experiment.
    run: ../tools/preseq-lc-extrap.cwl
    in:
      bam_file: clean_sam_headers_for_preseq/preseq_bam
      pe_mode:
        default: true
      extrapolation:
        default: 1000000000
    out: [estimates_file]

  mirdeep2:
    label: "Run pipeline for calling germline variants, and read, alignment, and variant stat generation"
    doc: |
      Calls shell wrapper for the SciDAP miRDeep2 pipeline.
    run: ../tools/mirna-mirdeep2-se.cwl
    in:
      threads: threads
      bt2dir: index_directory
      refgenome: refgenome
      genome:
        source: genome
        valueFrom: $(self)
      adapter: adapter
      fastq: extract_fastq/fastq_file
    out: [mirs_known, deseq_input_isoforms, deseq_input_genes, deseq_input_common_tss, mirs_novel, mirs_known_exocarta_deepmirs, mirs_known_gene_targets, known_mirs_mature, known_mirs_precursor, novel_mirs_mature, novel_mirs_precursor, overview, mirdeep2_result, log_file_stdout, log_file_stderr]


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "miRNA-Seq miRDeep2 pipeline"
label: "miRNA-Seq miRDeep2 pipeline"
s:alternateName: "miRNA-Seq workflow using miRDeep2 for discovering known and novel miRNA from RNA-Seq data."

s:downloadUrl: https://github.com/datirium/workflows/tree/master/workflows/workflows/mirna-mirdeep2-se.cwl
s:codeRepository: https://github.com/datirium/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Datirium LLC"
  s:location:
  - class: s:PostalAddress
    s:addressCountry: "USA"
    s:addressLocality: "Cincinnati"
    s:addressRegion: "OH"
    s:postalCode: ""
    s:streetAddress: ""
    s:telephone: ""
  s:logo: "https://avatars.githubusercontent.com/u/33202955?s=200&v=4"
  s:department:
  - class: s:Organization
    s:legalName: "Datirium LLC"
    s:department:
    - class: s:Organization
      s:legalName: "Bioinformatics"
      s:member:
      - class: s:Person
        s:name: Robert Player
        s:email: mailto:support@datirium.com
        s:sameAs:
        - id: https://orcid.org/0000-0001-5872-259X


doc: |
  A CWL workflow for discovering known or novel miRNAs from deep sequencing data using the miRDeep2 tool.
  The ExoCarta exosome database is also used for identifying exosome-related miRNAs, and TargetScan's
  organism-specific databases are used for identifying miRNA gene targets.

  ## __Outputs__
  #### Primary Output files:
    - mirs_known.tsv, detected known mature miRNAs, "Known miRNAs" tab
    - mirs_novel.tsv, detected novel mature miRNAs, "Novel miRNAs" tab
  #### Secondary Output files:
    - mirs_known_exocarta_deepmirs.tsv, list of detected miRNA also in ExoCarta's exosome database, "Detected Exosome miRNAs" tab
    - mirs_known_gene_targets.tsv, pre-computed gene targets of known mature mirs, downloadable
    - known_mirs_mature.fa, known mature mir sequences, downloadable
    - known_mirs_precursor.fa, known precursor mir sequences, downloadable
    - novel_mirs_mature.fa, novel mature mir sequences, downloadable
    - novel_mirs_precursor.fa, novel precursor mir sequences, downloadable
  #### Reports:
    - overview.md (input list, alignment & mir metrics), "Overview" tab
    - mirdeep2_result.html, summary of mirdeep2 results, "miRDeep2 Results" tab

  ## __Inputs__
  #### General Info
   - Sample short name/Alias: unique name for sample
   - Experimental condition: condition, variable, etc name (e.g. "control" or "20C 60min")
   - Cells: name of cells used for the sample
   - Catalog No.: vender catalog number if available
   - Bowtie2 index: Bowtie2 index directory of the reference genome.
   - Reference Genome FASTA: Reference genome FASTA file to be used for alignment.
   - Genome short name: Name used for setting organism name, genus, species, and tax ID.
   - Input FASTQ file: FASTQ file from a single-end miRNA sequencing run.
  #### Advanced
   - Adapter: Adapter sequence to be trimmed from miRNA sequence reads. (Default: TCGTAT)
   - Threads: Number of threads to use for steps that support multithreading (Default: 4).

  ## Hints & Tips:
  
  #### For the identification of novel miRNA candidates, the following may be used as a filtering guideline:
  1. miRDeep score > 4 (some authors use 1)
  2. not present a match with rfam
  3. should present a significant RNAfold ("yes")
  4. a number of mature reads > 10
  5. if applicable, novel mir must be expressed in multiple samples

  #### For filtering mirbase by organism.

  | genome | organism | division | name | tree | NCBI-taxid |
  | ---- | --- | --- | ----------- | ----------- | ----------- |
  | hg19 | hsa | HSA | Homo sapiens | Metazoa;Bilateria;Deuterostoma;Chordata;Vertebrata;Mammalia;Primates;Hominidae | 9606 |
  | hg38 | hsa | HSA | Homo sapiens | Metazoa;Bilateria;Deuterostoma;Chordata;Vertebrata;Mammalia;Primates;Hominidae | 9606 |
  | mm10 | mmu | MMU | Mus musculus | Metazoa;Bilateria;Deuterostoma;Chordata;Vertebrata;Mammalia;Rodentia | 10090 |
  | rn7  | rno | RNO | Rattus norvegicus | Metazoa;Bilateria;Deuterostoma;Chordata;Vertebrata;Mammalia;Rodentia | 10116 |
  | dm3	 | dme | DME | Drosophila melanogaster | Metazoa;Bilateria;Ecdysozoa;Arthropoda;Hexapoda | 7227 |

  ## __Data Analysis Steps__
  1. The miRDeep2 Mapper module processes Illumina FASTQ output and maps it to the reference genome.
  2. The miRDeep2 miRDeep2 module identifies known and novel (mature and precursor) miRNAs.
  3. The ExoCarta database of miRNA found in exosomes is then used to find overlap between mirs_known.tsv and exosome associated miRNAs.
  4. Finally, TargetScan organism-specific miRNA gene target database is used to find overlap between mirs_known.tsv and gene targets.

  ## __References__
  1. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3245920
  2. https://github.com/rajewsky-lab/mirdeep2
  3. https://biocontainers.pro/tools/mirdeep2
  4. https://www.mirbase.org/
  5. http://exocarta.org/index.html
  6. https://www.targetscan.org/vert_80/

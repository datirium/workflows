cwlVersion: v1.2.0-dev4
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: MultipleInputFeatureRequirement
  - class: InlineJavascriptRequirement
    expressionLib:
    - var get_root = function(basename) {
          return basename.split('.').slice(0,1).join('.');
      };


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

  bowtie_indices_folder:
    type: Directory
    label: "BowTie Ribosomal Indices"
    'sd:upstreamSource': "genome_indices/ribosomal_indices"
    doc: "Path to Bowtie generated indices"

  chrom_length_file:
    type: File
    label: "Chromosome length file"
    format: "http://edamontology.org/format_2330"
    'sd:upstreamSource': "genome_indices/chrom_length"
    doc: "Chromosome length file"

  annotation_file:
    type: File
    label: "Annotation file"
    format:
      - "http://edamontology.org/format_2306"
      - "http://edamontology.org/format_3475"
    'sd:upstreamSource': "genome_indices/annotation"
    doc: "GTF or TAB-separated annotation file"

  annotation_gtf_file:
    type: File
    label: "GTF annotation file"
    format: "http://edamontology.org/format_2306"
    'sd:upstreamSource': "genome_indices/annotation_gtf"
    doc: "GTF annotation file"

  fastq_file:
    type:
    - File
    - type: array
      items: File
    label: "FASTQ input file"
    format: "http://edamontology.org/format_1930"
    doc: "Reads data in a FASTQ format"

  use_umi:
    type: boolean?
    default: true
    'sd:layout':
      advanced: true
    label: "Use UMIs"
    doc: "Use UMIs (for FWD-UMI libraries)"

  min_length:
    type: int?
    default: 30
    'sd:layout':
      advanced: true
    label: "Set minimum length for trimmed reads when running FWD/REV pipeline. Shorter reads get discarded. Set 0 to disable"
    doc: |
      Set minimum length for trimmed reads when running FWD/REV (not UMI) pipeline.
      Shorter reads get discarded. Applied only when running trim_fastq step.
      For FWD-UMI pipeline we use cutadapt instead of TrimGalore, so this input is
      not used

  exclude_chr:
    type: string?
    default: ""
    'sd:layout':
      advanced: true
    label: "Coma-separated list of chromosomes to be excluded from gene expression calculation"
    doc: "Coma-separated list of chromosomes to be excluded from gene expression calculation"

  clip_3p_end:
    type: int?
    default: 0
    'sd:layout':
      advanced: true
    label: "Clip N bp from 3p end"
    doc: "Number of bp to clip from the 3p end"

  clip_5p_end:
    type: int?
    default: 0
    'sd:layout':
      advanced: true
    label: "Clip N bp from 5p end"
    doc: "Number of bp to clip from the 5p end"

  threads:
    type: int?
    default: 1
    'sd:layout':
      advanced: true
    label: "Number of threads"
    doc: "Number of threads for those steps that support multi-threading"


outputs:

  bigwig:
    type: File
    format: "http://edamontology.org/format_3006"
    label: "BigWig file"
    doc: "Generated BigWig file"
    outputSource: bam_to_bigwig/bigwig_file
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'wig'
        name: "BigWig Track"
        height: 120

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

  fastx_statistics:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "FASTQ statistics"
    doc: "fastx_quality_stats generated FASTQ file quality statistics file"
    outputSource: fastx_quality_stats/statistics_file
    'sd:visualPlugins':
    - line:
        tab: 'QC Plots'
        Title: 'Base frequency plot'
        xAxisTitle: 'Nucleotide position'
        yAxisTitle: 'Frequency'
        colors: ["#b3de69", "#888888", "#fb8072", "#fdc381", "#99c0db"]
        data: [$13, $14, $15, $16, $17]
    - boxplot:
        tab: 'QC Plots'
        Title: 'Quality Control'
        xAxisTitle: 'Nucleotide position'
        yAxisTitle: 'Quality score'
        colors: ["#b3de69", "#888888", "#fb8072", "#fdc381", "#99c0db"]
        data: [$11, $7, $8, $9, $12]

  bambai_pair:
    type: File
    format: "http://edamontology.org/format_2572"
    label: "Coordinate sorted BAM alignment file (+index BAI)"
    doc: "Coordinate sorted BAM file and BAI index file"
    outputSource: samtools_sort_index_2/bam_bai_pair
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        optional: true
        type: 'alignment'
        format: 'bam'
        name: "BAM Track"
        displayMode: "SQUISHED"

  bowtie_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "Bowtie alignment log"
    doc: "Bowtie alignment log file"
    outputSource: bowtie_aligner/log_file

  rpkm_genes:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "raw reads grouped by gene name"
    doc: "raw reads grouped by gene name"
    outputSource: group_isoforms/genes_file
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Gene Expression'
        Title: 'raw reads grouped by gene name'

  reads_per_gene_htseq_count:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "Gene expression from htseq-count (reads per gene)"
    doc: "Gene expression from htseq-count (reads per gene)"
    outputSource: htseq_calculate_expression/gene_expression_report

  rpkm_common_tss:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "raw reads grouped by common TSS"
    doc: "raw reads grouped by common TSS"
    outputSource: group_isoforms/common_tss_file

  get_stat_log:
    type: File?
    format: "http://edamontology.org/format_3750"
    label: "YAML formatted combined log"
    doc: "YAML formatted combined log"
    outputSource: get_stat/collected_statistics_yaml

  get_stat_markdown:
    type: File?
    label: "Markdown formatted combined log"
    format: "http://edamontology.org/format_3835"
    doc: "Markdown formatted combined log"
    outputSource: get_stat/collected_statistics_md
    'sd:visualPlugins':
    - markdownView:
        tab: 'Overview'

  get_formatted_stats:
    type: File?
    label: "Bowtie, STAR and GEEP mapping stats"
    format: "http://edamontology.org/format_2330"
    doc: "Processed and combined Bowtie & STAR aligner and GEEP logs"
    outputSource: get_stat/collected_statistics_tsv
    'sd:visualPlugins':
    - tableView:
        vertical: true
        tab: 'Overview'
    'sd:preview':
      'sd:visualPlugins':
      - pie:
          colors: ['#b3de69', '#99c0db', '#fdc381', '#fb8072', '#778899']
          data: [$2, $3, $4, $5, $6]

  bam_statistics_report:
    type: File
    label: "BAM statistics report"
    format: "http://edamontology.org/format_2330"
    doc: "BAM statistics report (right after alignment and sorting)"
    outputSource: get_bam_statistics/log_file

  trimgalore_report:
    type: File?
    format: "http://edamontology.org/format_2330"
    label: "Adapter trimming report from TrimGalore. Even if it was eventually bypassed"
    doc: "Adapter trimming report from TrimGalore. Even if it was eventually bypassed"
    outputSource: trim_fastq/report_file

  cutadapt_report:
    type: File?
    format: "http://edamontology.org/format_2330"
    label: "Adapter trimming report from Cutadapt"
    doc: "Adapter trimming report from Cutadapt"
    outputSource: umisep_cutadapt/report_file

  umi_tools_dedup_stdout:
    type: File?
    format: "http://edamontology.org/format_2330"
    label: "umi_tools dedup stdout log"
    doc: "umi_tools dedup stdout log"
    outputSource: umi_tools_dedup/stdout_log

  umi_tools_dedup_stderr:
    type: File?
    format: "http://edamontology.org/format_2330"
    label: "umi_tools dedup stderr log"
    doc: "umi_tools dedup stderr log"
    outputSource: umi_tools_dedup/stderr_log

  umi_tools_dedup_stats:
    type:
      - "null"
      - File[]
    label: "umi_tools dedup stats"
    doc: "umi_tools dedup stats"
    outputSource: umi_tools_dedup/output_stats


steps:

  extract_fastq:
    run: ../tools/extract-fastq.cwl
    in:
      compressed_file: fastq_file
    out:
    - fastq_file

  trim_fastq:
    when: $(!inputs.use_umi)                  # not sure if ! is valid syntax in this case
    run: ../tools/trimgalore.cwl
    in:
      use_umi: use_umi                        # need it for "when"
      input_file: extract_fastq/fastq_file
      dont_gzip:
        default: true                         # saves time
      length: min_length
    out:
    - trimmed_file
    - report_file

  bypass_trim:
    run: ../tools/bypass-trimgalore-se.cwl
    in:
      original_fastq_file: extract_fastq/fastq_file
      trimmed_fastq_file: trim_fastq/trimmed_file
      trimming_report_file: trim_fastq/report_file
      min_reads_count:
        default: 100                          # any small number should be good, as we are catching the case when TrimGalore discarded all reads
    out:
    - selected_fastq_file

  umisep_cutadapt:
    when: $(inputs.use_umi)
    in:
      use_umi: use_umi                        # need it for "when"
      input_file: extract_fastq/fastq_file
    out:
    - trimmed_file
    - report_file
    run:
      cwlVersion: v1.0
      class: CommandLineTool
      hints:
      - class: DockerRequirement
        dockerPull: scidap/trimgalore:v0.6.6
      inputs:
        bash_script:
          type: string?
          default: |
            #!/bin/bash
            FILE=$0
            BASENAME=$(basename "$FILE")
            cat ${FILE} | awk '
            NR%4==1{ rd_name=$1; rd_info=$2 }
            NR%4==2{ umi=substr($1,1,10); rd_seq=substr($1,11) }
            NR%4==0{ print rd_name"_"umi" "rd_info; print rd_seq; print "+"; print substr($1,11) }' |
            cutadapt -m 20 -O 20 -a "polyA=A{20}" -a "QUALITY=G{20}" -n 2 - |
            cutadapt  -m 20 -O 3 --nextseq-trim=10  -a "r1adapter=A{18}AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=3;max_error_rate=0.100000" - |
            cutadapt -m 20 -O 3 -a "r1polyA=A{18}" - |
            cutadapt -m 20 -O 20 -g "r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=20" --discard-trimmed -o trimmed_${BASENAME} -
          inputBinding:
            position: 1
          doc: "Bash function to run awk & cutadapt from Lexogen"
        input_file:
          type: File
          inputBinding:
            position: 2
          doc: "Input FASTQ file"
      outputs:
        trimmed_file:
          type: File
          outputBinding:
            glob: "trimmed_*"
        report_file:
          type: stderr
      baseCommand: [bash, '-c']
      stderr: umisep_cutadapt_report.txt

  rename:
    run: ../tools/rename.cwl
    in:
      source_file:
        source:
        - umisep_cutadapt/trimmed_file
        - bypass_trim/selected_fastq_file
        pickValue: the_only_non_null          # should be always only one non-null value
      target_filename:
        source: extract_fastq/fastq_file
        valueFrom: $(self.basename)
    out:
      - target_file

  star_aligner:
    run: ../tools/star-alignreads.cwl
    in:
      readFilesIn: rename/target_file
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

  fastx_quality_stats:
    run: ../tools/fastx-quality-stats.cwl
    in:
      input_file: rename/target_file
    out: [statistics_file]

  samtools_sort_index_1:
    run: ../tools/samtools-sort-index.cwl
    in:
      sort_input: star_aligner/aligned_file
      sort_output_filename:
        source: rename/target_file
        valueFrom: $(self.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+'.bam')
      threads: threads
    out: [bam_bai_pair]

  umi_tools_dedup:
    when: $(inputs.use_umi)
    run: ../tools/umi-tools-dedup.cwl
    in:
      use_umi: use_umi                               # need it for "when"
      bam_file: samtools_sort_index_1/bam_bai_pair
      multimapping_detection_method:
        default: "NH"
    out:
    - dedup_bam_file
    - output_stats
    - stdout_log
    - stderr_log
    
  samtools_sort_index_2:                              # easier to run it twice even if umi_tools_dedup was skipped
    run: ../tools/samtools-sort-index.cwl
    in:
      sort_input:
        source:
        - umi_tools_dedup/dedup_bam_file              # will be selected first (if not null)
        - samtools_sort_index_1/bam_bai_pair
        pickValue: first_non_null
      sort_output_filename:
        source: rename/target_file
        valueFrom: $(self.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+'.bam')
      threads: threads
    out: [bam_bai_pair]

  bam_to_bigwig:
    run: ../tools/bam-bedgraph-bigwig.cwl
    in:
      bam_file: samtools_sort_index_2/bam_bai_pair
      chrom_length_file: chrom_length_file
      mapped_reads_number: star_aligner/uniquely_mapped_reads_number
    out: [bigwig_file]

  bowtie_aligner:
    run: ../tools/bowtie-alignreads.cwl
    in:
      upstream_filelist: rename/target_file
      indices_folder: bowtie_indices_folder
      clip_3p_end: clip_3p_end
      clip_5p_end: clip_5p_end
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
    out: [log_file]

  calculate_expression:
    run: ../tools/geep.cwl
    in:
      bam_file: samtools_sort_index_2/bam_bai_pair
      annotation_file: annotation_file
      rpkm_threshold:
        default: 0
      max_cycles:
        default: 0
      exclude_chr: exclude_chr
      threads: threads
    out: [isoforms_file]

  group_isoforms:
    in:
      isoforms_file: calculate_expression/isoforms_file
    out:
    - genes_file
    - common_tss_file
    run:
      cwlVersion: v1.0
      class: CommandLineTool
      hints:
      - class: DockerRequirement
        dockerPull: biowardrobe2/scidap-deseq:v0.0.20
      inputs:
        bash_script:
          type: string?
          default: |
            #!/bin/bash
            FILE=$0
            BASENAME=$(basename "$FILE")
            get_gene_n_tss.R --isoforms "${FILE}" --gene grouped.genes.tsv --tss grouped.common_tss.tsv
            sed -ibak 's/[[:space:]]\{1,\}[^[:space:]]\{1,\}$//' grouped.genes.tsv
            sed -ibak 's/[[:space:]]\{1,\}[^[:space:]]\{1,\}$//' grouped.common_tss.tsv
            rm -f ./*bak
          inputBinding:
            position: 1
          doc: "Bash function to run R script to group expression by genes and common TSS"
        isoforms_file:
          type: File
          inputBinding:
            position: 5
      outputs:
        genes_file:
          type: File
          outputBinding:
            glob: $(inputs.genes_filename?inputs.genes_filename:"*genes.tsv")
          doc: "Output TSV gene expression file"
        common_tss_file:
          type: File
          outputBinding:
            glob: $(inputs.common_tss_file?inputs.common_tss_file:"*common_tss.tsv")
          doc: "Output TSV common tss expression file"
      baseCommand: [bash, '-c']

  get_bam_statistics:
    run: ../tools/samtools-stats.cwl
    in:
      bambai_pair: samtools_sort_index_2/bam_bai_pair
      output_filename:
        source: samtools_sort_index_2/bam_bai_pair
        valueFrom: $(get_root(self.basename)+"_bam_statistics_report.txt")
    out: [log_file]

  get_stat:
      run: ../tools/collect-statistics-rna-quantseq.cwl
      in:
        star_alignment_report: star_aligner/log_final
        bowtie_alignment_report: bowtie_aligner/log_file
        bam_statistics_report: get_bam_statistics/log_file
        isoforms_file: calculate_expression/isoforms_file
      out:
      - collected_statistics_yaml
      - collected_statistics_tsv
      - collected_statistics_md

  htseq_calculate_expression:
    run: ../tools/htseq-count.cwl
    in:
      alignment_bam_file: samtools_sort_index_2/bam_bai_pair
      annotation_gtf_file: annotation_gtf_file
    out:
    - gene_expression_report



$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "QuantSeq 3' FWD, FWD-UMI or REV for single-read mRNA-Seq data"
label: "QuantSeq 3' FWD, FWD-UMI or REV for single-read mRNA-Seq data"
s:alternateName: "Runs QuantSeq 3' FWD, FWD-UMI or REV analysis for single-read mRNA-Seq data"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/trim-quantseq-mrnaseq-se.cwl
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
        s:name: Artem Barski
        s:email: mailto:Artem.Barski@datirum.com
      - class: s:Person
        s:name: Andrey Kartashov
        s:email: mailto:Andrey.Kartashov@datirium.com
        s:sameAs:
        - id: http://orcid.org/0000-0001-9102-5681
      - class: s:Person
        s:name: Michael Kotliar
        s:email: mailto:misha.kotliar@gmail.com
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898


# doc:
#   $include: ../descriptions/trim-quantseq-mrnaseq-se.md


doc: |
  ### Devel version of QuantSeq 3' FWD, FWD-UMI or REV for single-read mRNA-Seq data 
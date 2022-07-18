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

  bowtie_indices_folder:
    type: Directory
    label: "BowTie Ribosomal Indices"
    'sd:upstreamSource': "genome_indices/ribosomal_indices"
    doc: "Path to Bowtie generated indices"

  annotation_file:
    type: File
    'sd:upstreamSource': "genome_indices/annotation"
    label: "Annotation file"
    format: "http://edamontology.org/format_3475"
    doc: "Tab-separated annotation file"

  chrom_length_file:
    type: File
    'sd:upstreamSource': "genome_indices/chrom_length"
    label: "Chromosomes length file"
    format: "http://edamontology.org/format_2330"
    doc: "Chromosomes length file"

  species:
    type: string?
    default: "mm10"
    label: "Species string for clipper (hg38, mm10)"
    doc: "species: one of ce10 ce11 dm3 hg19 GRCh38 mm9 mm10"

  fastq_file:
    type:
    - File
    - type: array
      items: File
    label: "FASTQ input file"
    format: "http://edamontology.org/format_1930"
    doc: "Reads data in a FASTQ format, received after single end sequencing"

  extract_method:
    type:
    - "null"
    - type: enum
      symbols: ["string", "regex"]
    default: "regex"
    'sd:layout':
      advanced: true
    label: "UMI extract method 'string' or 'regex'"
    doc: |
      How to extract the umi +/- cell barcodes, Choose from
      'string' or 'regex'

  bc_pattern:
    type: string?
    default: "(?P<umi_1>.{4})(?P<discard_1>G).*"
    'sd:layout':
      advanced: true
    label: "Barcode pattern"

  adapter:
    type: string?
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

  threads:
    type: int?
    default: 2
    'sd:layout':
      advanced: true
    doc: "Number of threads for those steps that support multithreading"
    label: "Number of threads"


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

  #  output:
  #    type: File
  #    label: "clipped file"
  #    format: "http://edamontology.org/format_1930"
  #    doc: "clipped fastq file"
  #    outputSource: extract_umi/umi_fastq_file_1

  rebosomal_bowtie_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "Bowtie alignment log"
    doc: "Bowtie alignment log file"
    outputSource: ribosomal_bowtie_aligner/log_file

  error_log:
    type: File
    label: "clipped error log file"
    doc: "clipped error log file"
    outputSource: extract_umi/stderr_log

  extract_log:
    type: File
    label: "clipped extract log file"
    doc: "clipped extract log file"
    outputSource: extract_umi/stdout_log

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
        tab: 'QC Plots'
        Title: 'Original Base frequency plot'
        xAxisTitle: 'Nucleotide position'
        yAxisTitle: 'Frequency'
        colors: ["#b3de69", "#888888", "#fb8072", "#fdc381", "#99c0db"]
        data: [$13, $14, $15, $16, $17]
    - boxplot:
        tab: 'QC Plots'
        Title: 'Original Quality Control'
        xAxisTitle: 'Nucleotide position'
        yAxisTitle: 'Quality score'
        colors: ["#b3de69", "#888888", "#fb8072", "#fdc381", "#99c0db"]
        data: [$11, $7, $8, $9, $12]

  fastx_statistics_after:
    type: File
    label: "FASTQ statistics"
    format: "http://edamontology.org/format_2330"
    doc: "fastx_quality_stats generated FASTQ file quality statistics file"
    outputSource: fastx_quality_stats_after/statistics_file
    'sd:visualPlugins':
    - line:
        tab: 'QC Plots'
        Title: 'After Clipper Base frequency plot'
        xAxisTitle: 'Nucleotide position'
        yAxisTitle: 'Frequency'
        colors: ["#b3de69", "#888888", "#fb8072", "#fdc381", "#99c0db"]
        data: [$13, $14, $15, $16, $17]
    - boxplot:
        tab: 'QC Plots'
        Title: 'After Clipper Quality Control'
        xAxisTitle: 'Nucleotide position'
        yAxisTitle: 'Quality score'
        colors: ["#b3de69", "#888888", "#fb8072", "#fdc381", "#99c0db"]
        data: [$11, $7, $8, $9, $12]

  trim_report:
    type: File
    label: "trimm report"
    format: "http://edamontology.org/format_2330"
    doc: "TrimGalore generated log"
    outputSource: trim_fastq/report_file

  bambai_pair:
    type: File
    format: "http://edamontology.org/format_2572"
    label: "Deduped BAM alignment file"
    doc: "Coordinate sorted BAM file and BAI index file (+index BAI)"
    outputSource: samtools_sort_index2/bam_bai_pair
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        optional: true
        type: 'alignment'
        format: 'bam'
        name: "BAM Track"
        displayMode: "SQUISHED"

  # @depricate
  dedup_output:
    type: File
    label: "deduped CLIP file"
    outputSource: dedup_umi/dedup_bam_file

  dedup_error_log:
    type: File
    label: "deduped CLIP error log file"
    doc: "deduped CLIP error log file"
    outputSource: dedup_umi/stderr_log

  dedup_log:
    type: File
    label: "deduped CLIP log file"
    doc: "deduped CLIP log file"
    outputSource: dedup_umi/stdout_log

  output_bed:
    type: File
    outputSource: bamtobed/bed_file

  peaks_bed:
    type: File
    outputSource: tagstopeak/peaks_bed

  get_stat_log:
    type: File?
    label: "Old Bowtie, STAR and GEEP combined log"
    format: "http://edamontology.org/format_2330"
    doc: "Processed and combined Bowtie & STAR aligner and GEEP logs"
    outputSource: stats_and_transformations/output_file

  get_formatted_stats:
    type: File?
    label: "Bowtie, STAR and GEEP mapping stats"
    format: "http://edamontology.org/format_2330"
    doc: "Processed and combined Bowtie & STAR aligner and GEEP logs"
    outputSource: stats_and_transformations/formatted_output_file
    'sd:visualPlugins':
    - tableView:
        vertical: true
        # tab: 'Overview'
    'sd:preview':
      'sd:visualPlugins':
        - pie:
            colors: ['#b3de69', '#99c0db', '#fdc381', '#fb8072']
            data: [$2, $3, $4, $5]

  clipper_bed:
    type: File
    outputSource: clipper/output_bed

  clipper_pickle:
    type: File
    outputSource: clipper/output_pickle

  # Remove in the future BioWardrobe plugs
  atdp_result:
    type: File
    label: "Fake ATDP results for BioWardrobe"
    format: "http://edamontology.org/format_3475"
    doc: "Average Tag Density generated results"
    outputSource: stats_and_transformations/fake_atdp_file

  transformed_peaks:
    type: File
    label: "Transformed peaks Mimics MACS2"
    format: "http://edamontology.org/format_3475"
    outputSource: stats_and_transformations/transformed_peaks

  iaintersect_result:
    type: File
    label: "Island intersect results"
    format: "http://edamontology.org/format_3475"
    doc: "Iaintersect generated results"
    outputSource: island_intersect/result_file



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
    run: ../tools/umi-tools-extract.cwl
    in:
      fastq_file_1: extract_fastq/fastq_file
      extract_method: extract_method
      bc_pattern_1: bc_pattern
      output_filename_1:
        source: extract_fastq/fastq_file
        valueFrom: |
          ${
            var _f = self.location.split('/').slice(-1)[0].split('.');
            return _f.slice(0,-1).join('.') + '_extracted.' + _f.slice(-1)[0];
          }
    out: [umi_fastq_file_1, stdout_log, stderr_log]

  trim_fastq:
    run: ../tools/trimgalore.cwl
    in:
      input_file: extract_umi/umi_fastq_file_1
      adapter: adapter
      dont_gzip:
        default: true
      length:
        default: 30
    out: [trimmed_file, report_file]

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

  ribosomal_bowtie_aligner:
    run: ../tools/bowtie-alignreads.cwl
    in:
      upstream_filelist: trim_fastq/trimmed_file
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

  samtools_sort_index1:
    run: ../tools/samtools-sort-index.cwl
    in:
      sort_input: star_aligner/aligned_file
      sort_output_filename:
        source: extract_fastq/fastq_file
        valueFrom: $(self.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+'.bam')
      threads: threads
    out: [bam_bai_pair]

  dedup_umi:
    run: ../tools/umi-tools-dedup.cwl
    in:
      bam_file: samtools_sort_index1/bam_bai_pair
    out: [dedup_bam_file, stdout_log, stderr_log]

  samtools_sort_index2:
    run: ../tools/samtools-sort-index.cwl
    in:
      sort_input: dedup_umi/dedup_bam_file
      sort_output_filename:
        source: extract_fastq/fastq_file
        valueFrom: $(self.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+'.bam')
      threads: threads
    out: [bam_bai_pair]

  bam_to_bigwig:
    run: ../tools/bam-bedgraph-bigwig.cwl
    in:
      bam_file: samtools_sort_index2/bam_bai_pair
      chrom_length_file: chrom_length_file
      mapped_reads_number: star_aligner/uniquely_mapped_reads_number
      bigwig_filename:
        source: extract_fastq/fastq_file
        valueFrom: $(self.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+'.bigWig')
    #     fragmentsize is not set (STAR gives only read length). It will be calculated automatically by bedtools genomecov.
    out: [bigwig_file]

  bamtobed:
    run: ../tools/bedtools-bamtobed.cwl
    in:
      bam_file: samtools_sort_index2/bam_bai_pair
    out: [bed_file]

  tagstopeak_transformations:
    in:
      annotation: annotation_file
    out: [transformed_annotation]
    run:
      cwlVersion: v1.0
      class: CommandLineTool
      requirements:
        - class: ShellCommandRequirement

      hints:
        - class: DockerRequirement
          dockerPull: biowardrobe2/scidap:v0.0.3

      inputs:
        script:
          type: string?
          default: |
            # !/usr/bin/env python
            import sys, re, math
            with open("transformed_annotation.tsv", 'w') as fof:
                with open(sys.argv[1], 'r') as afile:
                    next(afile) # header line
                    for line in afile:
                      al=line.split()
                      # orig                      3                         6       7         8                            11
                      # bin     name    chrom   strand  txStart txEnd   cdsStart  cdsEnd  exonCount exonStarts  exonEnds  score   name2   cdsStartStat    cdsEndStat      exonFrames
                      # req
                      # chrom chromStart chromEnd name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts
                      blkStarts= ','.join([str(int(x)-int(al[4])) for x in al[9].split(',') if x])
                      blkSizes=','.join([str(-int(e)+int(al[10].split(',')[i])) for i, e in enumerate([x for x in al[9].split(',') if x])])
                      fof.write(al[2]+"\t"+al[4]+"\t"+al[5]+"\t"+al[1]+"\t"+al[11]+"\t"+al[3]+"\t"+al[6]+"\t"+al[7]+"\t0\t"+al[8]+"\t"+blkSizes+"\t"+blkStarts+"\n")

          inputBinding:
            position: 2
        annotation:
          type: File
          inputBinding:
            position: 3

      outputs:
        transformed_annotation:
          type: File
          outputBinding:
            glob: "transformed_annotation.tsv"
      baseCommand: [python, '-c']

  tagstopeak:
    run: ../tools/clip-toolkit-tag2peak.cwl
    in:
      bed_file: bamtobed/bed_file
      big:
        default: true
      separate_strands:
        default: true
      valley_seeking:
        default: true
      gene_bed_file: tagstopeak_transformations/transformed_annotation
    out: [peaks_bed]

  clipper:
    run: ../tools/clipper.cwl
    in:
      input_file: samtools_sort_index2/bam_bai_pair
      species: species
    out: [output_tsv, output_bed, output_pickle]

  stats_and_transformations:
    in:
      star_log: star_aligner/log_final
      bowtie_log: ribosomal_bowtie_aligner/log_file
      dedup_log: dedup_umi/stdout_log
      peaks: clipper/output_bed
      # peaks: tagstopeak/peaks_bed
    out: [output_file, formatted_output_file, fake_atdp_file, transformed_peaks]
    run:
      cwlVersion: v1.0
      class: CommandLineTool
      requirements:
        - class: ShellCommandRequirement
        - class: InlineJavascriptRequirement
          expressionLib:
            - var get_output_filename = function() {
              return inputs.star_log.location.split('/').slice(-1)[0].replace(/_extracted_trimmed\.*Log\.final\.out$/i,'');
              }
      hints:
        - class: DockerRequirement
          dockerPull: biowardrobe2/scidap:v0.0.3

      inputs:
        script:
          type: string?
          default: |
            # !/usr/bin/env python
            import sys, re, math
            TOTAL, ALIGNED, RIBO, MULTIMAPPED, USED = 0, 0, 0, 0, 0
            with open(sys.argv[1], 'r') as star_log:
                for line in star_log:
                    if 'Number of input reads' in line:
                        TOTAL = int(line.split('|')[1])
                    if 'Uniquely mapped reads number' in line:
                        ALIGNED = int(line.split('|')[1])
                    if 'Number of reads mapped to too many loci' in line:
                        MULTIMAPPED = int(line.split('|')[1])
            with open(sys.argv[2], 'r') as bowtie_log:
                for line in bowtie_log:
                    if 'alignment:' in line:
                        RIBO = int(line.split('alignment:')[1].split()[0])
            with open(sys.argv[3], 'r') as dedup_log:
                for line in dedup_log:
                    if 'Number of reads out:' in line:
                        USED = int(line.split('Number of reads out:')[1])

            print TOTAL, ALIGNED, MULTIMAPPED, USED

            with open(sys.argv[4]+"_stats.tsv", 'w') as fof:
                fof.write("Reads total\tReads used\tMulti-mapped\tDuplicates\tUnmapped\tRibosomal contamination\n")
                fof.write(str(TOTAL) + "\t" + str(USED) + "\t" + str(MULTIMAPPED) + "\t" + str(ALIGNED-USED) + "\t" + str(TOTAL-ALIGNED-MULTIMAPPED) + "\t" + str(RIBO) + "\n")

            # TODO: Get rid of! No need without biowardrobe!
            with open(sys.argv[4]+"_atdp.tsv", 'w') as fof:
                fof.write("X\tY\n")
                for i in range(-5000, 5001):
                  fof.write(str(i) + "\t1\n")

            # TODO: Get rid of! No need with right iaintersect! clv toolkit
            #with open(sys.argv[4]+"_macs_peaks.tsv", 'w') as fof:
            #    fof.write("chr\tstart\tend\tlength\tabs_summit\tpileup\t-log10(pvalue)\tfold_enrichment\t-log10(qvalue)\tname\n")
            #    with open(sys.argv[5], 'r') as peak_file:
            #        for line in peak_file:
            #          tmpa=line.split()
            #          pis=[x.split('=')[1] for x in re.split(r'[\[\]]',tmpa[3])[1:] if x.strip()]
            #          fof.write(tmpa[0]+"\t"+tmpa[1]+"\t"+tmpa[2]+"\t"+str(int(tmpa[2])-int(tmpa[1]))+"\t0\t"+pis[1]+"\t"+str(-math.log10(float(pis[3])))+"\t"+tmpa[4]+"\t0\t"+tmpa[3]+"\n")

            # TODO: Get rid of! No need with right iaintersect! clv toolkit
            with open(sys.argv[4]+"_macs_peaks.tsv", 'w') as fof:
                fof.write("chr\tstart\tend\tlength\tabs_summit\tpileup\t-log10(pvalue)\tfold_enrichment\t-log10(qvalue)\tname\n")
                with open(sys.argv[5], 'r') as peak_file:
                    for line in peak_file:
                      tmpa=line.split()
                      fof.write(tmpa[0]+"\t"+tmpa[1]+"\t"+tmpa[2]+"\t"+str(int(tmpa[2])-int(tmpa[1]))+"\t0\t0\t"+str(-math.log10(float(tmpa[4])))+"\t"+tmpa[4]+"\t0\t"+tmpa[3]+"\n")

          inputBinding:
            position: 5

        star_log:
          type: File
          inputBinding:
            position: 6

        bowtie_log:
          type: File
          inputBinding:
            position: 7

        dedup_log:
          type: File
          inputBinding:
            position: 8

        output_filename:
          type:
            - string?
          inputBinding:
            position: 9
            valueFrom: $(get_output_filename())
          default: ""

        peaks:
          type: File
          inputBinding:
            position: 10

      outputs:
        output_file:
          type: stdout

        formatted_output_file:
          type: File
          outputBinding:
            glob: $(get_output_filename()+"_stats.tsv")

        fake_atdp_file:
          type: File
          outputBinding:
            glob: $(get_output_filename()+"_atdp.tsv")

        # faking MACS2 peaks name file
        transformed_peaks:
          type: File
          outputBinding:
            glob: $(get_output_filename()+"_macs_peaks.tsv")

      baseCommand: [python, '-c']

      stdout: $(get_output_filename()+".stat")

  island_intersect:
    run: ../tools/iaintersect.cwl
    in:
      input_filename: stats_and_transformations/transformed_peaks
      annotation_filename: annotation_file
      promoter_bp:
        default: 1000
    out: [result_file, log_file]

$namespaces:
  s: http://schema.org/

$schemas:
  - https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "CLIP-Seq pipeline for single-read experiment NNNNG"
label: "CLIP-Seq pipeline for single-read experiment NNNNG"
s:alternateName: "CLIP-Seq workflow for single-read experiment"

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
  Cross-Linking ImmunoPrecipitation
  =================================

  `CLIP` (`cross-linking immunoprecipitation`) is a method used in molecular biology that combines UV cross-linking with
  immunoprecipitation in order to analyse protein interactions with RNA or to precisely locate RNA modifications (e.g. m6A).
  (Uhl|Houwaart|Corrado|Wright|Backofen|2017)(Ule|Jensen|Ruggiu|Mele|2003)(Sugimoto|König|Hussain|Zupan|2012)(Zhang|Darnell|2011)
  (Ke| Alemu| Mertens| Gantman|2015) CLIP-based techniques can be used to map RNA binding protein binding sites or RNA modification
  sites (Ke| Alemu| Mertens| Gantman|2015)(Ke| Pandya-Jones| Saito| Fak|2017) of interest on a genome-wide scale,
  thereby increasing the understanding of post-transcriptional regulatory networks.

  The identification of sites where RNA-binding proteins (RNABPs) interact with target RNAs opens the door to understanding
  the vast complexity of RNA regulation. UV cross-linking and immunoprecipitation (CLIP) is a transformative technology in which RNAs
  purified from _in vivo_ cross-linked RNA-protein complexes are sequenced to reveal footprints of RNABP:RNA contacts.
  CLIP combined with high-throughput sequencing (HITS-CLIP) is a generalizable strategy to produce transcriptome-wide maps of RNA
  binding with higher accuracy and resolution than standard RNA immunoprecipitation (RIP) profiling or purely computational approaches.

  The application of CLIP to Argonaute proteins has expanded the utility of this approach to mapping binding sites for microRNAs
  and other small regulatory RNAs. Finally, recent advances in data analysis take advantage of cross-link–induced mutation sites
  (CIMS) to refine RNA-binding maps to single-nucleotide resolution. Once IP conditions are established, HITS-CLIP takes ~8 d to prepare
  RNA for sequencing. Established pipelines for data analysis, including those for CIMS, take 3–4 d.

  Workflow
  --------

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
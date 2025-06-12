cwlVersion: v1.0
class: Workflow
requirements:
- class: SubworkflowFeatureRequirement
- class: StepInputExpressionRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var get_root = function(basename) { return basename.split('.').slice(0,1).join('.'); };
sd:upstream:
  epi_sample:
  - chipseq-se.cwl
  - chipseq-pe.cwl
  - trim-chipseq-se.cwl
  - trim-chipseq-pe.cwl
  - trim-atacseq-se.cwl
  - trim-atacseq-pe.cwl
  - cutandrun-macs2-pe.cwl
  - cutandrun-seacr-pe.cwl
  filtered_experiment:
  - filter-peaks-for-heatmap.cwl
  - filter-deseq-for-heatmap.cwl
  - filter-diffbind-for-heatmap.cwl
  - genelists-sets.cwl
inputs:
  alias:
    type: string
    label: Experiment short name/Alias
    sd:preview:
      position: 1
  alignment_file:
    type: File[]
    format: http://edamontology.org/format_2572
    label: Epigenomic sample(s)
    doc: Array of alignment files in BAM format from epigenomic samples selected by user.
    sd:upstreamSource: epi_sample/bambai_pair
    sd:localLabel: true
  alignment_name:
    type: string[]
    label: Epigenomic sample(s)
    doc: Names for input alignment files from epigenomic samples selected by user. Order corresponds to the alignment_file
    sd:upstreamSource: epi_sample/alias
  regions_file:
    type: File
    format: http://edamontology.org/format_3003
    label: Filtered Peaks or DEGs sample
    doc: |
      "Regions of interest from a filtered epigenomic sample or filtered DEGs from a DESeq experiment. Formatted as headerless BED file with [chrom start end name score strand] for gene list and [chrom start end name] for peak file. [name] should be unique, [score] is ignored"
    sd:upstreamSource: filtered_experiment/filtered_file
    sd:localLabel: true
  recentering:
    type:
    - 'null'
    - type: enum
      symbols:
      - Gene TSS
      - Peak Center
    default: Gene TSS
    label: Re-center regions of interest. Choose [Gene TSS] for a gene list or [Peak Center] for a peak file
    doc: Re-center regions of interest. Choose [Gene TSS] for a gene list or [Peak Center] for a peak file
  fragment_size:
    type: int[]
    label: Epigenomic sample(s)
    doc: Array of fragment sizes for input BAM files, order corresponds to the alignment_file
    sd:upstreamSource: epi_sample/estimated_fragment_size
  mapped_reads_number:
    type: int[]
    label: Epigenomic sample(s)
    doc: Array of mapped read numbners for input BAM files, order corresponds to the alignment_file
    sd:upstreamSource: epi_sample/mapped_reads_number
  hist_width:
    type: int?
    default: 10000
    label: Histogram / heatmap width, bp
    doc: Histogram / heatmap width, bp
    sd:layout:
      advanced: true
  hist_bin_size:
    type: int?
    default: 50
    label: Histogram / heatmap bin size, bp
    doc: Histogram / heatmap bin size, bp
    sd:layout:
      advanced: true
  threads:
    type: int?
    default: 4
    label: Number of threads
    doc: Number of threads for steps that support multithreading
    sd:layout:
      advanced: true
outputs:
  heatmap_table:
    type: File
    format: http://edamontology.org/format_3475
    label: TSS or peak centered heatmap as TSV
    doc: TSS or peak centered heatmap as TSV
    outputSource: make_heatmap/histogram_file
  heatmap_plot:
    type: File?
    format: http://edamontology.org/format_3603
    label: TSS or peak centered heatmap as PNG
    doc: TSS or peak centered heatmap as PNG
    outputSource: preview_heatmap/heatmap_png
    sd:visualPlugins:
    - image:
        tab: Plots
        Caption: Tag Enrichment Heatmap
  histogram_table:
    type: File
    format: http://edamontology.org/format_3475
    label: TSS centered average tag density histogram as TSV
    doc: TSS centered average tag density histogram as TSV
    outputSource: preview_histogram/histogram_tsv
  histogram_plot:
    type: File
    format: http://edamontology.org/format_3603
    label: TSS centered average tag density histogram as PNG
    doc: TSS centered average tag density histogram as PNG
    outputSource: preview_histogram/histogram_png
    sd:visualPlugins:
    - image:
        tab: Plots
        Caption: Average Tag Density
  recentered_regions_file:
    type: File
    format: http://edamontology.org/format_3003
    label: Re-centered by [Gene TSS] or [Peak Center] regions of interest file
    doc: Re-centered by [Gene TSS] or [Peak Center] regions of interest file
    outputSource: recenter_regions/output_file
steps:
  make_tag_folders:
    run: ../tools/heatmap-prepare.cwl
    in:
      bam_file: alignment_file
      output_folder: alignment_name
      fragment_size: fragment_size
      total_reads: mapped_reads_number
    out:
    - tag_folder
  recenter_regions:
    run:
      cwlVersion: v1.0
      class: CommandLineTool
      hints:
      - class: DockerRequirement
        dockerPull: biowardrobe2/scidap:v0.0.3
      inputs:
        script:
          type: string?
          default: |
            if [ "$1" == "Gene TSS" ]
            then
              # BED for gene list
              # chrom  start  end  name  [score] strand
              echo "Recenter by the gene TSS"
              cat "$0" | awk '{tss=$2; if ($6=="-") tss=$3; print $1"\t"tss"\t"tss"\ts"$4"\t"$5"\t"$6}' > `basename $0`
            else
              # BED for peaks
              # chrom  start  end  name
              echo "Recenter by the peak center"
              cat "$0" | awk '{center=$2+int(($3-$2)/2); print $1"\t"center"\t"center"\ts"$4"\t"0"\t+"}' > `basename $0`
            fi
          inputBinding:
            position: 1
        input_file:
          type: File
          inputBinding:
            position: 2
        param:
          type:
          - 'null'
          - type: enum
            symbols:
            - Gene TSS
            - Peak Center
          inputBinding:
            position: 3
      outputs:
        output_file:
          type: File
          outputBinding:
            glob: '*'
      baseCommand:
      - bash
      - -c
    in:
      input_file: regions_file
      param: recentering
    out:
    - output_file
  make_heatmap:
    run: ../tools/homer-annotate-peaks-hist.cwl
    in:
      peak_file: recenter_regions/output_file
      tag_folders: make_tag_folders/tag_folder
      hist_width: hist_width
      hist_bin_size: hist_bin_size
      export_heatmap:
        default: true
      threads: threads
      histogram_filename:
        source: regions_file
        valueFrom: $(get_root(self.basename)+"_heatmap.cdt")
    out:
    - histogram_file
  make_histogram:
    run: ../tools/homer-annotate-peaks-hist.cwl
    in:
      peak_file: recenter_regions/output_file
      tag_folders: make_tag_folders/tag_folder
      hist_width: hist_width
      hist_bin_size: hist_bin_size
      export_heatmap:
        default: false
      threads: threads
      histogram_filename:
        source: regions_file
        valueFrom: $(get_root(self.basename)+"_histogram.tsv")
    out:
    - histogram_file
  preview_heatmap:
    in:
      heatmap_file: make_heatmap/histogram_file
      column_names: alignment_name
      output_name:
        source: regions_file
        valueFrom: $(get_root(self.basename)+"_heatmap.png")
    out:
    - heatmap_png
    run:
      cwlVersion: v1.0
      class: CommandLineTool
      requirements:
      - class: DockerRequirement
        dockerPull: biowardrobe2/hopach:v0.0.7
      - class: InitialWorkDirRequirement
        listing:
        - entryname: preview.R
          entry: "#!/usr/bin/env Rscript\noptions(warn=-1)\noptions(\"width\"=300)\nsuppressMessages(library(argparse))\nsuppressMessages(library(RColorBrewer))\nsuppressMessages(library(pheatmap))\nparser <- ArgumentParser(description='Heatmap')\nparser$add_argument(\"--input\",           help='Input CDT file', type=\"character\", required=\"True\")\nparser$add_argument(\"--name\",            help='Input aliases, the order and number corresponds to --input', type=\"character\", required=\"True\", nargs='+')\nparser$add_argument(\"--palette\",         help='Palette color names. Default: black, yellow, white',         type=\"character\", nargs='+', default=c(\"black\", \"yellow\", \"white\"))\nparser$add_argument(\"--output\",          help='Output prefix. Default: heatmap', type=\"character\", default=\"./heatmap.png\")\nargs <- parser$parse_args(commandArgs(trailingOnly = TRUE))\nraw_data <- read.table(args$input, sep=\"\\t\", header=TRUE, stringsAsFactors=FALSE)\ncorrected_data <- raw_data[,-1]\nrownames(corrected_data) <- raw_data[,1]\nprint(\"Centering by mean\")\ncorrected_data = corrected_data - rowMeans(corrected_data)    \nprint(\"Normalizing\")\nstd = sqrt(rowSums(corrected_data^2))\ncorrected_data = corrected_data/std\ncorrected_data = replace(corrected_data, is.na(corrected_data), 0)\ntryCatch(\n  expr = {\n    pheatmap(data.matrix(corrected_data),\n    cluster_row=FALSE,\n    cluster_cols=FALSE,\n    treeheight_col = 0,\n    main = \"Heatmap preview\",\n    color=colorRampPalette(args$palette)(n = 299),\n    scale=\"none\",\n    border_color=FALSE,\n    show_rownames=FALSE,\n    labels_col=args$name,\n    angle_col=90,\n    filename=args$output)\n    print(paste(\"Export heatmap to \", args$output, sep=\"\"))\n  },\n  error = function(e){ \n      print(\"Failed to export heatmap\")\n  }\n)\n"
      inputs:
        heatmap_file:
          type: File
          inputBinding:
            prefix: --input
            position: 1
        column_names:
          type: string[]
          inputBinding:
            prefix: --name
            position: 2
        output_name:
          type: string
          inputBinding:
            prefix: --output
            position: 3
      outputs:
        heatmap_png:
          type: File?
          outputBinding:
            glob: '*.png'
      baseCommand:
      - Rscript
      - preview.R
  preview_histogram:
    in:
      histogram_file: make_histogram/histogram_file
      column_names: alignment_name
      output_name:
        source: regions_file
        valueFrom: $(get_root(self.basename)+"_histogram")
    out:
    - histogram_png
    - histogram_tsv
    run:
      cwlVersion: v1.0
      class: CommandLineTool
      requirements:
      - class: DockerRequirement
        dockerPull: biowardrobe2/hopach:v0.0.7
      - class: InitialWorkDirRequirement
        listing:
        - entryname: preview.R
          entry: "#!/usr/bin/env Rscript\noptions(warn=-1)\noptions(\"width\"=300)\nsuppressMessages(library(argparse))\nsuppressMessages(library(ggplot2))\nsuppressMessages(library(reshape2))\nparser <- ArgumentParser(description='Heatmap')\nparser$add_argument(\"--input\",           help='Input TSV file', type=\"character\", required=\"True\")\nparser$add_argument(\"--name\",            help='Column aliases, the order corresponds to the every thierd column in --input. First is excluded', type=\"character\", required=\"True\", nargs='+')\nparser$add_argument(\"--output\",          help='Output prefix. Default: histogram', type=\"character\", default=\"./histogram\")\nargs <- parser$parse_args(commandArgs(trailingOnly = TRUE))\nraw_data <- read.table(args$input, sep=\"\\t\", header=TRUE, stringsAsFactors=FALSE)\nselected_columns = append(seq(2, ncol(raw_data), 3), 1, after=0)\nraw_data = raw_data[,selected_columns]\ncolnames(raw_data) = append(args$name, \"distance\", after=0)\nmelt_data <- melt(raw_data, id=\"distance\")\np = ggplot(data=melt_data,\n      aes(x=distance, y=value, colour=variable)) +\n      ggtitle(\"Average Tag Density Plot\") +\n      xlab(\"Distance from gene TSS or peak center, bp\") + \n      ylab(\"Density, tags\") +\n      labs(colour = \"Sample\") +\n      geom_line()\nggsave(paste(args$output, \"png\", sep=\".\"), plot = p)\nwrite.table(raw_data,\n            file=paste(args$output, \"tsv\", sep=\".\"),\n            sep=\"\\t\",\n            row.names=FALSE,\n            col.names=TRUE,\n            quote=FALSE)\n"
      inputs:
        histogram_file:
          type: File
          inputBinding:
            prefix: --input
            position: 1
        column_names:
          type: string[]
          inputBinding:
            prefix: --name
            position: 2
        output_name:
          type: string
          inputBinding:
            prefix: --output
            position: 3
      outputs:
        histogram_png:
          type: File
          outputBinding:
            glob: '*.png'
        histogram_tsv:
          type: File
          outputBinding:
            glob: '*.tsv'
      baseCommand:
      - Rscript
      - preview.R
label: Tag enrichment heatmap and density profile around regions of interest
doc: |
  Generates tag density heatmap and histogram for the centered list of features in a headerless regions file.

  - If provided regions file is a gene list with the following columns `chrom start end name score strand` set `Gene TSS` as a re-centering criteria.
  - If provided regions file is a peak list with the following columns `chrom start end name` set `Peak Center` as a re-centering criteria.

  `score` column is always ignored.
sd:version: 100

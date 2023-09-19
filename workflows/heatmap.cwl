cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
    expressionLib:
    - var get_root = function(basename) {
          return basename.split('.').slice(0,1).join('.');
      };


'sd:upstream':
  chipseq_sample:
  - "chipseq-se.cwl"
  - "chipseq-pe.cwl"
  - "trim-chipseq-se.cwl"
  - "trim-chipseq-pe.cwl"
  - "trim-atacseq-se.cwl"
  - "trim-atacseq-pe.cwl"
  filtered_experiment:
  - "filter-peaks-for-heatmap.cwl"
  - "filter-deseq-for-heatmap.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 1

  alignment_file:
    type: File[]
    format: "http://edamontology.org/format_2572"
    label: "ChIP-Seq experiment(s)"
    doc: "Array of alignment files in BAM format"
    'sd:upstreamSource': "chipseq_sample/bambai_pair"
    'sd:localLabel': true

  alignment_name:
    type: string[]
    label: "ChIP-Seq experiment(s)"
    doc: "Names for input alignment files. Order corresponds to the alignment_file"
    'sd:upstreamSource': "chipseq_sample/alias"

  regions_file:
    type: File
    format: "http://edamontology.org/format_3003"
    label: "Filter ChIP/ATAC peaks or filter DESeq genes experiment"
    doc: |
      "Regions of interest. Formatted as headerless BED file with [chrom start end name score strand] for gene list and
       [chrom start end name] for peak file. [name] should be unique, [score] is ignored"
    'sd:upstreamSource': "filtered_experiment/filtered_file"
    'sd:localLabel': true

  recentering:
    type:
      - "null"
      - type: enum
        symbols: ["Gene TSS", "Peak Center"]
    default: "Gene TSS"
    label: "Re-center regions of interest. Chose [Gene TSS] for a gene list or [Peak Center] for a peak file"
    doc: "Re-center regions of interest. Chose [Gene TSS] for a gene list or [Peak Center] for a peak file"

  fragment_size:
    type: int[]
    label: "ChIP-Seq experiment(s)"
    doc: "Array of fragment sizes for input BAM files, order corresponds to the alignment_file"
    'sd:upstreamSource': "chipseq_sample/estimated_fragment_size"

  mapped_reads_number:
    type: int[]
    label: "ChIP-Seq experiment(s)"
    doc: "Array of mapped reads numners for input BAM files, order corresponds to the alignment_file"
    'sd:upstreamSource': "chipseq_sample/mapped_reads_number"

  hist_width:
    type: int?
    default: 10000
    label: "Histogram / heatmap width, bp"
    doc: "Histogram / heatmap width, bp"
    'sd:layout':
      advanced: true

  hist_bin_size:
    type: int?
    default: 50
    label: "Histogram / heatmap bin size, bp"
    doc: "Histogram / heatmap bin size, bp"
    'sd:layout':
      advanced: true

  threads:
    type: int?
    default: 4
    label: "Number of threads"
    doc: "Number of threads for those steps that support multithreading"
    'sd:layout':
      advanced: true


outputs:

  heatmap_table:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "TSS centered heatmap as TSV"
    doc: "TSS centered heatmap as TSV"
    outputSource: make_heatmap/histogram_file
  
  heatmap_plot:
    type: File?
    format: "http://edamontology.org/format_3603"
    label: "TSS centered heatmap as PNG"
    doc: "TSS centered heatmap as PNG"
    outputSource: preview_heatmap/heatmap_png
    'sd:visualPlugins':
    - image:
        tab: 'Plots'
        Caption: 'TSS Centered Heatmap'

  histogram_table:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "TSS centered average tag density histogram as TSV"
    doc: "TSS centered average tag density histogram as TSV"
    outputSource: preview_histogram/histogram_tsv

  histogram_plot:
    type: File
    format: "http://edamontology.org/format_3603"
    label: "TSS centered average tag density histogram as PNG"
    doc: "TSS centered average tag density histogram as PNG"
    outputSource: preview_histogram/histogram_png
    'sd:visualPlugins':
    - image:
        tab: 'Plots'
        Caption: 'Average Tag Density Plot'

  recentered_regions_file:
    type: File
    format: "http://edamontology.org/format_3003"
    label: "Re-centered by [Gene TSS] or [Peak Center] regions of interest file"
    doc: "Re-centered by [Gene TSS] or [Peak Center] regions of interest file"
    outputSource: recenter_regions/output_file

steps:

  make_tag_folders:
    run: ../tools/heatmap-prepare.cwl
    in:
      bam_file: alignment_file
      output_folder: alignment_name
      fragment_size: fragment_size
      total_reads: mapped_reads_number
    out: [tag_folder]

  recenter_regions:
    run: ../tools/custom-bash.cwl
    in:
      input_file: regions_file
      param: recentering
      script:
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
    out: [output_file]

  make_heatmap:
    run: ../tools/homer-annotate-peaks-hist.cwl
    in:
      peak_file: recenter_regions/output_file
      tag_folders: make_tag_folders/tag_folder
      hist_width: hist_width
      hist_bin_size: hist_bin_size
      export_heatmap:
        default: True
      threads: threads
      histogram_filename:
        source: regions_file
        valueFrom: $(get_root(self.basename)+"_heatmap.cdt")
    out: [histogram_file]

  make_histogram:
    run: ../tools/homer-annotate-peaks-hist.cwl
    in:
      peak_file: recenter_regions/output_file
      tag_folders: make_tag_folders/tag_folder
      hist_width: hist_width
      hist_bin_size: hist_bin_size
      export_heatmap:
        default: False
      threads: threads
      histogram_filename:
        source: regions_file
        valueFrom: $(get_root(self.basename)+"_histogram.tsv")
    out: [histogram_file]

  preview_heatmap:
    in:
      heatmap_file: make_heatmap/histogram_file
      column_names: alignment_name
      output_name:
        source: regions_file
        valueFrom: $(get_root(self.basename)+"_heatmap.png")
    out: [heatmap_png]
    run:
      cwlVersion: v1.0
      class: CommandLineTool
      requirements:
      - class: DockerRequirement
        dockerPull: biowardrobe2/hopach:v0.0.7
      - class: InitialWorkDirRequirement
        listing:
          - entryname: preview.R
            entry: | 
              #!/usr/bin/env Rscript
              options(warn=-1)
              options("width"=300)
              suppressMessages(library(argparse))
              suppressMessages(library(RColorBrewer))
              suppressMessages(library(pheatmap))
              parser <- ArgumentParser(description='Heatmap')
              parser$add_argument("--input",           help='Input CDT file', type="character", required="True")
              parser$add_argument("--name",            help='Input aliases, the order and number corresponds to --input', type="character", required="True", nargs='+')
              parser$add_argument("--palette",         help='Palette color names. Default: black, yellow, white',         type="character", nargs='+', default=c("black", "yellow", "white"))
              parser$add_argument("--output",          help='Output prefix. Default: heatmap', type="character", default="./heatmap.png")
              args <- parser$parse_args(commandArgs(trailingOnly = TRUE))
              raw_data <- read.table(args$input, sep="\t", header=TRUE, stringsAsFactors=FALSE)
              corrected_data <- raw_data[,-1]
              rownames(corrected_data) <- raw_data[,1]
              print("Centering by mean")
              corrected_data = corrected_data - rowMeans(corrected_data)    
              print("Normalizing")
              std = sqrt(rowSums(corrected_data^2))
              corrected_data = corrected_data/std
              corrected_data = replace(corrected_data, is.na(corrected_data), 0)
              tryCatch(
                expr = {
                  pheatmap(data.matrix(corrected_data),
                  cluster_row=FALSE,
                  cluster_cols=FALSE,
                  treeheight_col = 0,
                  main = "Heatmap preview",
                  color=colorRampPalette(args$palette)(n = 299),
                  scale="none",
                  border_color=FALSE,
                  show_rownames=FALSE,
                  labels_col=args$name,
                  angle_col=90,
                  filename=args$output)
                  print(paste("Export heatmap to ", args$output, sep=""))
                },
                error = function(e){ 
                    print("Failed to export heatmap")
                }
              )
      inputs:
        heatmap_file:
          type: File
          inputBinding:
            prefix: "--input"
            position: 1
        column_names:
          type: string[]
          inputBinding:
            prefix: "--name"
            position: 2
        output_name:
          type: string
          inputBinding:
            prefix: "--output"
            position: 3
      outputs:
        heatmap_png:
          type: File?
          outputBinding:
            glob: "*.png"
      baseCommand: ["Rscript", "preview.R"]

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
            entry: |
              #!/usr/bin/env Rscript
              options(warn=-1)
              options("width"=300)
              suppressMessages(library(argparse))
              suppressMessages(library(ggplot2))
              suppressMessages(library(reshape2))
              parser <- ArgumentParser(description='Heatmap')
              parser$add_argument("--input",           help='Input TSV file', type="character", required="True")
              parser$add_argument("--name",            help='Column aliases, the order corresponds to the every thierd column in --input. First is excluded', type="character", required="True", nargs='+')
              parser$add_argument("--output",          help='Output prefix. Default: histogram', type="character", default="./histogram")
              args <- parser$parse_args(commandArgs(trailingOnly = TRUE))
              raw_data <- read.table(args$input, sep="\t", header=TRUE, stringsAsFactors=FALSE)
              selected_columns = append(seq(2, ncol(raw_data), 3), 1, after=0)
              raw_data = raw_data[,selected_columns]
              colnames(raw_data) = append(args$name, "distance", after=0)
              melt_data <- melt(raw_data, id="distance")
              p = ggplot(data=melt_data,
                    aes(x=distance, y=value, colour=variable)) +
                    ggtitle("Average Tag Density Plot") +
                    xlab("Distance from gene TSS or peak center, bp") + 
                    ylab("Density, tags") +
                    labs(colour = "Sample") +
                    geom_line()
              ggsave(paste(args$output, "png", sep="."), plot = p)
              write.table(raw_data,
                          file=paste(args$output, "tsv", sep="."),
                          sep="\t",
                          row.names=FALSE,
                          col.names=TRUE,
                          quote=FALSE)
      inputs:
        histogram_file:
          type: File
          inputBinding:
            prefix: "--input"
            position: 1
        column_names:
          type: string[]
          inputBinding:
            prefix: "--name"
            position: 2
        output_name:
          type: string
          inputBinding:
            prefix: "--output"
            position: 3
      outputs:
        histogram_png:
          type: File
          outputBinding:
            glob: "*.png"
        histogram_tsv:
          type: File
          outputBinding:
            glob: "*.tsv"
      baseCommand: ["Rscript", "preview.R"]     

$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "Tag density profile around regions of interest"
label: "Tag density profile around regions of interest"
s:alternateName: "Generate tag density heatmap and histogram around gene TSS or peak centers"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/heatmap.cwl
s:codeRepository: https://github.com/datirium/workflows
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
        s:email: mailto:misha.kotliar@gmail.com
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898


# doc:
#   $include: ../descriptions/heatmap.md


doc: |
  Generates tag density heatmap and histogram for the centered list of features in a headerless regions file.

  - If provided regions file is a gene list with the following columns `chrom start end name score strand` set `Gene TSS` as a re-centering criteria.
  - If provided regions file is a peak list with the following columns `chrom start end name` set `Peak Center` as a re-centering criteria.
  
  `score` column is always ignored.

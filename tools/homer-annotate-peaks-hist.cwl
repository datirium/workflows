cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/homer:v0.0.2
inputs:
  peak_file:
    type: File
    doc: Homer generated peak file or BED
  tag_folders:
    type:
    - Directory
    - Directory[]
    inputBinding:
      position: 7
      prefix: -d
    doc: Tag folders from homer-make-tag-directory tool
  hist_width:
    type:
    - int
    - string
    inputBinding:
      position: 8
      prefix: -size
    doc: |
      Possible values:
        "#" - performs analysis on the "#" bp surrounding the peak centers
        "#,#" - performs analysis from "#" to "#" relative to peak center
        "given" - set size to actual coordinates in peak/BED file
  hist_bin_size:
    type: int
    inputBinding:
      position: 9
      prefix: -hist
    doc: |
      Bin size, bp. If hist_width is "given" or skipped, this
      parameter will set the number of bins to divide each region into
  export_heatmap:
    type: boolean?
    inputBinding:
      position: 10
      prefix: -ghist
    doc: |
      Generate heatmap. Instead of averaging all of the data
      from each peak, keep data from each peak separate
  norm_fpkm:
    type: boolean?
    inputBinding:
      position: 11
      prefix: -fpkm
    doc: |
      Normalize read counts to million reads or fragments per kilobase mapped
  norm_raw:
    type: boolean?
    inputBinding:
      position: 12
      prefix: -raw
    doc: |
      Do not adjust the tag counts based on total tags sequenced.
      By default all tag counts will be normalized to norm_tag_count
  norm_tag_count:
    type: int?
    inputBinding:
      position: 13
      prefix: -norm
    doc: |
      Normalize tags to this tag count, default=1e7, 0=average tag count in all directories
  norm_fragment_size:
    type: int?
    inputBinding:
      position: 14
      prefix: -normLength
    doc: |
      Fragment length to normlize to for experiments with different lens. Default: 100bp
  strand:
    type: string?
    inputBinding:
      position: 15
      prefix: -strand
    doc: |
      Count tags on specific strands relative to peak. Default: both
      Possible values: +|-
  threads:
    type: int?
    inputBinding:
      position: 16
      prefix: -cpu
    doc: |
      Set the number of threads. Default: 1
  histogram_filename:
    type: string
    doc: Output histogram's filename
outputs:
  histogram_file:
    type: stdout
    doc: Output histogram file
stdout: ${return inputs.histogram_filename;}
baseCommand:
- annotatePeaks.pl
arguments:
- valueFrom: $(inputs.peak_file)
  position: 5
- valueFrom: $("none")
  position: 6
doc: "Tool is used to produce histogram or heatmaps only. Rest of the functionality is not implemented intentionally.\nIf TSS analysis needed, input peak_file should be centered on TSS, where the 'center' of the peak in the actual TSS.\nFor example:\n  1\tchr4\t978796\t978796\t-\n  2\tchr4\t1052109\t1052109\t+\n  3\tchr4\t1105422\t1105422\t-\n\nSkipped arguments:\n\n  Related to peaks annotation:\n    -organism\n    -gtf\n    -gff\n    -gff3\n    -gid\n    -ann\n    -mask\n    -p\n    -pdist\n    -pcount\n    -vcf\n    -editDistance\n    -individuals\n    -gene\n    -go\n    -genomeOntology\n    -ratio\n    -rlog\n    -vst\n    -CpG\n    -nfr\n    -nfrSize\n    -gwasCatalog\n    -map\n    -noann\n\n  Related to tss/tts/rna modes:\n    tss\n    tts\n    rna\n    -list\n    -cTSS\n\n  Related to motifs:\n    -m\n    -mscore\n    -nmotifs\n    -mdist\n    -mfasta\n    -fm\n    -rmrevopp\n    -matrix\n    -mbed\n    -mlogic\n    -norevopp\n\n  Related to peak centering:\n    -center\n    -mirror\n    -multi\n\n  Related to genome comparisons\n    -cmpGenome\n    -cmpLiftover\n\n  Currently not needed functionality:\n    -bedGraph\n    -wig\n    -nuc\n    -di\n    -histNorm\n    -rm\n    -log\n    -sqrt\n    -len\n    -pc\n    -noblanks\n    -homer1\n    -homer2\n"
label: homer-annotate-peaks-hist

cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/intervene:v0.0.1
inputs:
  figure_format:
    type:
    - 'null'
    - type: enum
      symbols:
      - pdf
      - svg
      - ps
      - tiff
      - png
    inputBinding:
      position: 5
      prefix: --figtype
    doc: |
      Format to export diagram figure
  intervals_colors:
    type:
    - 'null'
    - string[]
    doc: |
      Comma-separated list of matplotlib-valid colors for fill.
      E.g., –colors=r,b,k
  intervals_files:
    type: File[]
    inputBinding:
      position: 6
      prefix: --input
    doc: |
      Files with the input genomic regions in (BED/GTF/GFF) format
  intervals_aliases:
    type:
    - 'null'
    - string[]
    inputBinding:
      position: 7
      prefix: --names
      itemSeparator: ','
      valueFrom: |
        ${
          if (self !== null) {
            let cleaned = [];
            self.forEach(function (s, i) {
              cleaned.push(s.replace(/\t|\s|\[|\]|\>|\</g, "_"));
            });
            return cleaned;
          };
          return self;
        }
    doc: |
      Comma-separated list of names as labels for input files.
      For example: --names=A,B,C,D,E,F
      If it is not set file names will be used as labels.
  overlap_threshold:
    type: int?
    inputBinding:
      position: 8
      prefix: --overlap-thresh
      valueFrom: |
        ${
          if (inputs.diagram_type == "pairwise") {
            return null;
          };
          return self;
        }
    doc: |
      Minimum number of overlapped regions to save them into bed/txt file
      Default: 1
  diagram_type:
    type:
    - 'null'
    - type: enum
      symbols:
      - venn
      - upset
      - pairwise
    default: venn
    doc: |
      Diagram type to report
outputs:
  overlapped_intervals_files:
    type:
    - 'null'
    - File[]
    outputBinding:
      glob: ./results/sets/*
    doc: Overlapped intervals files
  overlapped_combinations:
    type: File?
    outputBinding:
      glob: ./results/*combinations.txt
    doc: Overlapped combinations file
  overlapped_fraction_matrix:
    type: File?
    outputBinding:
      glob: ./results/*frac_matrix.txt
    doc: Overlapped fraction matrix file
  overlapped_plot:
    type: File
    outputBinding:
      glob: |
        ${
          if (inputs.figure_format !== null) {
            return "./results/intervene*." + inputs.figure_format;
          };
          return "./results/intervene*.pdf";
        }
    doc: Diagram of the overlapped intervals
  stdout_log:
    type: stdout
  stderr_log:
    type: stderr
baseCommand:
- intervene
arguments:
- valueFrom: $(inputs.diagram_type)
- --type
- genomic
- valueFrom: $(inputs.diagram_type == "pairwise"?null:"--save-overlaps")
- valueFrom: $(inputs.diagram_type == "pairwise"?"--diagonal":null)
- valueFrom: $(inputs.diagram_type == "pairwise"?["--htype", "color"]:null)
- valueFrom: |
    ${
      if (inputs.diagram_type == "venn" && inputs.intervals_colors !== null) {
        return inputs.intervals_colors;
      }
      return null;
    }
  prefix: --colors
  itemSeparator: ','
- --project
- intervene
- -o
- results
stdout: intervene_stdout.log
stderr: intervene_stderr.log
doc: |
  Intervene: a tool for intersection and visualization of multiple genomic regions
  ================================================================================

  Hardcoded parameters:
  --type genomic       use genomic regions as inputs
  --save-overlaps      save overlapping regions/names for all the combinations as bed/txt
  --project intervene  for the guatanteed access to results location
  -o results           for the guatanteed access to results location

  Skipped parameters:
  --filenames as it seems not logical to have it when default
              for --names already uses file names

  When run in "pairwise" mode:
  --save-overlaps and --overlap-thresh set to null
  --htype and --diagonal added

  When run in "venn" mode:
  --colors aren't ignored
label: intervene

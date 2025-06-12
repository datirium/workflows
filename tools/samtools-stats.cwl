cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() { if (inputs.output_filename == ""){ var root = inputs.bambai_pair.basename.split('.').slice(0,-1).join('.'); return (root == "")?inputs.bambai_pair.basename+".log":root+".log"; } else { return inputs.output_filename; } };
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/samtools:v1.4
inputs:
  script:
    type: string?
    default: |
      #!/bin/bash
      samtools stats $0 > all_statistics.log
      cat all_statistics.log | grep SN | cut -f 2- > sn_section.log
      echo -e "insert size\tpairs total\tinward oriented pairs\toutward oriented pairs\tother pairs" > is_section.tsv
      cat all_statistics.log | grep ^IS | cut -f 2- | head -n 1501 >> is_section.tsv
      cat all_statistics.log
    inputBinding:
      position: 4
    doc: samtools stats with optional filtering criteria passes as $1
  bambai_pair:
    type: File
    inputBinding:
      position: 5
    secondaryFiles:
    - .bai
    doc: Coordinate sorted BAM alignment and index BAI files
  output_filename:
    type: string?
    default: ''
    doc: Output file name
outputs:
  log_file:
    type: File
    outputBinding:
      glob: $(default_output_filename())
    doc: BAM file statistics
  ext_is_section:
    type: File
    outputBinding:
      glob: is_section.tsv
    doc: BAM file statistics (IS section)
  raw_total_sequences:
    type: int
    outputBinding:
      loadContents: true
      glob: sn_section.log
      outputEval: |
        ${
          var s = self[0].contents.substring(self[0].contents.indexOf("raw total sequences"));
          return parseInt(s.substring(s.indexOf("raw total sequences")+21, s.indexOf("\n")));
        }
    doc: Raw Total Sequences
  reads_mapped:
    type: int
    outputBinding:
      loadContents: true
      glob: sn_section.log
      outputEval: |
        ${
          var s = self[0].contents.substring(self[0].contents.indexOf("reads mapped"));
          return parseInt(s.substring(s.indexOf("reads mapped")+14, s.indexOf("\n")));
        }
    doc: Reads Mapped
  average_length:
    type: int
    outputBinding:
      loadContents: true
      glob: sn_section.log
      outputEval: |
        ${
          var s = self[0].contents.substring(self[0].contents.indexOf("average length"));
          return parseInt(s.substring(s.indexOf("average length")+16, s.indexOf("\n")));
        }
    doc: Reads Average Length
baseCommand:
- bash
- -c
stdout: $(default_output_filename())
doc: |
  Generates statistics for the input BAM file.
label: samtools-stats

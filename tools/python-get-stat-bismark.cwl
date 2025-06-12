cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: biowardrobe2/scidap:v0.0.3
inputs:
  script:
    type: string?
    default: |
      #!/usr/bin/env python
      import sys, re
      TOTAL, UNIQUE, UNMAPPED, MULTIMAPPED, DISCARDED, CPG, CHG, CHH, ECPG, ECHG, ECHH = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
      with open(sys.argv[1], 'r') as alignment_report:
        for line in alignment_report:
          line = line.strip()
          if "in total:" in line:
            TOTAL = line.split("\t")[1]
          if 'unique best hit' in line:
            UNIQUE = line.split("\t")[1]
          if 'under any condition:' in line:
            UNMAPPED = line.split("\t")[1]
          if 'not map uniquely:' in line:
            MULTIMAPPED = line.split("\t")[1]
          if 'could not be extracted:' in line:
            DISCARDED = line.split("\t")[1]
          if 'C methylated in CpG context:' in line:
            CPG = line.split("\t")[1].strip("%")
          if 'C methylated in CHG context:' in line:
            CHG = line.split("\t")[1].strip("%")
          if 'C methylated in CHH context:' in line:
            CHH = line.split("\t")[1].strip("%")
      with open(sys.argv[2], 'r') as splitting_report:
        for line in splitting_report:
          line = line.strip()
          if 'C methylated in CpG context:' in line:
            ECPG = line.split("\t")[1].strip("%")
          if 'C methylated in CHG context:' in line:
            ECHG = line.split("\t")[1].strip("%")
          if 'C methylated in CHH context:' in line:
            ECHH = line.split("\t")[1].strip("%")
      print "Sequences analysed in total\tNumber of alignments with a unique best hit from the different alignments\tSequences with no alignments under any condition\tSequences did not map uniquely\tSequences which were discarded because genomic sequence could not be extracted\tC methylated in CpG context\tC methylated in CHG context\tC methylated in CHH context\tC methylated in CpG context after extraction\tC methylated in CHG context after extraction\tC methylated in CHH context after extraction"
      print TOTAL + "\t" + UNIQUE + "\t" + UNMAPPED + "\t" + MULTIMAPPED + "\t" + DISCARDED + "\t" + CPG + "\t" + CHG + "\t" + CHH + "\t" + ECPG + "\t" + ECHG + "\t" + ECHH
    inputBinding:
      position: 1
    doc: Python script to get TOTAL, ALIGNED, SUPPRESSED, USED values from log files
  alignment_report:
    type: File
    inputBinding:
      position: 6
    doc: Bismark generated alignment and methylation summary report
  splitting_report:
    type: File
    inputBinding:
      position: 7
    doc: Bismark generated splitting report
outputs:
  collected_report_formatted:
    type: stdout
    doc: Combined Bismark alignment and splitting reports
  mapped_reads:
    type: int
    outputBinding:
      loadContents: true
      glob: collected_report_formatted.tsv
      outputEval: $(parseInt(self[0].contents.split('\n')[1].split('\t')[1]))
stdout: collected_report_formatted.tsv
baseCommand:
- python
- -c
doc: |
  Tool refactores Bismark alignment report to be displayed as Pie Chart
label: python-get-stat-bismark

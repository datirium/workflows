cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: DockerRequirement
  dockerPull: biowardrobe2/scidap:v0.0.3


inputs:

  script:
    type: string?
    default: |
      #!/usr/bin/env python
      import sys, re
      TOTAL, UNIQUE, UNMAPPED, MULTIMAPPED, DISCARDED = 0, 0, 0, 0, 0
      with open(sys.argv[1], 'r') as alignment_report:
        for line in alignment_report:
          if "Sequences analysed in total:" in line:
            TOTAL = int(line.split("\t")[1])
          if 'different alignments:' in line:
            UNIQUE = int(line.split("\t")[1])
          if 'under any condition:' in line:
            UNMAPPED = int(line.split("\t")[1])
          if 'not map uniquely:' in line:
            MULTIMAPPED = int(line.split("\t")[1])
          if 'could not be extracted:' in line:
            DISCARDED = int(line.split("\t")[1])
      print "Sequences analysed in total\tNumber of alignments with a unique best hit from the different alignments\tSequences with no alignments under any condition\tSequences did not map uniquely\tSequences which were discarded because genomic sequence could not be extracted"
      print str(TOTAL) + "\t" + str(UNIQUE) + "\t" + str(UNMAPPED) + "\t" + str(MULTIMAPPED) + "\t" + str(DISCARDED)
    inputBinding:
      position: 1
    doc: "Python script to get TOTAL, ALIGNED, SUPPRESSED, USED values from log files"

  alignment_report:
    type: File
    inputBinding:
      position: 6
    doc: "Bismark generated alignment and methylation summary report"


outputs:

  alignment_report_formatted:
    type: stdout
    doc: "Refactored Bismark alignment report to be visualized as Pie Chart"

stdout: "alignment_report_formatted.tsv"
baseCommand: [python, '-c']


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "python-get-stat-bismark"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/python-get-stat-bismark.cwl
s:codeRepository: https://github.com/Barski-lab/workflows
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

doc: |
  Tool refactores Bismark alignment report to be displayed as Pie Chart

s:about: |
  Runs python code from the embeded script

cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var get_output_filename = function() {
        if (inputs.output_filename) {
          return inputs.output_filename;
        }
        return inputs.star_log.location.split('/').slice(-1)[0].replace(/\.*Log\.final\.out$/i,'')+".stat";
    };
  - var get_formatted_output_filename = function() {
      if (inputs.formatted_output_filename) {
            return inputs.formatted_output_filename;
         }
      return inputs.star_log.location.split('/').slice(-1)[0].replace(/\.*Log\.final\.out$/i,'')+"_stats.tsv";
    };


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/scidap:v0.0.3


inputs:

  script:
    type: string?
    default: |
      # !/usr/bin/env python
      import sys, re
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
      with open(sys.argv[3], 'r') as rpkm_isoforms_file:
          header = None
          key_index = None
          for line in rpkm_isoforms_file:
              if line == '':
                  continue
              if re.match ('.*RefseqId.*|.*GeneId.*|.*Chrom.*|.*TotalReads.*', line) and header is None:
                  header = line.split(',')
                  key_index = header.index('TotalReads')
                  continue
              line_splitted = line.split(',')
              USED += int(line_splitted[key_index])
      if len(sys.argv) > 6 and sys.argv[6] == "--pair":
          USED = USED/2
      with open(sys.argv[4], 'w') as fo:
          fo.write(str(TOTAL) + " " + str(ALIGNED) + " " + str(RIBO) + " " + str(MULTIMAPPED) + " " + str(USED))
      with open(sys.argv[5], 'w') as fof:
          fof.write("Tags total\tTranscriptome\tMulti-mapped\tOutside annotation\tUnmapped\tRibosomal contamination\n")
          fof.write(str(TOTAL) + "\t" + str(USED) + "\t" + str(MULTIMAPPED) + "\t" + str(ALIGNED-USED) + "\t" + str(TOTAL-ALIGNED-MULTIMAPPED) + "\t" + str(RIBO) + "\n")
    inputBinding:
      position: 5
    doc: "Python script to get TOTAL, ALIGNED, RIBO, MULTIMAPPED, USED values from log files"

  star_log:
    type: File
    inputBinding:
      position: 6
    doc: "Log file from STAR (Log.final.out)"

  bowtie_log:
    type: File
    inputBinding:
      position: 7
    doc: "Log file from Bowtie"

  rpkm_isoforms:
    type: File
    inputBinding:
      position: 8
    doc: "RPKM calcualted by GEEP, grouped by isoforms"
  
  output_filename:
    type: string?
    inputBinding:
      position: 9
      valueFrom: $(get_output_filename())
    default: ""
    doc: "Name for generated output file"
  
  formatted_output_filename:
    type: string?
    inputBinding:
      position: 10
      valueFrom: $(get_formatted_output_filename())
    default: ""
    doc: "Name for generated formatted output file"

  pair_end:
    type: boolean?
    inputBinding:
      position: 11
      prefix: --pair
    doc: "If true, USED values is divided on 2"


outputs:

  output_file:
    type: File
    outputBinding:
      glob: $(get_output_filename())

  formatted_output_file:
    type: File
    outputBinding:
      glob: $(get_formatted_output_filename())

  total_reads:
    type: int
    outputBinding:
      loadContents: true
      glob: $(get_output_filename())
      outputEval: $(parseInt(self[0].contents.split(' ')[0]))

  mapped_reads:
    type: int
    outputBinding:
      loadContents: true
      glob: $(get_output_filename())
      outputEval: $(parseInt(self[0].contents.split(' ')[1]))

  ribo_reads:
    type: int
    outputBinding:
      loadContents: true
      glob: $(get_output_filename())
      outputEval: $(parseInt(self[0].contents.split(' ')[2]))

  multimapped_reads:
    type: int
    outputBinding:
      loadContents: true
      glob: $(get_output_filename())
      outputEval: $(parseInt(self[0].contents.split(' ')[3]))

  used_reads:
    type: int
    outputBinding:
      loadContents: true
      glob: $(get_output_filename())
      outputEval: $(parseInt(self[0].contents.split(' ')[4]))


baseCommand: [python, '-c']


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "python-get-stat-rnaseq"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/python-get-stat-rnaseq.cwl
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
  Tool processes and combines log files generated by STAR/Bowtie aligners and GEEP rpkm results file.

  `get_output_filename` function returns output filename equal to `output_filename` (if input is provided) or
  generated on the base of STAR log basename with `.stat` extension.

  `get_formatted_output_filename` function returns output filename equal to `formatted_output_filename` (if input is provided) or
  generated on the base of STAR log basename with `_stats.tsv` extension.

s:about: |
  Runs python code from the `script` input

#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: ./metadata/envvar-global.yml
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/samtools:v1.4
  dockerFile: >
    $import: ./dockerfiles/samtools-Dockerfile

inputs:

  bash_script:
    type: string?
    default: |
      #!/usr/bin/env python
      import sys
      import os
      import subprocess as sub
      import magic
      import json
      import uuid

      def set_fail (conf_data, key):
          for step_name, value in conf_data['steps'].iteritems():
              if key in value[1]:
                  value[0] = "false"

      def compare (control_item, check_item):
          if not os.path.isfile(control_item): raise ValueError("missing file   {0}".format(control_item))
          if not os.path.isfile(check_item): raise ValueError("missing file   {0}".format(check_item))
          if magic.from_file(control_item, mime=True) == 'text/plain' and magic.from_file(check_item, mime=True) == 'text/plain':
              sub.check_output('diff -I "#.*" -B -q ' + value + ' ' + check_dict[key], shell=True)
          elif os.path.getsize(control_item) != os.path.getsize(check_item):
              raise ValueError('size comparison')

      def replace_bam_to_idxstats(files_dict):
          for key, value in files_dict.iteritems():
              if os.path.splitext(key)[1].lower() == '.bam'.lower():
                  print "     Found BAM file   {0}   ({1})".format(key,value)
                  try:
                      new_value = str(uuid.uuid4()) + '.idxstats'
                      sub.check_output('samtools idxstats %s > %s' % (value, new_value), shell=True)
                      files_dict[key] = os.path.abspath(new_value)
                      print "     Replaced   {0} --> {1}".format(value, files_dict[key])
                  except Exception, e:
                      print '    ',str(e)

      conf =  '{ \
                    "overall": "true", \
                    "steps": { \
                      "fastx_quality_stats": ["true",            ["SRR1198790.fastq.fastxstat"]], \
                      "bowtie_aligner": ["true",                 ["SRR1198790.sam.log", \
                                                                   "SRR1198790.sorted.bam", \
                                                                   "SRR1198790.sam.stat"]], \
                      "samtools_sort_index": ["true",            ["SRR1198790.sorted.bam", \
                                                                   "SRR1198790.sorted.bam.bai", \
                                                                   "SRR1198790.sam.stat"]], \
                      "samtools_rmdup": ["true",                 ["SRR1198790.sorted.bam", \
                                                                   "SRR1198790.sorted.bam.bai", \
                                                                   "SRR1198790.sam.stat"]], \
                      "samtools_sort_index_after_rmdup": ["true",["SRR1198790.sorted.bam", \
                                                                   "SRR1198790.sorted.bam.bai"]], \
                      "macs2_callpeak": ["true",                 ["SRR1198790.sorted_model.r", \
                                                                   "SRR1198790.sorted_peaks.narrowPeak", \
                                                                   "SRR1198790.sorted_peaks.xls", \
                                                                   "SRR1198790.sorted_summits.bed"]], \
                      "macs_island_count": ["true",              ["SRR1198790.sorted_peaks.xls"]], \
                      "macs2_callpeak_forced": ["true",          ["SRR1198790.sorted_model.r", \
                                                                   "SRR1198790.sorted_peaks.narrowPeak", \
                                                                   "SRR1198790.sorted_peaks.xls", \
                                                                   "SRR1198790.sorted_summits.bed"]], \
                      "bamtools_stats": ["true",                 ["SRR1198790.sorted.sorted.bigwig"]], \
                      "bam_to_bigwig": ["true",                  ["SRR1198790.sorted.sorted.bigwig"]] \
                     } \
                 }'

      conf_obj = json.loads(conf)

      control_dict = {key:os.path.abspath(os.path.join(sys.argv[1], key)) for key in os.listdir(sys.argv[1])}
      check_dict = {key:os.path.abspath(os.path.join(sys.argv[2], key)) for key in os.listdir(sys.argv[2])}
      print "LOAD control directory"
      replace_bam_to_idxstats (control_dict)
      print "LOAD output directory"
      replace_bam_to_idxstats (check_dict)

      print 'PROCESS'
      for key,value in control_dict.iteritems():
          print "   - {0}   ({1})".format(key,value)
          if key not in check_dict:
              print "     Fail: missing output file   {0}".format(key)
              set_fail (conf_obj, key)
              continue
          try:
              compare (value, check_dict[key])
              print "     Ok"
              continue
          except Exception, e:
              print '     Fail:', str(e)
              set_fail (conf_obj, key)

      for step_name in conf_obj['steps']:
          state = conf_obj['steps'][step_name][0]
          conf_obj['steps'][step_name] = state
          if state == 'false':
              conf_obj['overall'] = 'false'

      with open('results.json', 'w') as resutls_file:
          resutls_file.write(json.dumps(conf_obj, indent = 4))
    inputBinding:
      position: 5
    doc: |
      Bash function to run samtools rmdup with all input parameters or skip it if trigger is false

  control_dir:
    type: Directory
    inputBinding:
      position: 6
    doc: |
      Folder to include original files to check against to.

  check_dir:
    type: Directory
    inputBinding:
      position: 7
    doc: |
      Folder to include newly generated files.

outputs:
  results:
    type: File?
    outputBinding:
      glob: results.json
    doc: |
      Main results file
        {
          “overall” : “true|false”,
          “steps”: {
            “STEP_NAME” : “true|false”,
            ....
           }
         }
  log:
    type: File?
    outputBinding:
      glob: log.txt

baseCommand: [python, '-c']

arguments:
  - valueFrom: |
      ${
          return " > log.txt 2>&1";
      }
    position: 100000
    shellQuote: false

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "run-dna-single-end-validate"
s:downloadUrl: https://raw.githubusercontent.com/SciDAP/workflows/master/workflows/scidap/run-dna-single-end-validate.cwl
s:codeRepository: https://github.com/SciDAP/workflows
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
        s:email: mailto:michael.kotliar@cchmc.org
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898

doc: |
  Current workflow is used to validate output of run-dna-single-end-validate.cwl

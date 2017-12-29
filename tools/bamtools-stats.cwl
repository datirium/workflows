#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: ./metadata/envvar-global.yml
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/bamtools:v2.4.1


inputs:

  input_files:
    type:
      - File
      - type: array
        items: File
        inputBinding:
          prefix: -in
    inputBinding:
      position: 2
      valueFrom: |
        ${
          if ( Object.prototype.toString.call(inputs.input_files) === '[object Array]'){
            return null;
          } else {
            return ["-in", inputs.input_files.path];
          }
        }
    doc: |
      the input BAM file[s]
      NOTE: If cwl fix a bug https://github.com/common-workflow-language/common-workflow-language/issues/330
      we'll be able to use MultipleInputFeatureRequirement for single-item array and it will work
      even without additional File type and valueFrom field

outputs:
  stats_log:
    type: File
    outputBinding:
      glob: "stats.log"

  totalReads:
    type: int
    outputBinding:
      loadContents: true
      glob: "stats.log"
      outputEval: |
        ${
          var s = self[0].contents.replace(/ /g,'').replace(/ *\([^)]*\) */g,'');
          var totalReads = parseInt(s.substring ( s.indexOf("Totalreads")+11, s.indexOf("\t", (s.indexOf("Totalreads")))  ));
          return totalReads;
        }
  mappedreads:
    type: int
    outputBinding:
      loadContents: true
      glob: "stats.log"
      outputEval: |
        ${
          var s = self[0].contents.replace(/ /g,'').replace(/ *\([^)]*\) */g,'');
          var mappedreads = parseInt(s.substring ( s.indexOf("Mappedreads")+12, s.indexOf("\t", (s.indexOf("Mappedreads")))  ));
          return mappedreads;
        }
  forwardstrand:
    type: int
    outputBinding:
      loadContents: true
      glob: "stats.log"
      outputEval: |
        ${
          var s = self[0].contents.replace(/ /g,'').replace(/ *\([^)]*\) */g,'');
          var forwardstrand = parseInt(s.substring ( s.indexOf("Forwardstrand")+14, s.indexOf("\t", (s.indexOf("Forwardstrand")))  ));
          return forwardstrand;
        }
  reversestrand:
    type: int
    outputBinding:
      loadContents: true
      glob: "stats.log"
      outputEval: |
        ${
          var s = self[0].contents.replace(/ /g,'').replace(/ *\([^)]*\) */g,'');
          var reversestrand = parseInt(s.substring ( s.indexOf("Reversestrand")+14, s.indexOf("\t", (s.indexOf("Reversestrand")))  ));
          return reversestrand;
        }
  failedQC:
    type: int
    outputBinding:
      loadContents: true
      glob: "stats.log"
      outputEval: |
        ${
          var s = self[0].contents.replace(/ /g,'').replace(/ *\([^)]*\) */g,'');
          var failedQC = parseInt(s.substring ( s.indexOf("FailedQC")+9, s.indexOf("\t", (s.indexOf("FailedQC")))  ));
          return failedQC;
        }
  duplicates:
    type: int
    outputBinding:
      loadContents: true
      glob: "stats.log"
      outputEval: |
        ${
          var s = self[0].contents.replace(/ /g,'').replace(/ *\([^)]*\) */g,'');
          var duplicates = parseInt(s.substring ( s.indexOf("Duplicates")+11, s.indexOf("\t", (s.indexOf("Duplicates")))  ));
          return duplicates;
        }
  pairedendreads:
    type: int
    outputBinding:
      loadContents: true
      glob: "stats.log"
      outputEval: |
        ${
          var s = self[0].contents.replace(/ /g,'').replace(/ *\([^)]*\) */g,'');
          var pairedendreads = parseInt(s.substring ( s.indexOf("Paired-endreads")+16, s.indexOf("\t", (s.indexOf("Paired-endreads")))  ));
          return pairedendreads;
        }

baseCommand: [bamtools]
arguments:
  - stats
  - valueFrom: $('> ' + 'stats.log')
    position: 1000
    shellQuote: false

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: ./metadata/bamtools-metadata.yaml

s:name: "bamtools-stats"
s:downloadUrl: https://raw.githubusercontent.com/SciDAP/workflows/master/tools/bamtools-stats.cwl
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
  Tool is used to calculate general alignment statistics from the input BAM file

s:about: >
  Usage: bamtools stats [-in <filename> -in <filename> ... | -list <filelist>] [statsOptions]

  Input & Output:
    -in <BAM filename>                the input BAM file [stdin]
    -list <filename>                  the input BAM file list, one
                                      line per file

  Additional Stats:
    -insert                           summarize insert size data

  Help:
    --help, -h                        shows this help text

#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: ./metadata/envvar-global.yml
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement

hints:
- class: DockerRequirement
  dockerPull: scidap/bamtools:v2.4.1
  dockerFile: >
    $import: ./dockerfiles/bamtools-Dockerfile

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
    type: double
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
    type: double
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
    type: double
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
    type: double
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
    type: double
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
    type: double
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
    type: double
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

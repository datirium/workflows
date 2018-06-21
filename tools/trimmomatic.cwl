cwlVersion: v1.0
class: CommandLineTool


requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: InitialWorkDirRequirement
    listing: |
      ${
        var listing = []
        if (inputs.fastq_file_upstream){
          listing.push(inputs.fastq_file_upstream);
        }
        if (inputs.fastq_file_downstream){
          listing.push(inputs.fastq_file_downstream);
        }
        return listing;
      }

  - class: InlineJavascriptRequirement
    expressionLib:
    - var default_output_filename = function(input_file, ext) {
          let root = input_file.basename.split('.').slice(0,-1).join('.');
          return (root == "")?inputs.input_file.basename+ext:root+ext;
      };

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/trimmomatic:v0.35


inputs:

  bash_script:
    type: string?
    default: |
      #!/bin/bash
      if [ "$0" = "true" ]
      then
        echo "Run trimmomatic"
        java "${@:1}"
      else
        echo "Skip run trimmomatic"
      fi
    inputBinding:
      position: 1
    doc: |
      Bash function to run trimmomatic with all input parameters or skip it if trigger is false

  trigger:
    type: boolean?
    default: true
    inputBinding:
      position: 5
      valueFrom: |
        ${ return self ? "true" : "false" }
    doc: |
      If true - run trimmomatic, if false - return input fastq files, previously staged into output directory.
      Use valueFrom to return string instead of boolean, because if return boolean False, argument is not printed

  java_opts:
    type: string?
    inputBinding:
      position: 6
      shellQuote: false
    doc: JVM arguments should be a quoted, space separated list (e.g. "-Xms128m -Xmx512m")

  trimmomatic_jar_path:
    type: string
    default: '/usr/share/java/trimmomatic.jar'
    inputBinding:
      position: 7
      prefix: -jar

  lib_type:
    type:
      - type: enum
        name: "format"
        symbols: ['SE','PE']
    inputBinding:
      position: 8
    doc: |
      SE|PE
      Single End (SE) or Paired End (PE) mode

  phred:
    type:
      - "null"
      - type: enum
        name: "phred"
        symbols: ['33','64']
    inputBinding:
      prefix: -phred
      separate: false
      position: 9
    doc: |
      "33"|"64"
      -phred33 ("33") or -phred64 ("64") specifies the base quality encoding. Default: -phred64

  validate_pairs:
    type: boolean?
    inputBinding:
      prefix: "-validatePairs"
      position: 10
    doc: |
      Run -validatePairs. No official documentation

  threads:
    type: int
    default: 1
    inputBinding:
      position: 11
      prefix: -threads
    doc: |
      Number of threads

  log_filename:
    type: string?
    inputBinding:
      position: 12
      prefix: -trimlog
    doc: |
      <ouptut log file name>
      Specifying a trimlog file creates a log of all read trimmings, indicating the following details:
        the read name
        the surviving sequence length
        the location of the first surviving base, aka. the amount trimmed from the start
        the location of the last surviving base in the original read
        the amount trimmed from the end
      <ouptut log file name>: filename for the generated output log file.

  fastq_file_upstream:
    type: File
    inputBinding:
      position: 13
    doc: |
      FASTQ file for input read (read R1 in Paired End mode)

  fastq_file_downstream:
    type: File?
    inputBinding:
      position: 14
    doc: |
      FASTQ file for read R2 in Paired End mode

  adapters_file:
    type: File
    doc: |
      FASTA file containing adapters, PCR sequences, etc. It is used to search
      for and remove these sequences in the input FASTQ file(s)

  illuminaclip_step_param:
    type: string
    doc: |
      <fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>:<minAdapterLength>:<keepBothReads>
      Find and remove Illumina adapters.
      REQUIRED:
      <fastaWithAdaptersEtc>: specifies the path to a fasta file containing all the adapters, PCR sequences etc.
      The naming of the various sequences within this file determines how they are used. See below.
      <seedMismatches>: specifies the maximum mismatch count which will still allow a full match to be performed
      <palindromeClipThreshold>: specifies how accurate the match between the two 'adapter ligated' reads must be
      for PE palindrome read alignment.
      <simpleClipThreshold>: specifies how accurate the match between any adapter etc. sequence must be against a read
      OPTIONAL:
      <minAdapterLength>: In addition to the alignment score, palindrome mode can verify
      that a minimum length of adapter has been detected. If unspecified, this defaults to 8 bases,
      for historical reasons. However, since palindrome mode has a very low false positive rate, this
      can be safely reduced, even down to 1, to allow shorter adapter fragments to be removed.
      <keepBothReads>: After read-though has been detected by palindrome mode, and the
      adapter sequence removed, the reverse read contains the same sequence information as the
      forward read, albeit in reverse complement. For this reason, the default behaviour is to
      entirely drop the reverse read. By specifying „true‟ for this parameter, the reverse read will
      also be retained, which may be useful e.g. if the downstream tools cannot handle a
      combination of paired and unpaired reads.

  sliding_window_step:
    type: string?
    inputBinding:
      position: 101
      prefix: 'SLIDINGWINDOW:'
      separate: false
    doc: |
      <windowSize>:<requiredQuality>
      Perform a sliding window trimming, cutting once the average quality within the window falls
      below a threshold. By considering multiple bases, a single poor quality base will not cause the
      removal of high quality data later in the read.
      <windowSize>: specifies the number of bases to average across
      <requiredQuality>: specifies the average quality required

  leading_step:
    type: int?
    inputBinding:
      position: 102
      prefix: 'LEADING:'
      separate: false
    doc: |
      <quality>
      Remove low quality bases from the beginning. As long as a base has a value below this
      threshold the base is removed and the next base will be investigated.
      <quality>: Specifies the minimum quality required to keep a base.

  trailing_step:
    type: int?
    inputBinding:
      position: 103
      prefix: 'TRAILING:'
      separate: false
    doc: |
      <quality>
      Remove low quality bases from the end. As long as a base has a value below this threshold
      the base is removed and the next base (which as trimmomatic is starting from the 3‟ prime end
      would be base preceding the just removed base) will be investigated. This approach can be
      used removing the special illumina „low quality segment‟ regions (which are marked with
      quality score of 2), but we recommend Sliding Window or MaxInfo instead
      <quality>: Specifies the minimum quality required to keep a base.

  crop_step:
    type: int?
    inputBinding:
      position: 104
      prefix: 'CROP:'
      separate: false
    doc: |
      <length>
      Removes bases regardless of quality from the end of the read, so that the read has maximally
      the specified length after this step has been performed. Steps performed after CROP might of
      course further shorten the read.
      <length>: The number of bases to keep, from the start of the read.

  headcrop_step:
    type: int?
    inputBinding:
      position: 105
      prefix: 'HEADCROP:'
      separate: false
    doc: |
      <length>
      Removes the specified number of bases, regardless of quality, from the beginning of the read.
      <length>: The number of bases to keep, from the start of the read.

  minlen_step:
    type: int?
    inputBinding:
      position: 106
      prefix: 'MINLEN:'
      separate: false
    doc: |
      <length>
      This module removes reads that fall below the specified minimal length. If required, it should
      normally be after all other processing steps. Reads removed by this step will be counted and
      included in the „dropped reads‟ count presented in the trimmomatic summary.
      <length>:  Specifies the minimum length of reads to be kept

  tophred64_step:
    type: boolean?
    inputBinding:
      position: 107
      prefix: TOPHRED64
      separate: false
    doc: |
      This (re)encodes the quality part of the FASTQ file to base 64.

  tophred33_step:
    type: boolean?
    inputBinding:
      position: 107
      prefix: TOPHRED33
      separate: false
    doc: |
      This (re)encodes the quality part of the FASTQ file to base 33.


outputs:

  upstream_trimmed_file:
    type: File
    outputBinding:
      glob: |
        ${
          return inputs.trigger?default_output_filename(inputs.fastq_file_upstream, '.trimmed.fastq'):inputs.fastq_file_upstream.basename;
        }
  upstream_trimmed_unpaired_file:
    type: File?
    outputBinding:
      glob: |
        ${
          return inputs.lib_type=="PE"?default_output_filename(inputs.fastq_file_upstream, '.unpaired.trimmed.fastq'):null;
        }
  downstream_trimmed_file:
    type: File?
    outputBinding:
      glob: |
        ${
          if (inputs.lib_type=="PE"){
            return inputs.trigger?default_output_filename(inputs.fastq_file_downstream,'.trimmed.fastq'):inputs.fastq_file_downstream.basename;
          }
          return null;
        }

  downstream_trimmed_unpaired_file:
    type: File?
    outputBinding:
      glob: |
        ${
          return inputs.lib_type=="PE"?default_output_filename(inputs.fastq_file_downstream, '.unpaired.trimmed.fastq'):null;
        }

  log_file:
    type: File?
    outputBinding:
      glob: $(inputs.log_filename)
    doc: |
      Trimmomatic Log file.


baseCommand: [bash, '-c']

arguments:
# upstream_trimmed_file
- valueFrom: $(default_output_filename(inputs.fastq_file_upstream, '.trimmed.fastq'))
  position: 15
# upstream_trimmed_unpaired_file
- valueFrom: $(inputs.lib_type=="PE"?default_output_filename(inputs.fastq_file_upstream, '.trimmed.unpaired.fastq'):null)
  position: 16
# downstream_trimmed_file
- valueFrom: $(inputs.lib_type=="PE"?default_output_filename(inputs.fastq_file_downstream, '.trimmed.fastq'):null)
  position: 17
# downstream_trimmed_unpaired_file
- valueFrom: $(inputs.lib_type=="PE"?default_output_filename(inputs.fastq_file_downstream, '.trimmed.unpaired.fastq'):null)
  position: 18
- valueFrom: $("ILLUMINACLIP:" + inputs.adapters_file.path + ":"+ inputs.illuminaclip_step_param)
  position: 100


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: ./metadata/trimmomatic-metadata.yaml

s:name: "trimmomatic"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/trimmomatic.cwl
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
  Tool runs trimmomatic with ILLUMINACLIP step by default.

  `-basein` and `-baseout` inputs are skipped.

  If set `lib_type` to `PE`, both of the inputs `fastq_file_upstream` and `fastq_file_downstream` shoul be provided.

  If input `trigger` is set to `true` or isn't set at all (`true` is used by default), run `trimmomatic` and return
  FASTQ file[s] with trimmed adapters, alongside with the uppaired reads FASTQ files (if `lib_type` is set to `PE` and
  such files are present after running `trimmomatic`)
  If input `trigger` is set to `false`, return unchanged `fastq_file_upstream` and `fastq_file_downstream`, previously
  staged into output directory.

  Before execution `baseCommand`, `fastq_file_upstream` and `fastq_file_downstream` (if provided) are staged into directory
  set as docker parameter `--workdir` (tool's output directory), using `InitialWorkDirRequirement`. They are mount
  to docker container with `ro` mode as part of `--workdir`, because all generated files will have `trimmed` suffix
  in their names, so the staged files will not be overwritten.

  Trigger logic is implemented in a bash scripts set by default as `bash_script` input. If the first argument $0 (which is
  `trigger` input) is true, run `trimmomatic` with the rest of the arguments. If $0 is not true, skip `trimmomatic` and
  return `fastq_file_upstream` and `fastq_file_downstream` (if provided) staged into output directory.

  Input `trigger` is Boolean, but returns String, because of `valueFrom` field. The `valueFrom` is used, because if `trigger`
  is false, cwl-runner doesn't append this argument at all to the the `baseCommand` - new feature of CWL v1.0.2. Alternatively,
  `prefix` field could be used, but it causes changing in script logic.

  `default_output_name` function is used for generating output filename based on `input_file.basename` and provided
  extension.

s:about: |
  Usage:
         PE [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-quiet] [-validatePairs] [-basein <inputBase> | <inputFile1> <inputFile2>] [-baseout <outputBase> | <outputFile1P> <outputFile1U> <outputFile2P> <outputFile2U>] <trimmer1>...
     or:
         SE [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-quiet] <inputFile> <outputFile> <trimmer1>...

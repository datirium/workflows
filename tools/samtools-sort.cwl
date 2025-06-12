cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: ShellCommandRequirement
- class: InitialWorkDirRequirement
  listing: |
    ${
      return  [
                {
                  "entry": inputs.sort_input,
                  "entryname": inputs.sort_input.basename,
                  "writable": true
                }
              ]
    }
- class: InlineJavascriptRequirement
  expressionLib:
  - var ext = function() { if (inputs.out_format && inputs.out_format == 'SAM'){ return '.sam'; } else if (inputs.out_format && inputs.out_format == 'CRAM') { return '.cram'; } else { return '.bam'; } };
  - var default_output_name = function() { if (inputs.trigger == true){ return inputs.sort_input.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.') + ext(); } else { return inputs.sort_input.location.split('/').slice(-1)[0]; } };
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/samtools:v1.4
inputs:
  bash_script_sort:
    type: string?
    default: |
      #!/bin/bash
      if [ "$0" = "true" ]
      then
        samtools sort "${@:1}"
      else
        echo "Skip samtools sort " ${@:1}
      fi
    inputBinding:
      position: 5
    doc: |
      Bash function to run samtools sort with all input parameters or skip it if trigger is false
  trigger:
    type: boolean?
    default: true
    inputBinding:
      position: 6
      valueFrom: |
        ${ return self ? "true" : "false" }
    doc: |
      If true - run samtools, if false - return sort_input, previously staged into output directory
  sort_compression_level:
    type: int?
    inputBinding:
      position: 7
      prefix: -l
    doc: |
      SORT: desired compression level for the final output file, ranging from 0 (uncompressed)
      or 1 (fastest but minimal compression) to 9 (best compression but slowest to write),
      similarly to gzip(1)'s compression level setting.
      If -l is not used, the default compression level will apply.
  sort_by_name:
    type: boolean?
    inputBinding:
      position: 8
      prefix: -n
    doc: |
      Sort by read names (i.e., the QNAME field) rather than by chromosomal coordinates
  sort_output_filename:
    type: string?
    inputBinding:
      position: 9
      prefix: -o
      valueFrom: |
        ${
            if (self == "" || inputs.trigger == false){
              return default_output_name();
            } else {
              return self;
            }
        }
    default: ''
    doc: |
      Write the final sorted output to FILE, rather than to standard output
  out_format:
    type: string?
    inputBinding:
      position: 10
      prefix: -O
    doc: |
      Specify output format (SAM, BAM, CRAM)
  threads:
    type: int?
    inputBinding:
      position: 11
      prefix: -@
    doc: |
      Set number of sorting and compression threads
  sort_input:
    type: File
    inputBinding:
      position: 12
    doc: |
      Input only in.sam|in.bam
outputs:
  sorted_file:
    type: File
    outputBinding:
      glob: |
        ${
            if (inputs.sort_output_filename == "" || inputs.trigger == false){
              return default_output_name();
            } else {
              return inputs.sort_output_filename;
            }
        }
baseCommand:
- bash
- -c
doc: |
  Tool to sort BAM/SAM file (set as input `sort_input`).
  If input `trigger` is set to `true` or isn't set at all (`true` is used by default), run `samtools sort`, return
  newly generated sorted BAM/SAM/CRAM file (the actual format of the output file depends on `out_format` value and
  extension of a filename set in `sort_output_filename`).
  If input `trigger` is set to `false`, return unchanged BAM/SAM file, previously staged into output directory.

  Before `baseCommand` is executed, input BAM/SAM file is staged into output directory (docker parameter `--workdir`),
  using `InitialWorkDirRequirement`. Setting `writable: true` makes cwl-runner to copy input BAM/SAM file and
  mount it to docker container with `rw` mode as part of `--workdir` (if set to false, the file staged into output
  directory will be mounted to docker container separately with `ro` mode). Because `samtools sort` can overwrite input
  BAM/SAM file and save output with the same name, we don't need to rename it (as we did for samtools rmdup).

  Trigger logic is implemented in bash script set by default in input `bash_script`. If first argment $0 (which is `trigger` input)
  is true, run `samtools sort` with the rest of the arguments. If $0 is not true, skip `samtools sort` and return
  input BAM/SAM file, previously staged into output directory.

  Input `trigger` is Boolean, but returns String, because of `valueFrom` field. The `valueFrom` is used, because if `trigger`
  is false, cwl-runner doesn't append this argument at all to the the `baseCommand` - new feature of CWL v1.0.2. Alternatively,
  `prefix` field could be used, but it causes changing logic in bash script saved in `bash_script` input.

  `default_output_name` function is used for generating output filename if input `sort_output_filename` is not set or in
  case when `trigger` is false and we need to return original BAM/SAM file staged into output directory.

  `ext` function returns output filename extension on the base of `out_format` value. If input `out_format` is not set,
  use `.bam` by default. If input `trigger` is false, `out_format` is ignored.

  If `trigger` is set to true (or true is used by default), but both `sort_output_filename` and `out_format` are not set,
  the output file format will be BAM by default.

  If `trigger` is set to true (or true is used by default) and `sort_output_filename` is set to some value, but `out_format`
  is not set, the actual format of output file will be defined on the filename extension set in `sort_output_filename`
  (if filename doesn't have any extension, BAM is used by default). When `out_format` is also set, it overwrites extension
  generated on the base of `sort_output_filename` value.
label: samtools-sort

cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var ext = function() {
      if (inputs.csi && !inputs.bai){
        return '.csi';
      } else {
        return '.bai';
      }
    };
  - var default_bam = function() {
      if (inputs.trigger == true){
        return inputs.sort_input.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+".bam";
      } else {
        return inputs.sort_input.location.split('/').slice(-1)[0];
      }
    };


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
        echo "Run: samtools sort " ${@:1}
        samtools sort "${@:1}"
      else
        echo "Skip samtools sort " ${@:1}
        SOURCE=${@:(-1):1}
        TARGET=$(basename "$SOURCE")
        cp $SOURCE $TARGET
      fi
    inputBinding:
      position: 5
    doc: |
      Bash function to run samtools sort with all input parameters or skip it if trigger is false

  bash_script_index:
    type: string?
    default: |
      #!/bin/bash
      if [ "$0" = "true" ]
      then
        echo "Run: samtools index " ${@:1}
        samtools index "${@:1}"
      else
        echo "Skip samtools index " ${@:1}
      fi
    inputBinding:
      position: 20
    doc: |
      Bash function to run samtools index with all input parameters or skip it if trigger is false

  trigger:
    type: boolean?
    default: true
    doc: |
      If true - run samtools, if false - return sort_input and optional index file in secondaryFiles.

  sort_compression_level:
    type: int?
    inputBinding:
      position: 11
      prefix: -l
    doc: |
      SORT: desired compression level for the final output file, ranging from 0 (uncompressed)
      or 1 (fastest but minimal compression) to 9 (best compression but slowest to write),
      similarly to gzip(1)'s compression level setting.
      If -l is not used, the default compression level will apply.

  sort_output_filename:
    type: string?
    inputBinding:
      position: 12
      prefix: -o
      valueFrom: |
        ${
            if (self == "" || inputs.trigger == false){
              return default_bam();
            } else {
              return self;
            }
        }
    default: ""
    doc: |
      Write the final sorted output to FILE. Only out.bam|out.cram.
      If output file extension is set to SAM, tool will fail on the index step

  threads:
    type: int?
    doc: |
      Set number of sorting and compression threads [1] (Only for sorting)

  sort_input:
    type: File
    inputBinding:
      position: 16
    doc: |
      Input only in.sam|in.bam|in.cram. Optionally could be supplemented with index file in secondaryFiles

  csi_interval:
    type: int?
    inputBinding:
      position: 24
      prefix: -m
    doc: |
      Set minimum interval size for CSI indices to 2^INT [14]

  csi:
    type: boolean?
    doc: |
      Generate CSI-format index for BAM files. If input isn't cram.

  bai:
    type: boolean?
    doc: |
      Generate BAI-format index for BAM files [default]. If input isn't cram.


outputs:

  bam_bai_pair:
    type: File
    outputBinding:
      glob: |
        ${
            if (inputs.sort_output_filename == "" || inputs.trigger == false){
              return default_bam();
            } else {
              return inputs.sort_output_filename;
            }
        }
    secondaryFiles:
      ${
          if (inputs.trigger == true){
            return self.basename + ext();
          } else {
            return inputs.sort_input.secondaryFiles?inputs.sort_input.secondaryFiles:"null";
          }
      }


baseCommand: [bash, '-c']


arguments:
#   run script sort position 5
  - valueFrom: |
      ${ return inputs.trigger ? "true" : "false" }
    position: 6
  # -l - position 11
  # -o sort_output_filename - position 12
  - valueFrom: bam
    position: 13
    prefix: -O
    # -n - position 14
  - valueFrom: $(inputs.threads?inputs.threads:1)
    position: 15
    prefix: -@
  # sort_input - position 16
  - valueFrom: ";"
    position: 17
    shellQuote: false

  - valueFrom: "bash"
    position: 18
  - valueFrom: "-c"
    position: 19
#   run script index position 20
  - valueFrom: |
      ${ return inputs.trigger ? "true" : "false" }
    position: 21
  - valueFrom: $(inputs.bai?'-b':inputs.csi?'-c':[])
    position: 23
    # -m - position 24
  - valueFrom: $(inputs.threads?inputs.threads:1)
    position: 25
    prefix: -@
  - valueFrom: |
      ${
          if (inputs.sort_output_filename == "" || inputs.trigger == false){
            return default_bam();
          } else {
            return inputs.sort_output_filename;
          }
      }
    position: 26
  - valueFrom: |
      ${
          if (inputs.sort_output_filename == "" || inputs.trigger == false){
            return default_bam() + ext();
          } else {
            return inputs.sort_output_filename + ext();
          }
      }
    position: 27

$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:mainEntity:
  $import: ./metadata/samtools-metadata.yaml

s:name: "samtools-sort-index"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/samtools-sort-index.cwl
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
  Tool to sort and index input BAM/SAM/CRAM.
  If input `trigger` is set to `true` or isn't set at all (`true` is used by default), run `samtools sort` and
  `samtools index`, return sorted BAM and BAI/CSI index file.
  If input `trigger` is set to `false`, return unchanged `sort_input` (BAM/SAM/CRAM) and index (BAI/CSI, if provided in
  `secondaryFiles`) files.

  Trigger logic is implemented in two bash scripts set by default as `bash_script_sort` and `bash_script_index` inputs.
  For both of then, if the first argument $0 (which is `trigger` input) is true, run `samtools sort/index` with the rest
  of the arguments. If $0 is not true, skip `samtools sort/index` and return `sort_input` and `secondaryFiles`
  (if provided).

  Input `trigger` is Boolean, but returns String, because of `valueFrom` field. The `valueFrom` is used, because if `trigger`
  is false, cwl-runner doesn't append this argument at all to the the `baseCommand` - new feature of CWL v1.0.2. Alternatively,
  `prefix` field could be used, but it causes changing in script logic.

  If using `sort_output_filename`, the output file extension should be `*.bam`, because `samtools sort` defines the output
  file format on the base of the file extension. If `*.sam` is sed as output filename, it cannot be usefully indexed
  by `samtools index`.

  `default_bam` function is used to generate output filename for `samtools sort` if input `sort_output_filename` is not
  set or when `trigger` is false and we need to return `sort_input` and `secondaryFiles` (if provided) files. Output
  filename is generated on the base of `sort_input` basename with `.bam` extension by default.

  `ext` function is used to return the index file extension (BAI/CSI) based on `csi` and `bai` inputs according to the
  following logic
    `csi` &&  `bai`  => BAI
   !`csi` && !`bai ` => BAI
    `csi` && !`bai ` => CSI

s:about: |
  Usage: samtools sort [options...] [in.bam]
  Options:
    -l INT     Set compression level, from 0 (uncompressed) to 9 (best)
    -m INT     Set maximum memory per thread; suffix K/M/G recognized [768M]
    -n         Sort by read name
    -o FILE    Write final output to FILE rather than standard output
    -T PREFIX  Write temporary files to PREFIX.nnnn.bam
        --input-fmt-option OPT[=VAL]
                 Specify a single input file format option in the form
                 of OPTION or OPTION=VALUE
    -O, --output-fmt FORMAT[,OPT[=VAL]]...
                 Specify output format (SAM, BAM, CRAM)
        --output-fmt-option OPT[=VAL]
                 Specify a single output file format option in the form
                 of OPTION or OPTION=VALUE
        --reference FILE
                 Reference sequence FASTA FILE [null]
    -@, --threads INT
                 Number of additional threads to use [0]

  Usage: samtools index [-bc] [-m INT] <in.bam> [out.index]
  Options:
    -b       Generate BAI-format index for BAM files [default]
    -c       Generate CSI-format index for BAM files
    -m INT   Set minimum interval size for CSI indices to 2^INT [14]
    -@ INT   Sets the number of threads [none]

cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function(ext) {
        let root = inputs.input_file.basename.split('.').slice(0,-1).join('.');
        return (root == "")?inputs.input_file.basename+ext:root+ext;
    };
- class: InitialWorkDirRequirement
  listing: |
    ${
      return  [
                {
                  "entry": inputs.annotation_filename,
                  "entryname": inputs.annotation_filename.basename,
                  "writable": true
                }
              ]
    }


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/atdp:v0.0.1

inputs:

  script:
    type: string?
    default: |
        #!/bin/bash
        set -- "$0" "$@"
        echo "Original arguments:"
        for i in "$@";
            do echo $i;
        done;
        set -- "$1" --a=$(basename "${2:4}") "${@:3}"
        echo Updated arguments:
        for i in "$@";
            do echo $i;
        done;
        refgene-sort -i "${2:4}" -o "${2:4}" -s "ORDER BY chrom, strand, CASE strand WHEN '+' THEN txStart WHEN '-' THEN txEnd END"
        atdp "$@"
    inputBinding:
      position: 1
    doc: |
      Bash function to run samtools sort with all input parameters or skip it if trigger is false

  input_file:
    type:
      - File
    inputBinding:
      position: 2
      prefix: --in=
      separate: false
    doc: |
      Input indexed BAM file (optionally +BAI index file in secondaryFiles)

  annotation_filename:
    type:
      - File
    inputBinding:
      position: 3
      prefix: --a=
      separate: false
    doc: |
      Annotation file, tsv

  output_filename:
    type:
      - "null"
      - string
    inputBinding:
      position: 4
      prefix: --out=
      separate: false
      valueFrom: |
        ${
            if (self == null){
              return default_output_filename('_atdp.tsv');
            } else {
              return self;
            }
        }
    default: null
    doc: |
      Base output file name, tsv

  log_filename:
    type:
      - "null"
      - string
    inputBinding:
      position: 5
      prefix: --log=
      separate: false
      valueFrom: |
        ${
            if (self == null){
              return default_output_filename('_atdp.log');
            } else {
              return self;
            }
        }
    default: null
    doc: |
      Log filename

  fragmentsize_bp:
    type:
      - "null"
      - int
    inputBinding:
      position: 6
      prefix: --f=
      separate: false
    doc: |
      Fragmentsize, int [150]

  avd_window_bp:
    type:
      - "null"
      - int
    inputBinding:
      position: 7
      prefix: --avd_window=
      separate: false
    doc: |
      Average tag density window, int [5000]

  avd_smooth_bp:
    type:
      - "null"
      - int
    inputBinding:
      position: 8
      prefix: --avd_smooth=
      separate: false
    doc: |
      Average smooth window (odd), int [0]

  ignore_chr:
    type:
      - "null"
      - string
    inputBinding:
      position: 9
      prefix: --sam_ignorechr=
      separate: false
    doc: |
      The chromosomes to be ignored, string

  double_chr:
    type:
      - "null"
      - string
    inputBinding:
      position: 10
      prefix: --sam_twicechr=
      separate: false
    doc: |
      The chromosomes to be doubled, string

  avd_heat_window_bp:
    type:
      - "null"
      - int
    inputBinding:
      position: 11
      prefix: --avd_heat_window=
      separate: false
    doc: |
      Average tag density window for heatmap

  mapped_reads:
    type:
      - "null"
      - int
    inputBinding:
      position: 12
      prefix: --m=
      separate: false
    doc: |
      Mapped reads number

  index_file:
    type:
      - "null"
      - File
    inputBinding:
      position: 13
      prefix: --index=
      separate: false
    doc: |
      Index file

outputs:
  log_file:
    type: File
    outputBinding:
      glob: |
        ${
          if (inputs.log_filename == null){
            return default_output_filename('_atdp.log');
          } else {
            return inputs.log_filename;
          }
        }

  result_file:
    type: File
    outputBinding:
      glob: |
        ${
          if (inputs.output_filename == null){
            return default_output_filename('_atdp.tsv');
          } else {
            return inputs.output_filename;
          }
        }

baseCommand: [bash, '-c']

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: ./metadata/atdp-metadata.yaml

s:name: "atdp"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/atdp.cwl
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
  Tool calculates average tag density profile around all annotated TSS.

  `default_output_filename` function returns output filename with sufix set as `ext` argument. Function is called when
  either `output_filename` or `log_filename` inputs are not provided.

  Before running `baseCommand`, annotaion file `annotation_filename` is staged into output directory
  (Docker's `--workdir`) with `"writable": true` (to allow to overwrite it by `refgene-sort`).

  To run `atdp` index file should be provided (either in `secondaryFiles` of `input_file` or as separate input
  `index_file`)

  `baseCommand` runs bash script from `script` input. Script runs `refgene-sort` to sort annotaion file and then runs
  `atdp`. `refgene-sort` - utility to sort refgene annotation files, using MySQL syntax.

  Optionally (with cwltool==1.0.20171107133715), script can be simplified to
      #!/bin/bash
      set -- "$0" "$@"
      refgene-sort -i "${2:4}" -o "${2:4}" -s "ORDER BY chrom, strand, CASE strand WHEN '+' THEN txStart WHEN '-' THEN txEnd END"
      atdp "$@"
  because `set -- "$1" --a=$(basename "${2:4}") "${@:3}"` is not needed anymore.

s:about: |
  Usage:
  atdp [options] --in=pathToFile --a=pathtoFile --out=pathToFile
    --a                        	Tab-separated annotation file
    --avd_bsmooth              	Average smooth window for gene body
    --avd_guid                 	Genelist uid
    --avd_heat_window          	Average tag density window for heatmap
    --avd_smooth               	Average smooth window (odd)
    --avd_window               	Average tag density window
    --f                        	Fragmentsize, bp
    --in                       	Input BAM file
    --index                    	Input index bai file
    --log                      	Log file name (default is ./logfile_def.log)
    --out                      	Base output filename
    --sam_ignorechr            	Which chromosome to ignore
    --sam_twicechr             	Which chromosome to double
    --m                         Mapped reads number

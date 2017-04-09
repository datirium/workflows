#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: ./metadata/envvar-global.yml
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_name = function(input_staged, sufix) {
      input_staged = input_staged || false;
      sufix = sufix || "";
      if (inputs.trigger == false && input_staged){
        return input_staged.basename;
      } else {
        if (Object.prototype.toString.call(inputs.treatment) === '[object Array]'){
          return inputs.treatment[0].location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+sufix;
        } else {
          return inputs.treatment.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+sufix;
        }
      }
    }

- class: InitialWorkDirRequirement
  listing: |
    ${
      var listing = []
      if (inputs.peak_xls_file_staged){
        listing.push(inputs.peak_xls_file_staged);
      }
      if (inputs.narrow_peak_file_staged){
        listing.push(inputs.narrow_peak_file_staged);
      }
      if (inputs.broad_peak_file_staged){
        listing.push(inputs.broad_peak_file_staged);
      }
      if (inputs.gapped_peak_file_staged){
        listing.push(inputs.gapped_peak_file_staged);
      }
      if (inputs.peak_summits_file_staged){
        listing.push(inputs.peak_summits_file_staged);
      }
      if (inputs.moder_r_file_staged){
        listing.push(inputs.moder_r_file_staged);
      }
      if (inputs.treat_pileup_bdg_file_staged){
        listing.push(inputs.treat_pileup_bdg_file_staged);
      }
      if (inputs.control_lambda_bdg_file_staged){
        listing.push(inputs.control_lambda_bdg_file_staged);
      }
      return listing;
    }

hints:
- class: DockerRequirement
  dockerPull: scidap/macs2:v2.1.1
  dockerFile: >
    $import: ./dockerfiles/macs2-Dockerfile

inputs:

  script:
    type: string?
    default: |
      #!/bin/bash
      if [ "$0" = True ]
      then
        echo "Run: macs2 callpeak " ${@:1}
        ls | grep -v ${@: -1}.log | xargs rm
        macs2 callpeak "${@:1}"
      else
        echo "Skip macs2 callpeak " ${@:1}
      fi
    inputBinding:
      position: 1
    doc: |
      Bash function to run MACS2 callpeak with all input parameters or skip it if trigger is false

  trigger:
    type: boolean?
    default: true
    inputBinding:
      position: 2
    doc: |
      If true - run MACS2, if false - return staged files

  peak_xls_file_staged:
    type: File?
    doc: For staging in a case of trigger is set to false
  narrow_peak_file_staged:
    type: File?
    doc: For staging in a case of trigger is set to false
  broad_peak_file_staged:
    type: File?
    doc: For staging in a case of trigger is set to false
  gapped_peak_file_staged:
    type: File?
    doc: For staging in a case of trigger is set to false
  peak_summits_file_staged:
    type: File?
    doc: For staging in a case of trigger is set to false
  moder_r_file_staged:
    type: File?
    doc: For staging in a case of trigger is set to false
  treat_pileup_bdg_file_staged:
    type: File?
    doc: For staging in a case of trigger is set to false
  control_lambda_bdg_file_staged:
    type: File?
    doc: For staging in a case of trigger is set to false

  treatment:
    type:
      - File
      - type: array
        items: File
    inputBinding:
      position: 10
      prefix: -t
    doc: |
      This is the only REQUIRED parameter for MACS. File can be in any supported format specified by –format option.
      Check –format for detail. If you have more than one alignment files, you can specify them as `-t A B C`.
      MACS will pool up all these files together.

  name:
    type:
      - "null"
      - string
    inputBinding:
      position: 999
      prefix: -n
      valueFrom: |
        ${
            if (self == null || inputs.trigger == false){
              return default_name();
            } else {
              return self;
            }
        }
    default: null
    doc: |
      The name string of the experiment. MACS will use this string NAME to create output files like ‘NAME_peaks.xls’,
      ‘NAME_negative_peaks.xls’, ‘NAME_peaks.bed’ , ‘NAME_summits.bed’, ‘NAME_model.r’ and so on.
      So please avoid any confliction between these filenames and your existing files.
      DEFAULT: generated on the base of the treatment input

  control:
    type:
      - "null"
      - File
      - type: array
        items: File
    inputBinding:
      position: 12
      prefix: -c
    doc: |
      The control or mock data file. Please follow the same direction as for -t/–treatment.

  format:
    type:
      - "null"
      - string
    inputBinding:
      position: 13
      prefix: -f
    doc: |
      {AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,ELANDEXPORT,BOWTIE,BAMPE}, --format
      {AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,ELANDEXPORT,BOWTIE,BAMPE} Format of tag file,
      "AUTO", "BED" or "ELAND" or "ELANDMULTI" or "ELANDEXPORT" or "SAM" or "BAM"
      or "BOWTIE" or "BAMPE". The default AUTO option will let MACS decide which format
      the file is. Note that MACS can''t detect "BAMPE" or "BEDPE" format with "AUTO",
      and you have to implicitly specify the format for "BAMPE" and "BEDPE".
      DEFAULT: AUTO

  genome_size:
    type:
      - "null"
      - string
    inputBinding:
      position: 14
      prefix: -g
    doc: |
      It’s the mappable genome size or effective genome size which is defined as the genome size which can be sequenced.
      Because of the repetitive features on the chromsomes, the actual mappable genome size will be smaller than the
      original size, about 90% or 70% of the genome size. The default hs – 2.7e9 is recommended for UCSC human hg18
      assembly. Here are all precompiled parameters for effective genome size:
        hs:	2.7e9
        mm:	1.87e9
        ce:	9e7
        dm:	1.2e8
      DEFAULT: hs

  keep_dup:
    type:
      - "null"
      - string
    inputBinding:
      position: 15
      prefix: --keep-dup
    doc: |
      It controls the MACS behavior towards duplicate tags
      at the exact same location -- the same coordination
      and the same strand. The 'auto' option makes MACS
      calculate the maximum tags at the exact same location
      based on binomal distribution using 1e-5 as pvalue
      cutoff; and the 'all' option keeps every tags. If an
      integer is given, at most this number of tags will be
      kept at the same location. The default is to keep one
      tag at the same location.
      DEFAULT: 1

  buffer_size:
    type:
      - "null"
      - int
    inputBinding:
      position: 16
      prefix: --buffer-size
    doc: |
      Buffer size for incrementally increasing internal
      array size to store reads alignment information. In
      most cases, you don't have to change this parameter.
      However, if there are large number of
      chromosomes/contigs/scaffolds in your alignment, it's
      recommended to specify a smaller buffer size in order
      to decrease memory usage (but it will take longer time
      to read alignment files). Minimum memory requested for
      reading an alignment file is about # of CHROMOSOME *
      BUFFER_SIZE * 2 Bytes.
      DEFAULT: 100000

  bdg:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 17
      prefix: --bdg
    doc: |
      If this flag is on, MACS will store the fragment pileup, control lambda, -log10pvalue and -log10qvalue scores in bedGraph files.
      The bedGraph files will be stored in current directory named NAME+’_treat_pileup.bdg’ for treatment data, NAME+’_control_lambda.bdg’
      for local lambda values from control, NAME+’_treat_pvalue.bdg’ for Poisson pvalue scores (in -log10(pvalue) form),
      and NAME+’_treat_qvalue.bdg’ for q-value scores from
      Benjamini–Hochberg–Yekutieli procedure <http://en.wikipedia.org/wiki/False_discovery_rate#Dependent_tests>
      DEFAULT: False

  trackline:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 18
      prefix: --trackline
    doc: |
      Tells MACS to include trackline with bedGraph files. To include this trackline
      while displaying bedGraph at UCSC genome browser, can show name and description
      of the file as well. Require -B to be set.
      DEFAULT: False

  spmr:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 19
      prefix: --SPMR
    doc: |
      If True, MACS will save signal per million reads for fragment pileup profiles.
      Require --bdg to be set.
      DEFAULT: False

  tsize:
    type:
      - "null"
      - int
    inputBinding:
      position: 20
      prefix: --tsize
    doc: |
      Tag size. This will overide the auto detected tag size.
      DEFAULT: False

  bw:
    type:
      - "null"
      - int
    inputBinding:
      position: 21
      prefix: --bw
    doc: |
      Band width for picking regions to compute  fragment  size.  This
      value is only used while building the shifting model
      DEFAULT: 300

  mfold:
    type:
      - "null"
      - string
    inputBinding:
      position: 22
      prefix: -m
      valueFrom: |
        ${
          return self.replace(/\s+/g, ' ').split(' ');
        }
    doc: |
      Select the regions within MFOLD range of high-
      confidence enrichment ratio against background to
      build model. Fold-enrichment in regions must be lower
      than upper limit, and higher than the lower limit. Use
      as "-m 10 30"
      DEFAULT: 5 50

  fix_bimodal:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 23
      prefix: --fix-bimodal
    doc: |
      Whether turn on the auto pair model process. If set, when MACS failed to
      build paired model, it will use the nomodel settings, the --exsize parameter
      to extend each tags towards 3'' direction. Not to use this automate fixation
      is a default behavior now.
      DEFAULT: False

  nomodel:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 24
      prefix: --nomodel
    doc: |
      Whether or not to build the shifting model. If True,  MACS will not build
      model. by default it means  shifting size = 100, try to set extsize to change it.
      DEFAULT: False

  shift:
    type:
      - "null"
      - int
    inputBinding:
      position: 25
      prefix: --shift
    doc: |
      (NOT the legacy --shiftsize option!) The arbitrary shift in bp. Use discretion
      while setting it other than default value. When NOMODEL is set, MACS will use
      this value to move cutting ends (5'') towards 5''->3'' direction then apply
      EXTSIZE to extend them to fragments. When this value is negative, ends will
      be moved toward 3''->5'' direction. Recommended to keep it as default 0 for
      ChIP-Seq datasets, or -1 * half of EXTSIZE together with EXTSIZE option for
      detecting enriched cutting loci such as certain DNAseI-Seq datasets. Note, you
      can''t set values other than 0 if format is BAMPE for paired-end data.
      DEFAULT: 0

  extsize:
    type:
      - "null"
      - int
    inputBinding:
      position: 26
      prefix: --extsize
    doc: |
      The arbitrary extension size in bp. When nomodel is  true, MACS will use
      this value as fragment size to  extend each read towards 3'' end, then pile
      them up.  It''s exactly twice the number of obsolete SHIFTSIZE.  In previous
      language, each read is moved 5''->3''  direction to middle of fragment by 1/2
      d, then  extended to both direction with 1/2 d. This is  equivalent to say each
      read is extended towards 5''->3''  into a d size fragment.EXTSIZE
      and  SHIFT can be combined when necessary. Check SHIFT  option.
      DEFAULT: 200

  q_value:
    type:
      - "null"
      - float
    inputBinding:
      position: 27
      prefix: -q
    doc: |
      Minimum FDR (q-value) cutoff for peak detection. -q, and
      -p are mutually exclusive.
      DEFAULT: 0.05

  p_value:
    type:
      - "null"
      - float
    inputBinding:
      position: 28
      prefix: -p
    doc: |
      Pvalue cutoff for peak detection. DEFAULT: not set.  -q, and -p are mutually
      exclusive. If pvalue cutoff is set, qvalue will not be calculated and reported
      as -1  in the final .xls file.
      DEFAULT: null

  to_large:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 29
      prefix: --to-large
    doc: |
      When set, scale the small sample up to the bigger sample. By default, the
      bigger dataset will be scaled down towards the smaller dataset, which will lead
      to smaller p/qvalues and more specific results. Keep in mind that scaling down
      will bring down background noise more.
      DEFAULT: False

  ratio:
    type:
      - "null"
      - float
    inputBinding:
      position: 30
      prefix: --ratio
    doc: |
      When set, use a custom scaling ratio of ChIP/control
      (e.g. calculated using NCIS) for linear scaling.
      DEFAULT: null

  down_sample:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 31
      prefix: --down-sample
    doc: |
      When set, random sampling method will scale down the bigger sample. By default,
      MACS uses linear scaling. Warning: This option will make your result unstable
      and irreproducible since each time, random reads would be selected. Consider
      to use ''randsample'' script instead. <not implmented>If used together with
      --SPMR, 1 million unique reads will be randomly picked.</not implemented> Caution:
      due to the implementation, the final number of selected reads may not be as
      you expected!
      DEFAULT: False

  seed:
    type:
      - "null"
      - int
    inputBinding:
      position: 32
      prefix: --seed
    doc: |
      Set the random seed while down sampling data. Must be
      a non-negative integer in order to be effective.
      DEFAULT: null

  nolambda:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 33
      prefix: --nolambda
    doc: |
      If True, MACS will use fixed background lambda as local lambda for every
      peak region. Normally, MACS calculates a dynamic local lambda to reflect the
      local bias due to potential chromatin structure.
      DEFAULT: False

  slocal:
    type:
      - "null"
      - int
    inputBinding:
      position: 34
      prefix: --slocal
    doc: |
      The small nearby region in basepairs to calculate
      dynamic lambda. This is used to capture the bias near
      the peak summit region. Invalid if there is no control
      data. If you set this to 0, MACS will skip slocal
      lambda calculation. *Note* that MACS will always
      perform a d-size local lambda calculation. The final
      local bias should be the maximum of the lambda value
      from d, slocal, and llocal size windows.
      DEFAULT: 1000

  llocal:
    type:
      - "null"
      - int
    inputBinding:
      position: 35
      prefix: --llocal
    doc: |
      The large nearby region in basepairs to calculate
      dynamic lambda. This is used to capture the surround
      bias. If you set this to 0, MACS will skip llocal
      lambda calculation. *Note* that MACS will always
      perform a d-size local lambda calculation. The final
      local bias should be the maximum of the lambda value
      from d, slocal, and llocal size windows.
      DEFAULT: 10000.

  broad:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 36
      prefix: --broad
    doc: |
      If set, MACS will try to call broad peaks by linking nearby highly enriched
      regions. The linking region is controlled by another cutoff through --linking-cutoff.
      The maximum linking region length is 4 times of d from MACS.
      DEFAULT: False

  broad_cutoff:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 37
      prefix: --broad-cutoff
    doc: |
      Cutoff for broad region. This option is not available
      unless --broad is set. If -p is set, this is a pvalue
      cutoff, otherwise, it's a qvalue cutoff.
      DEFAULT: 0.1

  cutoff_analysis:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 38
      prefix: --cutoff-analysis
    doc: |
      While set, MACS2 will analyze number or total length of peaks that can be
      called by different p-value cutoff then output a summary table to help user
      decide a better cutoff. The table will be saved in NAME_cutoff_analysis.txt
      file. Note, minlen and maxgap may affect the results. WARNING: May take ~30
      folds longer time to finish.
      DEFAULT: False

  call_summits:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 39
      prefix: --call-summits
    doc: |
      If set, MACS will use a more sophisticated signal processing approach to
      find subpeak summits in each enriched peak region.
      DEFAULT: False

  fe-cutoff:
    type:
      - "null"
      - float
    inputBinding:
      position: 40
      prefix: --fe-cutoff
    doc: |
      When set, the value will be used to filter out peaks
      with low fold-enrichment. Note, MACS2 use 1.0 as
      pseudocount while calculating fold-enrichment.
      DEFAULT: 1.0

  verbose:
    type:
      - "null"
      - int
    inputBinding:
      position: 41
      prefix: --verbose
    doc: |
      Log level

outputs:
  peak_xls_file:
    type: File?
    outputBinding:
      glob:
        ${
            if (inputs.name == null || inputs.trigger == false){
              return default_name(inputs.peak_xls_file_staged, '_peaks.xls');
            } else {
              return inputs.name + '_peaks.xls';
            }
        }

  narrow_peak_file:
    type: File?
    outputBinding:
      glob:
        ${
            if (inputs.name == null || inputs.trigger == false){
              return default_name(inputs.narrow_peak_file_staged, '_peaks.narrowPeak');
            } else {
              return inputs.name + '_peaks.narrowPeak';
            }
        }

  broad_peak_file:
    type: File?
    outputBinding:
      glob:
        ${
            if (inputs.name == null || inputs.trigger == false){
              return default_name(inputs.broad_peak_file_staged, '_peaks.broadPeak');
            } else {
              return inputs.name + '_peaks.broadPeak';
            }
        }

  gapped_peak_file:
    type: File?
    outputBinding:
      glob:
        ${
            if (inputs.name == null || inputs.trigger == false){
              return default_name(inputs.gapped_peak_file_staged, '_peaks.gappedPeak');
            } else {
              return inputs.name + '_peaks.gappedPeak';
            }
        }

  peak_summits_file:
    type: File?
    outputBinding:
      glob:
        ${
            if (inputs.name == null || inputs.trigger == false){
              return default_name(inputs.peak_summits_file_staged, '_summits.bed');
            } else {
              return inputs.name + '_summits.bed';
            }
        }

  moder_r_file:
    type: File?
    outputBinding:
      glob:
        ${
            if (inputs.name == null || inputs.trigger == false){
              return default_name(inputs.moder_r_file_staged, '_model.r');
            } else {
              return inputs.name + '_model.r';
            }
        }

  treat_pileup_bdg_file:
    type: File?
    outputBinding:
      glob:
        ${
            if (inputs.name == null || inputs.trigger == false){
              return default_name(inputs.treat_pileup_bdg_file_staged, '_treat_pileup.bdg');
            } else {
              return inputs.name + '_treat_pileup.bdg';
            }
        }

  control_lambda_bdg_file:
    type: File?
    outputBinding:
      glob:
        ${
            if (inputs.name == null || inputs.trigger == false){
              return default_name(inputs.control_lambda_bdg_file_staged, '_control_lambda.bdg');
            } else {
              return inputs.name + '_control_lambda.bdg';
            }
        }

  macs_log:
    type: File?
    outputBinding:
      glob: |
        ${
            if (inputs.name == null || inputs.trigger == false){
              return default_name(null, '.log');
            } else {
              return inputs.name + '.log';
            }
        }

baseCommand: [bash, '-c']

arguments:
  - valueFrom:
      ${
          if (inputs.name == null || inputs.trigger == false ){
            return ' 2> ' + default_name(null, '.log');
          } else {
            return ' 2> ' + inputs.name + '.log';
          }
      }
    position: 100000
    shellQuote: false

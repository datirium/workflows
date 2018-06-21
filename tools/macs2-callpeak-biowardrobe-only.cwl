cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function(sufix) {
      sufix = sufix || "_macs";
      if (Object.prototype.toString.call(inputs.treatment_file) === '[object Array]'){
        let basename = inputs.treatment_file[0].basename;
        let root = basename.split('.').slice(0,-1).join('.');
        return (root == "")?basename+sufix:root+sufix;
      } else {
        let basename = inputs.treatment_file.basename;
        let root = basename.split('.').slice(0,-1).join('.');
        return (root == "")?basename+sufix:root+sufix;
      }
    }
- class: InitialWorkDirRequirement
  listing:
    - entryname: run.py
      entry: |
          #!/usr/bin/env python
          import subprocess
          import sys
          import argparse
          import re

          def arg_parser():
              """Returns argument parser. used only to get -n and --extsize values"""
              parser = argparse.ArgumentParser(description='MACS2 Biowardrobe parser')
              parser.add_argument("-n", help="MACS2 string name of the experiment")
              parser.add_argument("--extsize", help="MACS2 arbitrary extension size in bp")
              return parser

          def get_info (name):
              fragments, islands = 0, 0
              with open(name+'_peaks.xls', 'r') as infile:
                  for line in infile:
                      if re.match('^# d = ', line):
                          fragments = int(line.split('d = ')[1])
                          continue
                      if re.match('^#', line):
                          continue
                      if line.strip() != "":
                          islands = islands + 1
              islands = islands - 1
              return fragments, islands

          args, _ = arg_parser().parse_known_args()
          argv=' '.join(sys.argv[1:])

          extsize=args.extsize
          name=args.n

          runtime_error=False
          calculated_fragment_size=0
          expected_fragment_size=0
          island_count=0

          try:
              subprocess.check_output(argv, shell=True)
          except subprocess.CalledProcessError as e:
              if '--nomodel' in argv:
                  print "MACS2 critical error: ", str(e)
                  sys.exit(1) # stop the whole workflow
              print "MACS2 will be rerun to fix not critical error: ", str(e)
              runtime_error=True
              expected_fragment_size=extsize
          else:
              calculated_fragment_size,island_count=get_info(name)
              expected_fragment_size=calculated_fragment_size

          if runtime_error or (calculated_fragment_size < 80 and ('--nomodel' not in argv)):
              try:
                  subprocess.check_output(argv+' --nomodel', shell=True)
              except subprocess.CalledProcessError as e:
                  print "MACS2 forced error: ", str(e)
                  sys.exit(1) # stop the whole workflow
              else:
                  calculated_fragment_size,island_count=get_info(name)

          with open(name.replace('_macs','')+'_fragment_stat.tsv', 'w') as output_file:
              output_file.write(' '.join([str(calculated_fragment_size),str(expected_fragment_size),str(island_count)]))


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/macs2:v2.1.1


inputs:

  treatment_file:
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

  output_prefix:
    type:
      - "null"
      - string
    inputBinding:
      position: 999
      prefix: -n
      valueFrom: |
        ${
            if (self == null){
              return default_output_filename();
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

  control_file:
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

  format_mode:
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
    type: int
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

  fe_cutoff:
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
            if (inputs.output_prefix == null){
              return default_output_filename('_macs_peaks.xls');
            } else {
              return inputs.output_prefix + '_peaks.xls';
            }
        }

  narrow_peak_file:
    type: File?
    outputBinding:
      glob:
        ${
            if (inputs.output_prefix == null){
              return default_output_filename('_macs_peaks.narrowPeak');
            } else {
              return inputs.output_prefix + '_peaks.narrowPeak';
            }
        }

  broad_peak_file:
    type: File?
    outputBinding:
      glob:
        ${
            if (inputs.output_prefix == null){
              return default_output_filename('_macs_peaks.broadPeak');
            } else {
              return inputs.output_prefix + '_peaks.broadPeak';
            }
        }

  gapped_peak_file:
    type: File?
    outputBinding:
      glob:
        ${
            if (inputs.output_prefix == null){
              return default_output_filename('_macs_peaks.gappedPeak');
            } else {
              return inputs.output_prefix + '_peaks.gappedPeak';
            }
        }

  peak_summits_file:
    type: File?
    outputBinding:
      glob:
        ${
            if (inputs.output_prefix == null){
              return default_output_filename('_macs_summits.bed');
            } else {
              return inputs.output_prefix + '_summits.bed';
            }
        }

  moder_r_file:
    type: File?
    outputBinding:
      glob:
        ${
            if (inputs.output_prefix == null){
              return default_output_filename('_macs_model.r');
            } else {
              return inputs.output_prefix + '_model.r';
            }
        }

  treat_pileup_bdg_file:
    type: File?
    outputBinding:
      glob:
        ${
            if (inputs.output_prefix == null){
              return default_output_filename('_macs_treat_pileup.bdg');
            } else {
              return inputs.output_prefix + '_treat_pileup.bdg';
            }
        }

  control_lambda_bdg_file:
    type: File?
    outputBinding:
      glob:
        ${
            if (inputs.output_prefix == null){
              return default_output_filename('_macs_control_lambda.bdg');
            } else {
              return inputs.output_prefix + '_control_lambda.bdg';
            }
        }

  macs_log:
    type: File?
    outputBinding:
      glob: |
        ${
            if (inputs.output_prefix == null){
              return default_output_filename('_macs.log');
            } else {
              return inputs.output_prefix + '.log';
            }
        }

  macs2_stat_file:
    type: File
    outputBinding:
      glob: |
        ${
            if (inputs.output_prefix == null){
              return default_output_filename('_fragment_stat.tsv');
            } else {
              return inputs.output_prefix + '_fragment_stat.tsv';
            }
        }

  macs2_fragments_calculated:
    type: int
    outputBinding:
      loadContents: true
      glob: |
        ${
            if (inputs.output_prefix == null){
              return default_output_filename('_fragment_stat.tsv');
            } else {
              return inputs.output_prefix + '_fragment_stat.tsv';
            }
        }
      outputEval: $(parseInt(self[0].contents.split(' ')[0]))

  macs2_fragments_expected:
    type: int
    outputBinding:
      loadContents: true
      glob: |
        ${
            if (inputs.output_prefix == null){
              return default_output_filename('_fragment_stat.tsv');
            } else {
              return inputs.output_prefix + '_fragment_stat.tsv';
            }
        }
      outputEval: $(parseInt(self[0].contents.split(' ')[1]))

  macs2_islands:
    type: int
    outputBinding:
      loadContents: true
      glob: |
        ${
            if (inputs.output_prefix == null){
              return default_output_filename('_fragment_stat.tsv');
            } else {
              return inputs.output_prefix + '_fragment_stat.tsv';
            }
        }
      outputEval: $(parseInt(self[0].contents.split(' ')[2]))




baseCommand: [python, 'run.py', 'macs2 callpeak']

arguments:
  - valueFrom:
      ${
          if (inputs.output_prefix == null){
            return ' 2> ' + default_output_filename('_macs.log');
          } else {
            return ' 2> ' + inputs.output_prefix + '.log';
          }
      }
    position: 100000
    shellQuote: true

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: ./metadata/macs2-metadata.yaml

s:name: "macs2-callpeak"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/macs2-callpeak.cwl
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
  Tool runs peak calling using MACS2 within the custom script `run.py`

  `default_output_filename` function returns default output prefix for all generated output files. Function is called
  when `output_prefix` input is not set. Optional argument `sufix` is equal to `_macs2` by default.
  `sufix` is hardcoded and should not be changed (if necessary, update `_macs2` in all tool outputs and
  `run.py` script too). `sufix` is included in the filenames of all output files except those, who ends
  with `_fragment_stat.tsv`. For the rest of the output files `_macs2` is already included as a part of the
  `sufix` argument when `default_output_filename` function is called. This is mostly done because of the
  naming conventions for output files in supplemental project `https://github.com/Barski-lab/biowardrobe-analysis`
  (if `_fragment_stat.tsv` is updated to `_macs2_fragment_stat.tsv` in `biowardrobe-analysis`, the way we use
  `default_output_filename` function can be simplified). The basename for output filename is taken from
  `treatment_file` input (if array - the first item is used).

  `baseCommand` runs Python script `run.py` staged into the output directory. Within this script MACS2 is called one
  ortwo times as a subprocess with all provided input parameters.

  `extsize` is mandatory parameter. It will be ignored by MASC2, when `--nomodel` is not set.

  The aim of `run.py` script is to return `calculated_fragment_size`, `expected_fragment_size` and `island_count` as
  a text file.

  Scenarios:
  1. Run MACS2 with all input parameters.
     If MACS2 first run failed and `--nomodel` wasn't set (it means that `--extsize` wasn't used), set
     `expected_fragment_size` to `--extsize` value and try to rerun MACS2 with `--nomodel` parameter. If MACS2 second
     run failed too - exit workflow, if it didn't fail - get values for `calculated_fragment_size` and `island_count`
     from generated output file.
     Returned data:
      `calculated_fragment_size` - the value of fragment size from MACS2 second run
      `expected_fragment_size` - the value of fragment size equal to `--extsize` (set when MACS2 first run failed)
      `island_count` - the number of islands from MACS2 second run

  2. Run MACS2 with all input parameters.
     If MACS2 first run failed and `--nomodel` was already set, exit workflow, because there is no need to rerun
     MACS2, it will fail too.

  3. Run MACS2 with all input parameters.
     If MACS2 first run didn't fail, we get `calculated_fragment_size` and `island_count` from generated output file,
     and set `expected_fragment_size` to be equal to `calculated_fragment_size`. If `calculated_fragment_size` is
     bigger or equal to 80, there is no need to rerun MACS2.
     Returned data:
      `calculated_fragment_size` - the value of fragment size from MACS2 first run
      `expected_fragment_size` - the value of fragment size equal to `calculated_fragment_size`
      `island_count` - the number of islands from MACS2 first run

  4. Run MACS2 with all input parameters.
     If MACS2 first run didn't fail, we get `calculated_fragment_size` and `island_count` from generated output file,
     and set `expected_fragment_size` to be equal to `calculated_fragment_size`. If `calculated_fragment_size` is less
     then 80 and MACS2 was run without '--nomodel' (it means we didn't force it to use the fixed fragment size
     `--extsize`), we rerun MACS2 with `--nomodel` argument. If MACS2 second run failed too - exit workflow, if not -
     get values for `calculated_fragment_size` and `island_count` from generated output file.
     Returned data:
      `calculated_fragment_size` - the value of fragment size from MACS2 second run
      `expected_fragment_size` - the value of fragment size from MACS2 first run
      `island_count` - the number of islands from MACS2 second run

  5. Run MACS2 with all input parameters.
     If MACS2 first run didn't fail and was already run with '--nomodel' parameter (forced to use fixed fragment size
     `--extsize`), we get `calculated_fragment_size` and `island_count` from generated output file, and set
     `expected_fragment_size` to be equal to `calculated_fragment_size`. No mater which value we got for
     `calculated_fragment_size`, we don't rerun MACS2.
     Returned data:
      `calculated_fragment_size` - the value of fragment size got from MACS2 first run
      `expected_fragment_size` - the value of fragment size equal to `calculated_fragment_size`
      `island_count` - the number of islands got from MACS2 first run of MACS2

s:about: |
  usage: macs2 callpeak [-h] -t TFILE [TFILE ...] [-c [CFILE [CFILE ...]]]
                        [-f {AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,ELANDEXPORT,BOWTIE,BAMPE,BEDPE}]
                        [-g GSIZE] [--keep-dup KEEPDUPLICATES]
                        [--buffer-size BUFFER_SIZE] [--outdir OUTDIR] [-n NAME]
                        [-B] [--verbose VERBOSE] [--trackline] [--SPMR]
                        [-s TSIZE] [--bw BW] [-m MFOLD MFOLD] [--fix-bimodal]
                        [--nomodel] [--shift SHIFT] [--extsize EXTSIZE]
                        [-q QVALUE | -p PVALUE] [--to-large] [--ratio RATIO]
                        [--down-sample] [--seed SEED] [--tempdir TEMPDIR]
                        [--nolambda] [--slocal SMALLLOCAL] [--llocal LARGELOCAL]
                        [--broad] [--broad-cutoff BROADCUTOFF]
                        [--cutoff-analysis] [--call-summits]
                        [--fe-cutoff FECUTOFF]

  optional arguments:
    -h, --help            show this help message and exit

  Input files arguments:
    -t TFILE [TFILE ...], --treatment TFILE [TFILE ...]
                          ChIP-seq treatment file. If multiple files are given
                          as '-t A B C', then they will all be read and pooled
                          together. REQUIRED.
    -c [CFILE [CFILE ...]], --control [CFILE [CFILE ...]]
                          Control file. If multiple files are given as '-c A B
                          C', they will be pooled to estimate ChIP-seq
                          background noise.
    -f {AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,ELANDEXPORT,BOWTIE,BAMPE,BEDPE}, --format {AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,ELANDEXPORT,BOWTIE,BAMPE,BEDPE}
                          Format of tag file, "AUTO", "BED" or "ELAND" or
                          "ELANDMULTI" or "ELANDEXPORT" or "SAM" or "BAM" or
                          "BOWTIE" or "BAMPE" or "BEDPE". The default AUTO
                          option will let MACS decide which format (except for
                          BAMPE and BEDPE which should be implicitly set) the
                          file is. Please check the definition in README. Please
                          note that if the format is set as BAMPE or BEDPE,
                          MACS2 will call its special Paired-end mode to call
                          peaks by piling up the actual ChIPed fragments defined
                          by both aligned ends, instead of predicting the
                          fragment size first and extending reads. Also please
                          note that the BEDPE only contains three columns, and
                          is NOT the same BEDPE format used by BEDTOOLS.
                          DEFAULT: "AUTO"
    -g GSIZE, --gsize GSIZE
                          Effective genome size. It can be 1.0e+9 or 1000000000,
                          or shortcuts:'hs' for human (2.7e9), 'mm' for mouse
                          (1.87e9), 'ce' for C. elegans (9e7) and 'dm' for
                          fruitfly (1.2e8), Default:hs
    --keep-dup KEEPDUPLICATES
                          It controls the MACS behavior towards duplicate tags
                          at the exact same location -- the same coordination
                          and the same strand. The 'auto' option makes MACS
                          calculate the maximum tags at the exact same location
                          based on binomal distribution using 1e-5 as pvalue
                          cutoff; and the 'all' option keeps every tags. If an
                          integer is given, at most this number of tags will be
                          kept at the same location. Note, if you've used
                          samtools or picard to flag reads as 'PCR/Optical
                          duplicate' in bit 1024, MACS2 will still read them
                          although the reads may be decided by MACS2 as
                          duplicate later. The default is to keep one tag at the
                          same location. Default: 1
    --buffer-size BUFFER_SIZE
                          Buffer size for incrementally increasing internal
                          array size to store reads alignment information. In
                          most cases, you don't have to change this parameter.
                          However, if there are large number of
                          chromosomes/contigs/scaffolds in your alignment, it's
                          recommended to specify a smaller buffer size in order
                          to decrease memory usage (but it will take longer time
                          to read alignment files). Minimum memory requested for
                          reading an alignment file is about # of CHROMOSOME *
                          BUFFER_SIZE * 2 Bytes. DEFAULT: 100000

  Output arguments:
    --outdir OUTDIR       If specified all output files will be written to that
                          directory. Default: the current working directory
    -n NAME, --name NAME  Experiment name, which will be used to generate output
                          file names. DEFAULT: "NA"
    -B, --bdg             Whether or not to save extended fragment pileup, and
                          local lambda tracks (two files) at every bp into a
                          bedGraph file. DEFAULT: False
    --verbose VERBOSE     Set verbose level of runtime message. 0: only show
                          critical message, 1: show additional warning message,
                          2: show process information, 3: show debug messages.
                          DEFAULT:2
    --trackline           Tells MACS to include trackline with bedGraph files.
                          To include this trackline while displaying bedGraph at
                          UCSC genome browser, can show name and description of
                          the file as well. However my suggestion is to convert
                          bedGraph to bigWig, then show the smaller and faster
                          binary bigWig file at UCSC genome browser, as well as
                          downstream analysis. Require -B to be set. Default:
                          Not include trackline.
    --SPMR                If True, MACS will save signal per million reads for
                          fragment pileup profiles. Require -B to be set.
                          Default: False

  Shifting model arguments:
    -s TSIZE, --tsize TSIZE
                          Tag size. This will override the auto detected tag
                          size. DEFAULT: Not set
    --bw BW               Band width for picking regions to compute fragment
                          size. This value is only used while building the
                          shifting model. DEFAULT: 300
    -m MFOLD MFOLD, --mfold MFOLD MFOLD
                          Select the regions within MFOLD range of high-
                          confidence enrichment ratio against background to
                          build model. Fold-enrichment in regions must be lower
                          than upper limit, and higher than the lower limit. Use
                          as "-m 10 30". DEFAULT:5 50
    --fix-bimodal         Whether turn on the auto pair model process. If set,
                          when MACS failed to build paired model, it will use
                          the nomodel settings, the --exsize parameter to extend
                          each tags towards 3' direction. Not to use this
                          automate fixation is a default behavior now. DEFAULT:
                          False
    --nomodel             Whether or not to build the shifting model. If True,
                          MACS will not build model. by default it means
                          shifting size = 100, try to set extsize to change it.
                          DEFAULT: False
    --shift SHIFT         (NOT the legacy --shiftsize option!) The arbitrary
                          shift in bp. Use discretion while setting it other
                          than default value. When NOMODEL is set, MACS will use
                          this value to move cutting ends (5') towards 5'->3'
                          direction then apply EXTSIZE to extend them to
                          fragments. When this value is negative, ends will be
                          moved toward 3'->5' direction. Recommended to keep it
                          as default 0 for ChIP-Seq datasets, or -1 * half of
                          EXTSIZE together with EXTSIZE option for detecting
                          enriched cutting loci such as certain DNAseI-Seq
                          datasets. Note, you can't set values other than 0 if
                          format is BAMPE or BEDPE for paired-end data. DEFAULT:
                          0.
    --extsize EXTSIZE     The arbitrary extension size in bp. When nomodel is
                          true, MACS will use this value as fragment size to
                          extend each read towards 3' end, then pile them up.
                          It's exactly twice the number of obsolete SHIFTSIZE.
                          In previous language, each read is moved 5'->3'
                          direction to middle of fragment by 1/2 d, then
                          extended to both direction with 1/2 d. This is
                          equivalent to say each read is extended towards 5'->3'
                          into a d size fragment. DEFAULT: 200. EXTSIZE and
                          SHIFT can be combined when necessary. Check SHIFT
                          option.

  Peak calling arguments:
    -q QVALUE, --qvalue QVALUE
                          Minimum FDR (q-value) cutoff for peak detection.
                          DEFAULT: 0.05. -q, and -p are mutually exclusive.
    -p PVALUE, --pvalue PVALUE
                          Pvalue cutoff for peak detection. DEFAULT: not set.
                          -q, and -p are mutually exclusive. If pvalue cutoff is
                          set, qvalue will not be calculated and reported as -1
                          in the final .xls file.
    --to-large            When set, scale the small sample up to the bigger
                          sample. By default, the bigger dataset will be scaled
                          down towards the smaller dataset, which will lead to
                          smaller p/qvalues and more specific results. Keep in
                          mind that scaling down will bring down background
                          noise more. DEFAULT: False
    --ratio RATIO         When set, use a custom scaling ratio of ChIP/control
                          (e.g. calculated using NCIS) for linear scaling.
                          DEFAULT: ingore
    --down-sample         When set, random sampling method will scale down the
                          bigger sample. By default, MACS uses linear scaling.
                          Warning: This option will make your result unstable
                          and irreproducible since each time, random reads would
                          be selected. Consider to use 'randsample' script
                          instead. <not implmented>If used together with --SPMR,
                          1 million unique reads will be randomly picked.</not
                          implemented> Caution: due to the implementation, the
                          final number of selected reads may not be as you
                          expected! DEFAULT: False
    --seed SEED           Set the random seed while down sampling data. Must be
                          a non-negative integer in order to be effective.
                          DEFAULT: not set
    --tempdir TEMPDIR     Optional directory to store temp files. DEFAULT: /tmp
    --nolambda            If True, MACS will use fixed background lambda as
                          local lambda for every peak region. Normally, MACS
                          calculates a dynamic local lambda to reflect the local
                          bias due to potential chromatin structure.
    --slocal SMALLLOCAL   The small nearby region in basepairs to calculate
                          dynamic lambda. This is used to capture the bias near
                          the peak summit region. Invalid if there is no control
                          data. If you set this to 0, MACS will skip slocal
                          lambda calculation. *Note* that MACS will always
                          perform a d-size local lambda calculation. The final
                          local bias should be the maximum of the lambda value
                          from d, slocal, and llocal size windows. DEFAULT: 1000
    --llocal LARGELOCAL   The large nearby region in basepairs to calculate
                          dynamic lambda. This is used to capture the surround
                          bias. If you set this to 0, MACS will skip llocal
                          lambda calculation. *Note* that MACS will always
                          perform a d-size local lambda calculation. The final
                          local bias should be the maximum of the lambda value
                          from d, slocal, and llocal size windows. DEFAULT:
                          10000.
    --broad               If set, MACS will try to call broad peaks by linking
                          nearby highly enriched regions. The linking region is
                          controlled by another cutoff through --linking-cutoff.
                          The maximum linking region length is 4 times of d from
                          MACS. DEFAULT: False
    --broad-cutoff BROADCUTOFF
                          Cutoff for broad region. This option is not available
                          unless --broad is set. If -p is set, this is a pvalue
                          cutoff, otherwise, it's a qvalue cutoff. DEFAULT: 0.1
    --cutoff-analysis     While set, MACS2 will analyze number or total length
                          of peaks that can be called by different p-value
                          cutoff then output a summary table to help user decide
                          a better cutoff. The table will be saved in
                          NAME_cutoff_analysis.txt file. Note, minlen and maxgap
                          may affect the results. WARNING: May take ~30 folds
                          longer time to finish. DEFAULT: False

  Post-processing options:
    --call-summits        If set, MACS will use a more sophisticated signal
                          processing approach to find subpeak summits in each
                          enriched peak region. DEFAULT: False
    --fe-cutoff FECUTOFF  When set, the value will be used to filter out peaks
                          with low fold-enrichment. Note, MACS2 use 1.0 as
                          pseudocount while calculating fold-enrichment.
                          DEFAULT: 1.0

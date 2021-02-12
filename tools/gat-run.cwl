cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() {
          if (inputs.output_filename == ""){
            var root = inputs.segment_file.basename.split('.').slice(0,-1).join('.');
            return (root == "")?inputs.segment_file.basename+".tsv":root+".tsv";
          } else {
            return inputs.output_filename;
          }
        };


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/gat:v0.0.1


inputs:

  segment_file:
    type: File
    inputBinding:
      position: 5
      prefix: "-s"
    doc: |
      BED file (strictly 3 columns) with sets of intervals whose association is tested with annotation_file

  annotation_file:
    type: File
    inputBinding:
      position: 6
      prefix: "-a"
    doc: |
      BED file (strictly 3 columns) with sets of intervals that are used for testing association of segment_file

  workspace_file:
    type: File
    inputBinding:
      position: 7
      prefix: "-w"
    doc: |
      BED file (strictly 3 columns) with genomic regions accessible for simulation

  output_filename:
    type: string?
    inputBinding:
      position: 8
      valueFrom: $(default_output_filename())
      prefix: "-S"
    default: ""
    doc: |
      Output report file name

  iterations:
    type: int?
    inputBinding:
      position: 9
      prefix: "-n"
    doc: |
      Number of iterations

  counter:
    type:
    - "null"
    - type: enum
      name: "counter"
      symbols:
      - "nucleotide-overlap"
      - "nucleotide-density"
      - "segment-overlap"
      - "annotation-overlap"
      - "segment-midoverlap"
      - "annotation-midoverlap"
    inputBinding:
      position: 10
      prefix: "-c"
    doc: |
      Set the measure of association that is tested.
      Default: nucleotide-overlap

  threads:
    type: int?
    inputBinding:
      position: 11
      prefix: "-t"
    doc: |
      Threads number

  seed:
    type: int?
    inputBinding:
      position: 12
      prefix: "--random-seed="
      separate: false
    doc: |
      Random seed to initialize random number generator with


outputs:

  report_file:
    type: File
    outputBinding:
      glob: $(default_output_filename())
    doc: |
      Report file

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["gat-run.py", "--ignore-segment-tracks"]
stdout: gat_stdout.log
stderr: gat_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


s:name: "gat-run"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/gat-run.cwl
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
  A common question in genomic analysis is whether two sets of genomic intervals overlap significantly.
  This question arises, for example, in the interpretation of ChIP-Seq or RNA-Seq data. The Genomic
  Association Tester (GAT) is a tool for computing the significance of overlap between multiple sets of
  genomic intervals. GAT estimates significance based on simulation.

  Gat implemements a sampling algorithm. Given a chromosome (workspace) and segments of interest, for
  example from a ChIP-Seq experiment, gat creates randomized version of the segments of interest falling
  into the workspace. These sampled segments are then compared to existing genomic annotations.

  Note:
  --ignore-segment-tracks parameter is hardcoded

s:about: |
  Usage
  -----

  Example::

    python gat-run.py
        --segment-file=segments.bed.gz
        --workspace-file=workspace.bed.gz
        --annotation-file=annotations_architecture.bed.gz

  Type::

    python gat-run.py --help

  for command line help.

  Documentation
  -------------

  Code
  ----



  Options:
    --version             show program's version number and exit
    -h, --help            show this help message and exit

    Input options:
      -a ANNOTATION_FILES, --annotation-bed-file=ANNOTATION_FILES, --annotations=ANNOTATION_FILES, --annotation-file=ANNOTATION_FILES
                          filename with annotations [default=[]].
      -s SEGMENT_FILES, --segment-bed-file=SEGMENT_FILES, --segments=SEGMENT_FILES, --segment-file=SEGMENT_FILES
                          filename with segments. Also accepts a glob in
                          parentheses [default=[]].
      -w WORKSPACE_FILES, --workspace-bed-file=WORKSPACE_FILES, --workspace=WORKSPACE_FILES, --workspace-file=WORKSPACE_FILES
                          filename with workspace segments. Also accepts a glob
                          in parentheses [default=[]].
      -i ISOCHORE_FILES, --isochore-bed-file=ISOCHORE_FILES, --isochores=ISOCHORE_FILES, --isochore-file=ISOCHORE_FILES
                          filename with isochore segments. Also accepts a glob
                          in parentheses [default=none].
      -l SAMPLE_FILES, --sample-file=SAMPLE_FILES
                          filename with sample files. Start processing from
                          samples [default=[]].
      --input-counts-file=INPUT_FILENAME_COUNTS
                          start processing from counts - no segments required
                          [default=none].
      --input-results-file=INPUT_FILENAME_RESULTS
                          start processing from results - no segments required
                          [default=none].
      --ignore-segment-tracks
                          ignore segment tracks - all segments belong to one
                          track and called 'merged' [default=True]
      --with-segment-tracks
                          the segments data file is arranged in tracks.
                          [default=True]
      --enable-split-tracks
                          permit the same track to be in multiple files
                          [default=False]
      --overlapping-annotations
                          the annotations within a track are overlapping and
                          should not be merged. This is useful for working with
                          short-read data. [default=default]
      --annotations-label=ANNOTATIONS_LABEL
                          ignore tracks in annotations and instead set them to
                          label [default=default]
      --annotations-to-points=ANNOTATIONS_TO_POINTS
                          convert annotations from segments to positions.
                          Available methods are 'midpoint', 'start' or 'end'.
                          [default=default]

    Output options:
      -o OUTPUT_ORDER, --order=OUTPUT_ORDER
                          order results in output by fold, track, etc.
                          [default=fold].
      --output-tables-pattern=OUTPUT_TABLES_PATTERN
                          output pattern for result tables. Used if there are
                          multiple counters used [default=%s.tsv.gz].
      --output-counts-pattern=OUTPUT_COUNTS_PATTERN
                          output pattern for counts [default=none].
      --output-plots-pattern=OUTPUT_PLOTS_PATTERN
                          output pattern for plots [default=none]
      --output-samples-pattern=OUTPUT_SAMPLES_PATTERN
                          output pattern for samples. Samples are stored in bed
                          format, one for  each segment [default=none]
      --output-stats=OUTPUT_STATS
                          output overlap summary stats [default=[]].
      --output-bed=OUTPUT_BED
                          output bed files [default=[]].
      --descriptions=INPUT_FILENAME_DESCRIPTIONS
                          filename mapping annotation terms to descriptions. if
                          given, the output table will contain additional
                          columns [default=none]

    Sampling algorithm options:
      -c COUNTERS, --counter=COUNTERS
                          quantity to use for estimating enrichment
                          [default=[]].
      -m SAMPLER, --sampler=SAMPLER
                          quantity to test [default=annotator].
      -n NUM_SAMPLES, --num-samples=NUM_SAMPLES
                          number of samples to compute [default=1000].
      --shift-extension=SHIFT_EXTENSION
                          if the sampling method is 'shift', create a segment of
                          size # anound the segment to determine the size of the
                          region for shifthing [default=0].
      --shift-expansion=SHIFT_EXPANSION
                          if the sampling method is 'shift', multiply each
                          segment by # to determine the size of the region for
                          shifthing [default=2.0].
      --bucket-size=BUCKET_SIZE
                          size of a bin for histogram of segment lengths. If 0,
                          it will be automatically scaled to fit nbuckets
                          [default=0]
      --nbuckets=NBUCKETS
                          number of bins for histogram of segment lengths
                          [default=100000]

    Statistics options:
      -p PVALUE_METHOD, --pvalue-method=PVALUE_METHOD
                          type of pvalue reported [default=empirical].
      -q QVALUE_METHOD, --qvalue-method=QVALUE_METHOD
                          method to perform multiple testing correction by
                          controlling the fdr [default=BH].
      --qvalue-lambda=QVALUE_LAMBDA
                          fdr computation: lambda [default=none].
      --qvalue-pi0-method=QVALUE_PI0_METHOD
                          fdr computation: method for estimating pi0
                          [default=smoother].
      --pseudo-count=PSEUDO_COUNT
                          pseudo count. The pseudo count is added to both the
                          observed and expected overlap. Using a pseudo-count
                          avoids gat reporting fold changes of 0 [default=1.0].
      --null=NULL         null hypothesis. The default is to test categories for
                          enrichment/depletion. If a filename with gat output is
                          given, gat will test for the difference in fold change
                          between the segments supplied and in the other file
                          [default=default].

    Processing options:
      -e CACHE, --cache=CACHE
                          filename for caching samples [default=none].
      -t NUM_THREADS, --num-threads=NUM_THREADS
                          number of threads to use for sampling [default=0]
      --random-seed=RANDOM_SEED
                          random seed to initialize number generator with
                          [none].

    Workspace manipulation (experimental):
      --conditional=CONDITIONAL
                          conditional workspace creation
                          [default=unconditional]*cooccurance* - compute
                          enrichment only within workspace segments that contain
                          both segments and annotations, *annotation-centered* -
                          workspace centered around annotations. See
                          --conditional-extension,segment-centered - workspace
                          centered around segments. See --conditional-extension
      --conditional-extension=CONDITIONAL_EXTENSION
                          if workspace is created conditional, extend by this
                          amount (in bp) [default=none].
      --conditional-expansion=CONDITIONAL_EXPANSION
                          if workspace is created conditional, expand by this
                          amount (ratio) [default=none].
      --restrict-workspace
                          restrict workspace to those segments that contain both
                          track and annotations [default=False]
      --truncate-workspace-to-annotations
                          truncate workspace with annotations [default=False]
      --truncate-segments-to-workspace
                          truncate segments to workspace before sampling
                          [default=False]

    Script timing options:
      --timeit=TIMEIT_FILE
                          store timeing information in file [none].
      --timeit-name=TIMEIT_NAME
                          name in timing file for this class of jobs [all].
      --timeit-header     add header for timing information [none].

    Common options:
      -v LOGLEVEL, --verbose=LOGLEVEL
                          loglevel [1]. The higher, the more output.
      -?                  output short help (command line options only.

    Input/output options:
      -P OUTPUT_FILENAME_PATTERN, --output-filename-pattern=OUTPUT_FILENAME_PATTERN
                          OUTPUT filename pattern for various methods [%s].
      -F, --force-output  force over-writing of existing files.
      -I FILE, --stdin=FILE
                          file to read stdin from [default = stdin].
      -L FILE, --log=FILE
                          file with logging information [default = stdout].
      -E FILE, --error=FILE
                          file with error information [default = stderr].
      -S FILE, --stdout=FILE
                          file where output is to go [default = stdout].
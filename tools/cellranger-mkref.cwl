cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var get_output_folder_name = function() {
          if (inputs.output_folder_name == ""){
            var root = inputs.genome_fasta_file.basename.split('.').slice(0,-1).join('.');
            return (root == "")?inputs.genome_fasta_file.basename:root;
          } else {
            return inputs.output_folder_name;
          }
        };

hints:
- class: DockerRequirement
  dockerPull: cumulusprod/cellranger:8.0.1


inputs:
  
  genome_fasta_file:
    type: File
    inputBinding:
      position: 5
      prefix: "--fasta"
    doc: |
      Genome FASTA file

  annotation_gtf_file:
    type: File
    inputBinding:
      position: 6
      prefix: "--genes"
    doc: |
      GTF annotation file

  output_folder_name:
    type: string?
    inputBinding:
      position: 7
      prefix: "--genome"
      valueFrom: $(get_output_folder_name())
    default: ""
    doc: |
      Unique genome name, used to name output folder

  threads:
    type: int?
    inputBinding:
      valueFrom: $(["--nthreads", self, "--localcores", self])
      position: 8
    doc: |
      Number of threads used during STAR
      genome index. And the max cores the
      pipeline may request at one time.
      Default: 1 for --nthreads and all
      available for --localcores

  memory_limit:
    type: int?
    inputBinding:
      position: 9
      prefix: "--memgb"
    doc: |
      Maximum memory (GB) used when aligning reads with STAR
      Defaults: 16


outputs:

  indices_folder:
    type: Directory
    outputBinding:
      glob: $(get_output_folder_name())
    doc: |
      Cellranger-compatible reference folder that includes
      STAR indices and some additional files

  chrom_length_file:
    type: File
    outputBinding:
      glob: $(get_output_folder_name() + "/star/chrNameLength.txt")
    doc: |
      Chromosome length file in TSV format

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["cellranger", "mkref"]


stdout: cellranger_mkref_stdout.log
stderr: cellranger_mkref_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Cell Ranger Reference (RNA)"
s:name: "Cell Ranger Reference (RNA)"
s:alternateName: "Builds a reference genome of a selected species for quantifying gene expression"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/cellranger-mkref.cwl
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
  Cell Ranger Reference (RNA)

  Builds a reference genome of a selected species
  for quantifying gene expression.

  Both --nthreads and --localcores parameters are
  configured through "threads" input.


s:about: |
  Prepare a reference for use with 10x analysis software. Requires a GTF and FASTA

  Usage: cellranger mkref [OPTIONS] --genome <GENOME_NAMES> --fasta <FASTA_FILES> --genes <GTF_FILES>

  Options:
        --genome <GENOME_NAMES>     Unique genome name, used to name output folder [a-zA-Z0-9_-]+. Specify
                                    multiple genomes by specifying this argument multiple times; the output
                                    folder will be <name1>_and_<name2>
        --fasta <FASTA_FILES>       Path to FASTA file containing your genome reference. Specify multiple
                                    genomes by specifying this argument multiple times
        --genes <GTF_FILES>         Path to genes GTF file containing annotated genes for your genome
                                    reference. Specify multiple genomes by specifying this argument multiple
                                    times
        --nthreads <NUM_THREADS>    Number of threads used during STAR genome index generation. Defaults to 1
                                    [default: 1]
        --memgb <MEM_GB>            Maximum memory (GB) used [default: 16]
        --ref-version <REF_VERSION> Optional reference version string to include with reference
        --dry                       Do not execute the pipeline. Generate a pipeline invocation (.mro) file
                                    and stop
        --jobmode <MODE>            Job manager to use. Valid options: local (default), sge, lsf, slurm or
                                    path to a .template file. Search for help on "Cluster Mode" at
                                    support.10xgenomics.com for more details on configuring the pipeline to
                                    use a compute cluster
        --localcores <NUM>          Set max cores the pipeline may request at one time. Only applies to local
                                    jobs
        --localmem <NUM>            Set max GB the pipeline may request at one time. Only applies to local
                                    jobs
        --localvmem <NUM>           Set max virtual address space in GB for the pipeline. Only applies to
                                    local jobs
        --mempercore <NUM>          Reserve enough threads for each job to ensure enough memory will be
                                    available, assuming each core on your cluster has at least this much
                                    memory available. Only applies to cluster jobmodes
        --maxjobs <NUM>             Set max jobs submitted to cluster at one time. Only applies to cluster
                                    jobmodes
        --jobinterval <NUM>         Set delay between submitting jobs to cluster, in ms. Only applies to
                                    cluster jobmodes
        --overrides <PATH>          The path to a JSON file that specifies stage-level overrides for cores
                                    and memory. Finer-grained than --localcores, --mempercore and --localmem.
                                    Consult https://support.10xgenomics.com/ for an example override file
        --output-dir <PATH>         Output the results to this directory
        --uiport <PORT>             Serve web UI at http://localhost:PORT
        --disable-ui                Do not serve the web UI
        --noexit                    Keep web UI running after pipestance completes or fails
        --nopreflight               Skip preflight checks
    -h, --help                      Print help
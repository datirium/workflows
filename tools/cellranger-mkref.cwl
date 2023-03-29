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
  dockerPull: cumulusprod/cellranger:7.0.0


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
      position: 8
      prefix: "--nthreads"
    doc: |
      Number of threads used during STAR genome index
      Default: 1

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

label: "Cell Ranger Build Reference Indices"
s:name: "Cell Ranger Build Reference Indices"
s:alternateName: "Builds Cell Ranger compatible reference folder from the custom genome FASTA and gene GTF annotation files"

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
  Cell Ranger Build Reference Indices

  Builds Cell Ranger compatible reference folder from
  the custom genome FASTA and gene GTF annotation files.


s:about: |
  Build a Cell Ranger-compatible reference folder from user-supplied genome
  FASTA and gene GTF files. Creates a new folder named after the genome.

  The commands below should be preceded by 'cellranger':

  Usage:
      mkref
          --genome=NAME ...
          --fasta=PATH ...
          --genes=PATH ...
          [options]
      mkref -h | --help | --version

  Arguments:
      genome              Unique genome name(s), used to name output folder
                              [a-zA-Z0-9_-]+. Specify multiple genomes by
                              specifying the --genome argument multiple times; the
                              output folder will be <name1>_and_<name2>.
      fasta               Path(s) to FASTA file containing your genome reference.
                              Specify multiple genomes by specifying the --fasta
                              argument multiple times.
      genes               Path(s) to genes GTF file(S) containing annotated genes
                              for your genome reference. Specify multiple genomes
                              by specifying the --genes argument multiple times.

  Options:
      --nthreads=<num>    Number of threads used during STAR genome index
                              generation. Defaults to 1.
      --memgb=<num>       Maximum memory (GB) used when aligning reads with STAR.
                              Defaults to 16.
      --ref-version=<str> Optional reference version string to include with
                              reference.
      -h --help           Show this message.
      --version           Show version.
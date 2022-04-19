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
- class: InitialWorkDirRequirement
  listing: |
    ${
      var exclude_chr = "[]";
      if (inputs.exclude_chr && inputs.exclude_chr.length > 0){
        exclude_chr = '["' + inputs.exclude_chr.join('", "') + '"]'
      }
      var entry = `
      {
          genome: ["${get_output_folder_name()}"]
          input_fasta: ["${inputs.genome_fasta_file.path}"]
          input_gtf: ["${inputs.annotation_gtf_file.path}"]
          non_nuclear_contigs: ${exclude_chr}
      }`
      return [{
        "entry": entry,
        "entryname": "config.txt"
      }];
    }


hints:
- class: DockerRequirement
  dockerPull: cumulusprod/cellranger-arc:2.0.0


inputs:
  
  genome_fasta_file:
    type: File
    doc: |
      Genome FASTA file

  annotation_gtf_file:
    type: File
    doc: |
      GTF annotation file.

  exclude_chr:
    type:
      - "null"
      - string[]
    doc: |
      Contigs that do not have any chromatin structure, for example,
      mitochondria or plastids. These contigs are excluded from peak
      calling since the entire contig will be "open" due to a lack of
      chromatin structure

  output_folder_name:
    type: string?
    default: ""
    doc: |
      Unique genome name, used to name output folder

  threads:
    type: int?
    inputBinding:
      position: 5
      prefix: "--nthreads"
    doc: |
      Number of threads used during STAR genome indexing
      Default: 1

  memory_limit:
    type: int?
    inputBinding:
      position: 6
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
      Compatible with Cell Ranger ARC reference folder that includes
      STAR and BWA indices

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["cellranger-arc", "mkref", "--config", "config.txt"]


stdout: cellranger_arc_mkref_stdout.log
stderr: cellranger_arc_mkref_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Cell Ranger ARC Build Reference Indices"
s:name: "Cell Ranger ARC Build Reference Indices"
s:alternateName: "Builds Cell Ranger ARC compatible reference folder from the custom genome FASTA and gene GTF annotation files"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/cellranger-arc-mkref.cwl
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
  Cell Ranger ARC Build Reference Indices
  ====================================================================
  
  Builds Cell Ranger ARC compatible reference folder from the custom
  genome FASTA and gene GTF annotation files.
  
  Notes:
  - `input_motifs` parameter in the `config.txt` file is not
    implemented.
  - if GTF file provided in `annotation_gtf_file` has records with
    duplicate gene_id, they should be grouped together. Applicable to
    USCS RefGene annotations.


s:about: |
  Reference preparation tool for 10x Genomics Cell Ranger Multiome ATAC + Gene Expression.

  Build a reference package from a user-supplied genome FASTA and gene GTF file.
  Creates a new folder named after the genome.

  NOTE: Multi-species references are not supported by cellranger-arc. If you
  construct a multi-species reference and run 'cellranger-arc count' you will not
  be able to generate all the outputs of the pipeline.

  The commands below should be preceded by 'cellranger-arc':

  Usage:
      mkref
          --config=PATH
          [options]
      mkref -h | --help | --version

  Arguments:
      config              Path to configuration file containing additional
                              information about the reference. See online
                              documentation for more details. The following is an
                              example of a config file:
                              {
                                  organism: "human"
                                  genome: ["GRCh38"]
                                  input_fasta: ["/path/to/GRCh38/assembly.fa"]
                                  input_gtf: ["/path/to/gencode/annotation.gtf"]
                                  non_nuclear_contigs: ["chrM"]
                                  input_motifs: "/path/to/jaspar/motifs.pfm"
                              }
                              Parameters:
                                  - organism: (optional; string) name of the
                                      organism
                                  - genome: (required; list of strings) name(s) of
                                      the genome(s) that comprise the organism
                                  - input_fasta: (required; list of paths) path(s)
                                      to the assembly fasta file(s) for each
                                      genome
                                  - input_gtf: (required; list of paths) path(s)
                                      to the gene annotation GTF file(s) for each
                                      genome
                                  - non_nuclear_contigs:
                                      (optional; list of strings) contigs in the
                                      assembly that are not nuclear and have no
                                      chromatin structure (e.g., mitochondria)
                                  - input_motifs: (optional; path) path to a
                                      motif annotations file in the JASPAR format
                              The above config file would create a reference
                              package in "$(pwd)/GRCh38".

  Options:
      --nthreads=<num>    Number of threads used during STAR genome index
                              generation. Defaults to 1.
      --memgb=<num>       Maximum memory (GB) used when aligning reads with STAR.
                              Defaults to 16.
      --ref-version=<str> Optional reference version string to include with
                              reference.
      -h --help           Show this message.
      --version           Show version.

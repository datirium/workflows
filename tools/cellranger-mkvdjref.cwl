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
      Genome FASTA file. Hard/soft-masked files are not allowed.

  annotation_gtf_file:
    type: File
    inputBinding:
      position: 6
      prefix: "--genes"
    doc: |
      GTF annotation file. Should include gene_biotype/transcript_biotype fields

  output_folder_name:
    type: string?
    inputBinding:
      position: 7
      prefix: "--genome"
      valueFrom: $(get_output_folder_name())
    default: ""
    doc: |
      Unique genome name, used to name output folder


outputs:

  indices_folder:
    type: Directory
    outputBinding:
      glob: $(get_output_folder_name())
    doc: |
      Cell Ranger V(D)J-compatible reference folder.
      This folder will include V(D)J segment FASTA file.

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["cellranger", "mkvdjref"]


stdout: cellranger_mkvdjref_stdout.log
stderr: cellranger_mkvdjref_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Cell Ranger Build V(D)J Reference Indices"
s:name: "Cell Ranger Build V(D)J Reference Indices"
s:alternateName: "Build a Cell Ranger V(D)J-compatible reference folder from a user-supplied genome FASTA and gene GTF files"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/cellranger-mkvdjref.cwl
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
  Cell Ranger Build V(D)J Reference Indices
  
  Build a Cell Ranger V(D)J-compatible reference folder from:
  1) A user-supplied genome FASTA and gene GTF files.
        For example, using files from ENSEMBL.
  2) A FASTA file containing V(D)J segments as per the mkvdjref spec.
        For example, using files from IMGT.

  For simplicity purpose only option 1) is supported - user need to
  provide GTF annotation file, input --seqs is not implemented.

  Chromosome names in GTF file should correspond to the chromosome
  names in FASTA file.


s:about: |
  Reference preparation tool for 10x Genomics Cell Ranger V(D)J assembler.

  Build a Cell Ranger V(D)J-compatible reference folder from:
  1) A user-supplied genome FASTA and gene GTF files.
        For example, using files from ENSEMBL.
      OR
  2) A FASTA file containing V(D)J segments as per the mkvdjref spec.
        For example, using files from IMGT.

  Creates a new folder named after the genome.

  The commands below should be preceded by 'cellranger':

  Usage:
      mkvdjref --genome=NAME --fasta=PATH --genes=PATH ...[options]
      mkvdjref --genome=NAME --seqs=PATH [options]
      mkvdjref -h | --help | --version

  Arguments:
      genome              A unique genome name, used to name output folder
                              [a-zA-Z0-9_-]+.
      fasta               Path to FASTA file containing your genome reference.
      genes               One or more GTF files containing annotated genes for
                              your genome reference. Specify multiple files by
                              specifying the --genes argument multiple times. The
                              files will be concatenated.
      seqs                A FASTA file that directly specifies V(D)J sequences.
                              This is mutually exclusive with the the "fasta" and
                              "genes" args above.

  Options:
      --ref-version=<str>
                          Optional reference version string to include.
      --rm-transcripts=PATH
                          Path to text file with transcript IDs to ignore. This
                              file should have one transcript ID per line where
                              the IDs correspond to the "transcript_id" key in the
                              GTF info column.
      -h --help           Show this message.
      --version           Show version.

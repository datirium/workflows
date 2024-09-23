cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: robertplayer/scidap-tgif:v1.0.0


inputs:

  insertions_filtered:
    type: File
    inputBinding:
      prefix: "-i"
    doc: "FASTQ file of WGS or NCATS enriched library (ideally long read, short read single-end only)"

  reference_fasta:
    type: File
    inputBinding:
      prefix: "-r"
    doc: "reference fasta file for mapping"

  threads_count:
    type: int?
    inputBinding:
      prefix: "-t"
    doc: "Threads number"


outputs:

  primer3_output:
    type: File?
    outputBinding:
      glob: "primer3.tar"

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["/TgIF/tgif_primer3.sh"]

stdout: tgif-primer3_stdout.log
stderr: tgif-primer3_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "tgif-primer3"
s:downloadUrl: https://github.com/datirium/workflows/blob/master/tools/tgif-primer3.cwl
s:codeRepository: https://github.com/datirium/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Datirium LLC"
  s:location:
  - class: s:PostalAddress
    s:addressCountry: "USA"
    s:addressLocality: "Cincinnati"
    s:addressRegion: "OH"
    s:postalCode: ""
    s:streetAddress: ""
    s:telephone: ""
  s:logo: "https://avatars.githubusercontent.com/u/33202955?s=200&v=4"
  s:department:
  - class: s:Organization
    s:legalName: "Datirium LLC"
    s:department:
    - class: s:Organization
      s:legalName: "Bioinformatics"
      s:member:
      - class: s:Person
        s:name: Robert Player
        s:email: mailto:support@datirium.com
        s:sameAs:
        - id: https://orcid.org/0000-0001-5872-259X

      
doc: |
  TgIF Primer3 module
  ==============================================

  This tgif module is used to design primers for validation of potential insertion sites identified by the primary tgif algorithm (tgif_ncats.sh).

  DEPENDENCIES:
          GNU Parallel
                  O. Tange (2011): GNU Parallel - The Command-Line Power Tool,
                  ;login: The USENIX Magazine, February 2011:42-47.
          Primer3
                  please symlink 'primer3_core' to 'tgif/bin/' (see Installation section of README)
                  http://primer3.org/manual.html

  HELP/OUTFMT:
          -h      help    show this message
          -f        format        format of input
  REQUIRED:
          -t      INT     number of threads to GNU parallel over
          -i      TGIF    full path to output from primary tgif algorithm (insertions_filtered.tgif)
          -r      FASTA   fasta reference of insert target organism (used to generate -i)
  OPTIONAL:
          -g      INT             if GAP_LENGTH of insertion site is greater than INT bp (-g), no primers will be generated for the site [default 5000]

  USAGE:
  i="/data/project/tgif_ncats-reads.fastq/insertions_filtered.tgif"
  r="/data/project/org_reference.fna"
  bash tgif_primer3.sh -t 10 -i "$i" -r "$r"
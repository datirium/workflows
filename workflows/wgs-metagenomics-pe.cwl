cwlVersion: v1.0
class: Workflow


requirements:
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement


'sd:upstream':
  database: "kraken2-databases.cwl"


inputs:

  alias:
    type: string
    label: "Sample short name/Alias:"
    sd:preview:
      position: 1

  condition:
    type: string?
    label: "Experimental condition:"
    sd:preview:
      position: 2

  fastq_file_R1:
    type:
      - File
      - type: array
        items: File
    label: "Read 1 file:"
    'sd:localLabel': true
    format: "http://edamontology.org/format_1930"
    doc: "Read 1 FASTQ file from a paired-end sequencing run"
    sd:preview:
      position: 4

  fastq_file_R2:
    type:
      - File
      - type: array
        items: File
    label: "Read 2 file:"
    'sd:localLabel': true
    format: "http://edamontology.org/format_1930"
    doc: "Read 2 FASTQ file that pairs with the input R1 file"
    sd:preview:
      position: 5

  k2db:
    type: Directory
    'sd:upstreamSource': "database/k2db"
    label: "Reference genome database for metagenomic sequence classification:"
    'sd:localLabel': true
    doc: "Pre-built kraken2 reference genome database for taxonomic classification of sequencing reads. A 'database' sample needs to be added to your project that will populate this dropdown."
    sd:preview:
      position: 3


outputs:

  fastx_statistics_upstream:
    type: File
    label: "FASTQ 1 statistics"
    format: "http://edamontology.org/format_2330"
    doc: "fastx_quality_stats generated FASTQ 1 quality statistics file"
    outputSource: fastx_quality_stats_upstream/statistics_file
    'sd:visualPlugins':
    - line:
        tab: 'QC Plots'
        Title: 'FASTQ 1 Base frequency plot'
        xAxisTitle: 'Nucleotide position'
        yAxisTitle: 'Frequency'
        colors: ["#b3de69", "#888888", "#fb8072", "#fdc381", "#99c0db"]
        data: [$13, $14, $15, $16, $17]
    - boxplot:
        tab: 'QC Plots'
        Title: 'FASTQ 1 Quality Control'
        xAxisTitle: 'Nucleotide position'
        yAxisTitle: 'Quality score'
        colors: ["#b3de69", "#888888", "#fb8072", "#fdc381", "#99c0db"]
        data: [$11, $7, $8, $9, $12]

  fastx_statistics_downstream:
    type: File
    label: "FASTQ 2 statistics"
    format: "http://edamontology.org/format_2330"
    doc: "fastx_quality_stats generated FASTQ 2 quality statistics file"
    outputSource: fastx_quality_stats_downstream/statistics_file
    'sd:visualPlugins':
    - line:
        tab: 'QC Plots'
        Title: 'FASTQ 2 Base frequency plot'
        xAxisTitle: 'Nucleotide position'
        yAxisTitle: 'Frequency'
        colors: ["#b3de69", "#888888", "#fb8072", "#fdc381", "#99c0db"]
        data: [$13, $14, $15, $16, $17]
    - boxplot:
        tab: 'QC Plots'
        Title: 'FASTQ 2 Quality Control'
        xAxisTitle: 'Nucleotide position'
        yAxisTitle: 'Quality score'
        colors: ["#b3de69", "#888888", "#fb8072", "#fdc381", "#99c0db"]
        data: [$11, $7, $8, $9, $12]

  trim_report_upstream:
    type: File
    label: "TrimGalore report Upstream"
    doc: "TrimGalore generated log for FASTQ 1"
    outputSource: trim_fastq/report_file

  trim_report_downstream:
    type: File
    label: "TrimGalore report Downstream"
    doc: "TrimGalore generated log for FASTQ 2"
    outputSource: trim_fastq/report_file_pair

  decontaminated_kneaddata_reads_R1:
    type:
      - "null"
      - File
    format: "http://edamontology.org/format_1930"
    label: "decontaminated reads using kneaddata r1 FASTQ file"
    doc: "decontaminated reads using kneaddata r1 FASTQ file"
    outputSource: decontaminate_with_kneaddata/kneaddata_cleaned_R1

  decontaminated_kneaddata_reads_R2:
    type:
      - "null"
      - File
    format: "http://edamontology.org/format_1930"
    label: "decontaminated reads using kneaddata r2 FASTQ file"
    doc: "decontaminated reads using kneaddata r2 FASTQ file"
    outputSource: decontaminate_with_kneaddata/kneaddata_cleaned_R2

  decontaminate_with_kneaddata_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "tool log file"
    doc: "captures kneaddata command log file"
    outputSource: decontaminate_with_kneaddata/kneaddata_log

  decontaminate_with_kneaddata_stdout:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "stdout logfile"
    doc: "captures standard output from wgs-kneaddata-pe.cwl"
    outputSource: decontaminate_with_kneaddata/stdout_log

  decontaminate_with_kneaddata_stderr:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "stderr logfile"
    doc: "captures standard error from wgs-kneaddata-pe.cwl"
    outputSource: decontaminate_with_kneaddata/stderr_log

  k2_output:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "kraken2 raw output file"
    doc: "raw per read taxonomic classifications from kraken2"
    outputSource: kraken2_classify/k2_output

  k2_report_file:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "kraken2 report file"
    doc: "summary of all read taxonomic classifications from kraken2"
    outputSource: kraken2_classify/k2_report

  k2_report_tsv:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "kraken2 report file tsv"
    doc: "summary of all read taxonomic classifications from kraken2 formatted as a tsv"
    outputSource: kraken2_classify/k2_report_tsv
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Kraken Report'
        Title: 'Summary of Taxonomic Classifications from Kraken2'

  k2_stderr:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "parsed k2 stderr"
    doc: "markdown parsed standard error captured directly from kraken2 classify command in k2-classify-pe.cwl"
    outputSource: kraken2_classify/k2_stderr

  kraken2_classify_stdout:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "stdout logfile"
    doc: "captures standard output from k2-classify-pe.cwl"
    outputSource: kraken2_classify/stdout_log

  kraken2_classify_stderr:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "stderr logfile"
    doc: "captures standard error from k2-classify-pe.cwl"
    outputSource: kraken2_classify/stderr_log

  krona_plot_link:
    type: File
    format: "http://edamontology.org/format_2331"
    label: "Krona plot - hierarchical visualization of taxonomic classifications"
    doc: "hierarchical visualization of taxonomic classifications"
    outputSource: kraken2_classify/krona_html
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  metaphlan_k2_unclassified_reads_profile:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "metagenomic abundance profile"
    doc: "metagenomic abundance profile of reads that were left unclassified by kraken2 step"
    outputSource: classify_unclassified_k2_reads_with_metaphlan/abundance_profile

  metaphlan_cleaned_reads_profile:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "metagenomic abundance profile"
    doc: "metagenomic abundance profile of kneaddata decontaminated reads"
    outputSource: classify_cleaned_reads_with_metaphlan/abundance_profile

  metaphlan_cleaned_reads_profile_scidap_tab:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "metagenomic abundance profile"
    doc: "metagenomic abundance profile of kneaddata decontaminated reads headers cleaned for scidap output tab"
    outputSource: classify_cleaned_reads_with_metaphlan/abundance_profile_scidap
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'MetaPhlAn Profile'
        Title: 'Summary of Taxonomic Classifications from MetaPhlAn'

  metaphlan_cleaned_reads_table:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "metagenomic abundance table"
    doc: "metagenomic abundance table of kneaddata decontaminated reads"
    outputSource: classify_cleaned_reads_with_metaphlan/abundance_table

  metaphlan_cleaned_reads_table_species:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "metagenomic abundance table at species level only"
    doc: "metagenomic abundance table at species level only of kneaddata decontaminated reads"
    outputSource: classify_cleaned_reads_with_metaphlan/abundance_table_species

  metaphlan_cleaned_reads_stdout:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "stdout logfile"
    doc: "captures stdout from wgs-metaphlan-pe.cwl"
    outputSource: classify_cleaned_reads_with_metaphlan/stdout_log

  metaphlan_cleaned_reads_stderr:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "stderr logfile"
    doc: "captures stderr from wgs-metaphlan-pe.cwl"
    outputSource: classify_cleaned_reads_with_metaphlan/stderr_log

  humann_cleaned_reads_genefamilies_rpk:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "gene families profile in rpk units"
    doc: "gene families profile of kneaddata decontaminated reads in rpk units"
    outputSource: functional_assignment_with_humann/genefamilies_rpk

  humann_cleaned_reads_genefamilies_cpm:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "gene families profile in cpm units"
    doc: "gene families profile of kneaddata decontaminated reads in cpm units"
    outputSource: functional_assignment_with_humann/genefamilies_cpm

  humann_cleaned_reads_regroup_rxn_cpm:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "regrouped genes to other functional categories in rxn database in cpm units"
    doc: "regrouped genes to other functional categories in rxn database in cpm units"
    outputSource: functional_assignment_with_humann/regroup_to_rxn_cpm

  humann_cleaned_reads_regroup_rxn_cpm_named:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "regrouped genes to other functional categories in rxn database in cpm units with human-readable names"
    doc: "regrouped genes to other functional categories in rxn database in cpm units with human-readable names"
    outputSource: functional_assignment_with_humann/regroup_to_rxn_cpm_named
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'HUMAnN Functional Catagories'
        Title: 'Functional Catagories from HUMAnN'

  humannn_cleaned_reads_stdout_log:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "stdout logfile"
    doc: "captures stdout from wgs-humann-pe.cwl"
    outputSource: functional_assignment_with_humann/stdout_log

  humann_cleaned_reads_stderr_log:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "stderr logfile"
    doc: "captures stderr from wgs-humann-pe.cwl"
    outputSource: functional_assignment_with_humann/stderr_log


steps:

  extract_fastq_R1:
    label: "Loading unmapped sequence data for read 1"
    doc: |
      Most DNA cores and commercial NGS companies return unmapped sequence data in FASTQ format.
      The data can be uploaded from users computer, downloaded directly from an ftp server of
      the core facility by providing a URL or from GEO by providing SRA accession number.
    run: ../tools/extract-fastq.cwl
    in:
      compressed_file: fastq_file_R1
      output_prefix:
        default: "merged_R1"
    out: [fastq_file]

  extract_fastq_R2:
    label: "Loading unmapped sequence data for read 2"
    doc: |
      Most DNA cores and commercial NGS companies return unmapped sequence data in FASTQ format.
      The data can be uploaded from users computer, downloaded directly from an ftp server of
      the core facility by providing a URL or from GEO by providing SRA accession number.
    run: ../tools/extract-fastq.cwl
    in:
      compressed_file: fastq_file_R2
      output_prefix:
        default: "merged_R2"
    out: [fastq_file]

  trim_fastq:
    label: "Adapter trimming"
    doc: |
      For libraries sequenced on the Illumina platform it’s recommended to remove adapter sequences
      from the reads. If adapters are not trimmed there is a high risk of reads being unmapped to a
      reference genome. This becomes particularly important when the reads are long and the fragments
      are short - resulting in sequencing adapters at the end of read. If adapter trimming will cause
      all the reads become too short (<30bp), this step will be skipped.
    run: ../tools/trimgalore.cwl
    in:
      input_file: extract_fastq_R1/fastq_file
      input_file_pair: extract_fastq_R2/fastq_file
      dont_gzip:
        default: true
      length:
        default: 30
      trim1:
        default: false
      paired:
        default: true
    out:
      - trimmed_file
      - trimmed_file_pair
      - report_file
      - report_file_pair

  bypass_trim:
    run: ../tools/bypass-trimgalore-pe.cwl
    in:
      original_fastq_file_1: extract_fastq_R1/fastq_file
      trimmed_fastq_file_1: trim_fastq/trimmed_file
      trimming_report_file_1: trim_fastq/report_file
      original_fastq_file_2: extract_fastq_R2/fastq_file
      trimmed_fastq_file_2: trim_fastq/trimmed_file_pair
      trimming_report_file_2: trim_fastq/report_file_pair
      min_reads_count:
        default: 100  # any small number should be good, as we are catching the case when trimgalore discarded all reads
    out:
      - selected_fastq_file_1
      - selected_report_file_1
      - selected_fastq_file_2
      - selected_report_file_2

  rename_upstream:
    run: ../tools/rename.cwl
    in:
      source_file: bypass_trim/selected_fastq_file_1
      target_filename:
        source: extract_fastq_R1/fastq_file
        valueFrom: $(self.basename)
    out:
      - target_file

  rename_downstream:
    run: ../tools/rename.cwl
    in:
      source_file: bypass_trim/selected_fastq_file_2
      target_filename:
        source: extract_fastq_R2/fastq_file
        valueFrom: $(self.basename)
    out:
      - target_file

  fastx_quality_stats_upstream:
    label: "Quality control of unmapped sequence data for read 1"
    doc: |
      Evaluates the quality of your sequence data. Provides per base quality scores as well as
      base frequencies along the reads. These metrics can be used to identify whether your data
      has any problems that should be taken into account in the subsequent analysis steps.
    run: ../tools/fastx-quality-stats.cwl
    in:
      input_file: rename_upstream/target_file
    out: [statistics_file]

  fastx_quality_stats_downstream:
    label: "Quality control of unmapped sequence data for read 2"
    doc: |
      Evaluates the quality of your sequence data. Provides per base quality scores as well as
      base frequencies along the reads. These metrics can be used to identify whether your data
      has any problems that should be taken into account in the subsequent analysis steps.
    run: ../tools/fastx-quality-stats.cwl
    in:
      input_file: rename_downstream/target_file
    out: [statistics_file]

  decontaminate_with_kneaddata:
    label: "Removal of human contamination from read data using kneaddata"
    doc: |
      Tool filters out human contamination from input read files. Outputs both contaminate labeled read pairs and cleaned read pairs.
    run: ../tools/wgs-kneaddata-pe.cwl
    in:
      read1file: rename_upstream/target_file
      read2file: rename_downstream/target_file
    out:
      - kneaddata_cleaned_R1
      - kneaddata_cleaned_R2
      - kneaddata_contaminated_R1
      - kneaddata_contaminated_R2
      - kneaddata_log
      - stdout_log
      - stderr_log

  kraken2_classify:
    label: "Kraken2 taxonomic classification of cleaned sequence reads"
    doc: |
      Assigns taxonomy to each sequence in the input paired end read files, and reports raw
      classificaiton as well as a summary report.
    run: ../tools/k2-classify-pe.cwl
    in:
      k2db: k2db
      read1file: decontaminate_with_kneaddata/kneaddata_cleaned_R1
      read2file: decontaminate_with_kneaddata/kneaddata_cleaned_R2
    out:
      - k2_classified_R1
      - k2_classified_R2
      - k2_unclassified_R1
      - k2_unclassified_R2
      - k2_output
      - k2_report
      - k2_report_tsv
      - k2_stderr
      - krona_html
      - stdout_log
      - stderr_log

  classify_unclassified_k2_reads_with_metaphlan:
    label: "attempt to classify the reads left unclassified by kraken2 with metaphlan latest db"
    doc: |
      Reports abundance tables for input paired end read files. Inputs comes from unclassifed kraken2 output.
    run: ../tools/wgs-metaphlan-pe.cwl
    in:
      read1file: kraken2_classify/k2_unclassified_R1
      read2file: kraken2_classify/k2_unclassified_R2
    out:
      - classification_alignments_bowtie2
      - abundance_profile
      - abundance_profile_scidap
      - abundance_table
      - abundance_table_species
      - stdout_log
      - stderr_log

  classify_cleaned_reads_with_metaphlan:
    label: "classify all cleaned reads from kneaddata with metaphlan latest db"
    doc: |
      Reports abundance tables for input paired end read files. Input comes from kneaddata decontamination output.
    run: ../tools/wgs-metaphlan-pe.cwl
    in:
      read1file: decontaminate_with_kneaddata/kneaddata_cleaned_R1
      read2file: decontaminate_with_kneaddata/kneaddata_cleaned_R2
    out:
      - classification_alignments_bowtie2
      - abundance_profile
      - abundance_profile_scidap
      - abundance_table
      - abundance_table_species
      - stdout_log
      - stderr_log

  functional_assignment_with_humann:
    label: "assign gene families and functions to input data set"
    doc: |
      Reports RPK and CPM for gene families, and renamed human-readable functional assignment for input data set.
    run: ../tools/wgs-humann-pe.cwl
    in:
      read1file: decontaminate_with_kneaddata/kneaddata_cleaned_R1
      read2file: decontaminate_with_kneaddata/kneaddata_cleaned_R2
    out:
      - genefamilies_rpk
      - genefamilies_cpm
      - regroup_to_rxn_cpm
      - regroup_to_rxn_cpm_named
      - stdout_log
      - stderr_log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "WGS Metagenomic pipeline paired-end"
label: "WGS Metagenomic pipeline paired-end"
s:alternateName: "WGS metagenomic pipeline for PE reads using kraken2 kneaddata metaphlan humann"

s:downloadUrl: https://github.com/datirium/workflows/tree/master/workflows/workflows/kraken2-classify-pe.cwl
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
  This workflow taxonomically classifies paired-end sequencing reads in FASTQ format for a SINGLE sample.
  Reads are first adapter trimmed with trimgalore and filtered using kneaddata with a bmtagger database.
  The resulting cleaned reads are classified using Kraken2 and a user-selected pre-built database from a list of
  [genomic index files](https://benlangmead.github.io/aws-indexes/k2). Unaligned reads are then classified using
  metaphlan4 with the mpa_vJan21_CHOCOPhlAnSGB_202103 database. The kraken2 report is used to generate a krona plot
  visualization of the abundance profile. Cleaned reads are also run through HUMANN3 using the uniref90 diamond
  databaseto produce a gene abundance report and metabolic pathway file. The latter is used for abundance coverage
  and functional assignment.

  ### __Inputs__
  Kraken2 database for taxonomic classification:
    - Standard is recommended
  Read 1 file:
    - FASTA/Q input R1 from a paired end library
  Read 2 file:
    - FASTA/Q input R2 from a paired end library
  Number of threads for steps that support multithreading:
     - Number of threads for steps that support multithreading - default set to `4`
  Advanced Inputs Tab (Optional):
     - Number of bases to clip from the 3p end
     - Number of bases to clip from the 5p end

  ### __Outputs__
   - kraken2 report (abundance profile)
   - krona plot (hierarchical visualization of taxonomic classifications)
   - various log files
   - metabolic pathway file
   - functional assignment

  ### __Data Analysis Steps__
  1. QC raw FASTQ files with fastQC and trimmomatic
     - OUTPUT1: trimmed FASTQ files
  2. Filter human reads out of OUTPUT1 with the KneadData tool ()
     - OUTPUT2: filtered FASTQ files
  3. Classify OUTPUT2 with kraken2 using “Standard” database (Refeq archaea, bacteria, viral, plasmid, human, UniVec_Core)
     - OUTPUT3: taxonomic abundance profile
     - OUTPUT4: FASTQ files of unclassified reads
     - VISUALIZATION1: krakenreport to kronaplot
  4. Attempt to classify OUTPUT4 with MetaPhlAn using “latest” database
     - OUTPUT5: taxonomic abundance profile of unclassified kraken2 reads
  5. Classify OUTPUT2 with MetaPhlAn using “latest” database
     - OUTPUT6: final computed taxon abundances (listed one clade per line, tab-separated from the clade's relative abundance in percent)
       - format: https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-Workshop-on-Genomics-2023#13-metaphlan-output-files
       - used in the multi-sample workflow (https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-Workshop-on-Genomics-2023#15-analyzing-multiple-samples)
  6. Use OUTPUT2 in Metagenome functional profiling/assignment with HUMAnN using “uniref : uniref90_diamond” database
     - database link: http://huttenhower.sph.harvard.edu/humann_data/uniprot/uniref_annotated/uniref90_annotated_v201901b_full.tar.gz
     -   OUTPUT7: *_genefamilies.tsv, contains the abundances of each gene family in the community in reads per kilobase (RPK) units
     -   OUTPUT8: *_pathabundance.tsv, lists the abundances of each pathway in the community, also in RPK units as described for gene families
     -   OUTPUT9: normalized_genefamilies-cpm.tsv, contains the normalized abundances of each gene family in counts per million (CPM) units
     -   OUTPUT10: rxn-cpm.tsv, regroup our CPM-normalized gene family abundance values to MetaCyc reaction (RXN) abundances
     -   https://github.com/biobakery/MetaPhlAn/wiki/HUMAnN-Workshop-on-Genomics-2023#3-manipulating-humann-output-tables

  ### __References__
    - McIver LJ, Abu-Ali G, Franzosa EA, Schwager R, Morgan XC, Waldron L, Segata N, Huttenhower C. bioBakery: a meta'omic analysis environment. Bioinformatics. 2018 Apr 1;34(7):1235-1237. PMID: 29194469
    - Benson G. Tandem repeats finder: a program to analyze DNA sequences. Nucleic Acids Res. 1999; 27(2):573–580. doi:10.1093/nar/27.2.573
    - [**Extending and improving metagenomic taxonomic profiling with uncharacterized species using MetaPhlAn 4.**](https://doi.org/10.1038/s41587-023-01688-w) Aitor Blanco-Miguez, Francesco Beghini, Fabio Cumbo, Lauren J. McIver, Kelsey N. Thompson, Moreno Zolfo, Paolo Manghi, Leonard Dubois, Kun D. Huang, Andrew Maltez Thomas, Gianmarco Piccinno, Elisa Piperni, Michal Punčochář, Mireia Valles-Colomer, Adrian Tett, Francesca Giordano, Richard Davies, Jonathan Wolf, Sarah E. Berry, Tim D. Spector, Eric A. Franzosa, Edoardo Pasolli, Francesco Asnicar, Curtis Huttenhower, Nicola Segata. Nature Biotechnology (2023)
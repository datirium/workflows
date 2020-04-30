cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement


'sd:upstream':
  deseq_experiment:
    - "deseq.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 1

  read_counts_file:
    type: File
    format: "http://edamontology.org/format_3709"
    label: "DESeq experiment"
    doc: "Input gene expression dataset file in txt or gct format. Same with GSEA"
    'sd:upstreamSource': "deseq_experiment/read_counts_file"
    'sd:localLabel': true

  phenotypes_file:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "DESeq experiment"
    doc: "Input class vector (phenotype) file in CLS format. Same with GSEA"
    'sd:upstreamSource': "deseq_experiment/phenotypes_file"
    'sd:localLabel': true

  gene_set_database:
    type:
    - type: enum
      name: "genesetdatabase"
      symbols:
      - ARCHS4_Cell-lines
      - ARCHS4_IDG_Coexp
      - ARCHS4_Kinases_Coexp
      - ARCHS4_TFs_Coexp
      - ARCHS4_Tissues
      - Achilles_fitness_decrease
      - Achilles_fitness_increase
      - Aging_Perturbations_from_GEO_down
      - Aging_Perturbations_from_GEO_up
      - Allen_Brain_Atlas_down
      - Allen_Brain_Atlas_up
      - BioCarta_2013
      - BioCarta_2015
      - BioCarta_2016
      - BioPlanet_2019
      - BioPlex_2017
      - CCLE_Proteomics_2020
      - CORUM
      - Cancer_Cell_Line_Encyclopedia
      - ChEA_2013
      - ChEA_2015
      - ChEA_2016
      - Chromosome_Location
      - Chromosome_Location_hg19
      - ClinVar_2019
      - DSigDB
      - Data_Acquisition_Method_Most_Popular_Genes
      - DepMap_WG_CRISPR_Screens_Broad_CellLines_2019
      - DepMap_WG_CRISPR_Screens_Sanger_CellLines_2019
      - DisGeNET
      - Disease_Perturbations_from_GEO_down
      - Disease_Perturbations_from_GEO_up
      - Disease_Signatures_from_GEO_down_2014
      - Disease_Signatures_from_GEO_up_2014
      - DrugMatrix
      - Drug_Perturbations_from_GEO_2014
      - Drug_Perturbations_from_GEO_down
      - Drug_Perturbations_from_GEO_up
      - ENCODE_Histone_Modifications_2013
      - ENCODE_Histone_Modifications_2015
      - ENCODE_TF_ChIP-seq_2014
      - ENCODE_TF_ChIP-seq_2015
      - ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X
      - ESCAPE
      - Elsevier_Pathway_Collection
      - Enrichr_Libraries_Most_Popular_Genes
      - Enrichr_Submissions_TF-Gene_Coocurrence
      - Epigenomics_Roadmap_HM_ChIP-seq
      - GO_Biological_Process_2013
      - GO_Biological_Process_2015
      - GO_Biological_Process_2017
      - GO_Biological_Process_2017b
      - GO_Biological_Process_2018
      - GO_Cellular_Component_2013
      - GO_Cellular_Component_2015
      - GO_Cellular_Component_2017
      - GO_Cellular_Component_2017b
      - GO_Cellular_Component_2018
      - GO_Molecular_Function_2013
      - GO_Molecular_Function_2015
      - GO_Molecular_Function_2017
      - GO_Molecular_Function_2017b
      - GO_Molecular_Function_2018
      - GTEx_Tissue_Sample_Gene_Expression_Profiles_down
      - GTEx_Tissue_Sample_Gene_Expression_Profiles_up
      - GWAS_Catalog_2019
      - GeneSigDB
      - Gene_Perturbations_from_GEO_down
      - Gene_Perturbations_from_GEO_up
      - Genes_Associated_with_NIH_Grants
      - Genome_Browser_PWMs
      - HMDB_Metabolites
      - HMS_LINCS_KinomeScan
      - HomoloGene
      - HumanCyc_2015
      - HumanCyc_2016
      - Human_Gene_Atlas
      - Human_Phenotype_Ontology
      - InterPro_Domains_2019
      - Jensen_COMPARTMENTS
      - Jensen_DISEASES
      - Jensen_TISSUES
      - KEA_2013
      - KEA_2015
      - KEGG_2013
      - KEGG_2015
      - KEGG_2016
      - KEGG_2019_Human
      - KEGG_2019_Mouse
      - Kinase_Perturbations_from_GEO_down
      - Kinase_Perturbations_from_GEO_up
      - L1000_Kinase_and_GPCR_Perturbations_down
      - L1000_Kinase_and_GPCR_Perturbations_up
      - LINCS_L1000_Chem_Pert_down
      - LINCS_L1000_Chem_Pert_up
      - LINCS_L1000_Ligand_Perturbations_down
      - LINCS_L1000_Ligand_Perturbations_up
      - Ligand_Perturbations_from_GEO_down
      - Ligand_Perturbations_from_GEO_up
      - MCF7_Perturbations_from_GEO_down
      - MCF7_Perturbations_from_GEO_up
      - MGI_Mammalian_Phenotype_2013
      - MGI_Mammalian_Phenotype_2017
      - MGI_Mammalian_Phenotype_Level_3
      - MGI_Mammalian_Phenotype_Level_4
      - MGI_Mammalian_Phenotype_Level_4_2019
      - MSigDB_Computational
      - MSigDB_Oncogenic_Signatures
      - Microbe_Perturbations_from_GEO_down
      - Microbe_Perturbations_from_GEO_up
      - Mouse_Gene_Atlas
      - NCI-60_Cancer_Cell_Lines
      - NCI-Nature_2015
      - NCI-Nature_2016
      - NIH_Funded_PIs_2017_AutoRIF_ARCHS4_Predictions
      - NIH_Funded_PIs_2017_GeneRIF_ARCHS4_Predictions
      - NIH_Funded_PIs_2017_Human_AutoRIF
      - NIH_Funded_PIs_2017_Human_GeneRIF
      - NURSA_Human_Endogenous_Complexome
      - OMIM_Disease
      - OMIM_Expanded
      - Old_CMAP_down
      - Old_CMAP_up
      - PPI_Hub_Proteins
      - Panther_2015
      - Panther_2016
      - Pfam_Domains_2019
      - Pfam_InterPro_Domains
      - PheWeb_2019
      - Phosphatase_Substrates_from_DEPOD
      - ProteomicsDB_2020
      - RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO
      - Rare_Diseases_AutoRIF_ARCHS4_Predictions
      - Rare_Diseases_AutoRIF_Gene_Lists
      - Rare_Diseases_GeneRIF_ARCHS4_Predictions
      - Rare_Diseases_GeneRIF_Gene_Lists
      - Reactome_2013
      - Reactome_2015
      - Reactome_2016
      - SILAC_Phosphoproteomics
      - SubCell_BarCode
      - SysMyo_Muscle_Gene_Sets
      - TF-LOF_Expression_from_GEO
      - TF_Perturbations_Followed_by_Expression
      - TRANSFAC_and_JASPAR_PWMs
      - TRRUST_Transcription_Factors_2019
      - Table_Mining_of_CRISPR_Studies
      - TargetScan_microRNA
      - TargetScan_microRNA_2017
      - Tissue_Protein_Expression_from_Human_Proteome_Map
      - Tissue_Protein_Expression_from_ProteomicsDB
      - Transcription_Factor_PPIs
      - UK_Biobank_GWAS_v1
      - Virus-Host_PPI_P-HIPSTer_2020
      - VirusMINT
      - Virus_Perturbations_from_GEO_down
      - Virus_Perturbations_from_GEO_up
      - WikiPathways_2013
      - WikiPathways_2015
      - WikiPathways_2016
      - WikiPathways_2019_Human
      - WikiPathways_2019_Mouse
      - dbGaP
      - huMAP
      - lncHUB_lncRNA_Co-Expression
      - miRTarBase_2017
    default: "GeneSigDB"
    label: "Gene set database"
    doc: "Gene set database"

  permutation_type:
    type:
    - "null"
    - type: enum
      name: "permutationtype"
      symbols:
      - gene_set
      - phenotype
    default: "gene_set"
    label: "Permutation type"
    doc: "Permutation type. Default: gene_set"

  permutation_count:
    type: int?
    default: 1000
    label: "Number of random permutations"
    doc: "Number of random permutations. For calculating esnulls. Default: 1000"

  min_gene_set_size:
    type: int?
    default: 15
    label: "Min size of input genes presented in Gene Sets"
    doc: "Min size of input genes presented in Gene Sets. Default: 15"
    'sd:layout':
      advanced: true

  max_gene_set_size:
    type: int?
    default: 500
    label: "Max size of input genes presented in Gene Sets"
    doc: "Max size of input genes presented in Gene Sets. Default: 500"
    'sd:layout':
      advanced: true

  ranking_metrics:
    type:
    - "null"
    - type: enum
      name: "rankingmetrics"
      symbols:
      - signal_to_noise
      - t_test
      - ratio_of_classes
      - diff_of_classes
      - log2_ratio_of_classes
    default: "signal_to_noise"
    label: "Methods to calculate correlations of ranking metrics"
    doc: "Methods to calculate correlations of ranking metrics. Default: log2_ratio_of_classes"

  ascending_rank_sorting:
    type: boolean?
    default: false
    label: "Ascending rank metric sorting order"
    doc: "Ascending rank metric sorting order. Default: False"

  graphs_count:
    type: int?
    default: 20
    label: "Numbers of top graphs produced"
    doc: "Numbers of top graphs produced. Default: 20"
    'sd:layout':
      advanced: true

  seed:
    type: int?
    label: "Number of random seed. Default: None"
    doc: "Number of random seed. Default: None"
    'sd:layout':
      advanced: true

  threads:
    type: int?
    default: 4
    label: "Number of threads"
    doc: "Number of threads for those steps that support multithreading"
    'sd:layout':
      advanced: true


outputs:

  gseapy_enrichment_report:
    type: File?
    format: "http://edamontology.org/format_3475"
    label: "Enrichment report"
    doc: "Enrichment report"
    outputSource: convert_to_tsv/output_file
    "sd:visualPlugins":
    - syncfusiongrid:
        tab: "Gene Set Enrichment"
        Title: "Gene Set Enrichment"

  gseapy_enrichment_plots:
    type: File[]
    format: "http://edamontology.org/format_3508"
    label: "Enrichment plots"
    doc: "Enrichment plots"
    outputSource: run_gseapy/enrichment_plots
    "sd:visualPlugins":
    - image:
        tab: "Plots"
        Caption: "Enrichment plot"

  gseapy_enrichment_heatmaps:
    type: File[]
    format: "http://edamontology.org/format_3508"
    label: "Enrichment heatmaps"
    doc: "Enrichment heatmaps"
    outputSource: run_gseapy/enrichment_heatmaps
    "sd:visualPlugins":
    - image:
        tab: "Plots"
        Caption: "Enrichment heatmap"

  gseapy_stdout_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "GSEApy stdout log"
    doc: "GSEApy stdout log"
    outputSource: run_gseapy/stdout_log

  gseapy_stderr_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "GSEApy stderr log"
    doc: "GSEApy stderr log"
    outputSource: run_gseapy/stderr_log


steps:

  run_gseapy:
    run: ../tools/gseapy.cwl
    in:
      read_counts_file: read_counts_file
      phenotypes_file: phenotypes_file
      gene_set_database: gene_set_database
      permutation_type: permutation_type
      permutation_count: permutation_count
      min_gene_set_size: min_gene_set_size
      max_gene_set_size: max_gene_set_size
      ranking_metrics: ranking_metrics
      ascending_rank_sorting: ascending_rank_sorting
      graphs_count: graphs_count
      seed: seed
      threads: threads
    out:
      - enrichment_report
      - enrichment_plots
      - enrichment_heatmaps
      - stdout_log
      - stderr_log

  convert_to_tsv:
    run: ../tools/custom-bash.cwl
    in:
      input_file: run_gseapy/enrichment_report
      script:
        default: |
          cat "$0" | tr "," "\t" > `basename $0 csv`tsv
    out: [output_file]


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "GSEApy - Gene Set Enrichment Analysis in Python"
label: "GSEApy - Gene Set Enrichment Analysis in Python"
s:alternateName: "GSEApy - Gene Set Enrichment Analysis in Python"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/gseapy.cwl
s:codeRepository: https://github.com/datirium/workflows
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


# doc:
#   $include: ../descriptions/gseapy.md


doc: |
  GSEAPY: Gene Set Enrichment Analysis in Python
  ==============================================

  Gene Set Enrichment Analysis is a computational method that determines whether an a priori
  defined set of genes shows statistically significant, concordant differences between two
  biological states (e.g. phenotypes).

  GSEA requires as input an expression dataset, which contains expression profiles for multiple samples.
  While the software supports multiple input file formats for these datasets, the tab-delimited GCT format
  is the most common. The first column of the GCT file contains feature identifiers (gene ids or symbols in
  the case of data derived from RNA-Seq experiments). The second column contains a description of the feature;
  this column is ignored by GSEA and may be filled with “NA”s. Subsequent columns contain the expression
  values for each feature, with one sample's expression value per column. It is important to note that there
  are no hard and fast rules regarding how a GCT file's expression values are derived. The important point is
  that they are comparable to one another across features within a sample and comparable to one another
  across samples. Tools such as DESeq2 can be made to produce properly normalized data (normalized counts)
  which are compatible with GSEA.
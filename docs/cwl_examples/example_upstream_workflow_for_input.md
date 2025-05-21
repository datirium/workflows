# Example of workflow using upstream workflows for inputs


> include ```sd:upstream```'s at the root of the cwl file (above inputs)
```yaml
...
'sd:upstream':
  rnaseq_sample_untreated:
    - "rnaseq-se.cwl"
    - "rnaseq-pe.cwl"
    - "rnaseq-se-dutp.cwl"
    - "rnaseq-pe-dutp.cwl"
    - "rnaseq-se-dutp-mitochondrial.cwl"
    - "rnaseq-pe-dutp-mitochondrial.cwl"
    - "trim-rnaseq-pe.cwl"
    - "trim-rnaseq-se.cwl"
    - "trim-rnaseq-pe-dutp.cwl"
    - "trim-rnaseq-se-dutp.cwl"
  rnaseq_sample_treated:
    - "rnaseq-se.cwl"
    - "rnaseq-pe.cwl"
    - "rnaseq-se-dutp.cwl"
    - "rnaseq-pe-dutp.cwl"
    - "rnaseq-se-dutp-mitochondrial.cwl"
    - "rnaseq-pe-dutp-mitochondrial.cwl"
    - "trim-rnaseq-pe.cwl"
    - "trim-rnaseq-se.cwl"
    - "trim-rnaseq-pe-dutp.cwl"
    - "trim-rnaseq-se-dutp.cwl"

inputs: 
...
```

> in order to use upstream outputs as inputs, include a reference to the upstream workflows output files using the ```sd:upstreamSource``` keyword
```yaml
... 
inputs:

  untreated_files:
    type: File[]
    format:
     - "http://edamontology.org/format_3752"
     - "http://edamontology.org/format_3475"
    label: "Untreated input CSV/TSV files"
    doc: "Untreated input CSV/TSV files"
    'sd:upstreamSource': "rnaseq_sample_untreated/rpkm_common_tss"
    'sd:localLabel': true
...
```

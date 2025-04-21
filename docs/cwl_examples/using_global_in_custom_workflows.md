# Example of how to use Global workflows in Custom workflows

When working with custom workflows, there is often a need to allow these custom workflows to use  global (datirium) workflows as their upstream.

This can be done by adding an html reference to the list of upstream workflows usable.

EX:

```yaml

'sd:upstream':
  genome_indices:
    - "genome-indices.cwl"
    # NOTE: not actual link to workflow, but follows URL format
    # actual link would include "/tree/" or "/head/" somewhere
    - "https://github.com/datirium/workflows/workflows/genome-indoces.cwl"

```
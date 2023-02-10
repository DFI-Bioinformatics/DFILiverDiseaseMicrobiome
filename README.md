# Data Repo for Prebiotic activity of lactulose optimizes gut metabolites and prevents systemic infection in liver disease patients

## Folder structure

1. **code**: contains all R scripts for generating heatmap, volcano plots and/or taxonomy barplots, etc. Scripts are organized by number. Resultant plots are labeled with same number as the script. 
2. **results**: contains all unmodified graphs in PDF format. Graphs published along with the paper have been slightly rearranged and beautified in Adobe Illustrator.
3. **data**: contains raw or derived data or template files
- `LD850.meta.quant.metabolomics.csv`: the most important file. It contains clinical variables, such as stool consistency, lactulose, Bifidobacterium group (10% is the cutoff), SBP, bacteremia status, or metagenomically derived toxin RPKM values, quant metabolomics readings from SCFA (unit: mM) and bile acid (unit: ug/ml) panels, etc.
- `LD850.metaphlan.csv`: taxonomy annotation from [metaphlan4](https://github.com/biobakery/MetaPhlAn)
- `LD850.qual.metabolomics.csv`: a long format of qualitative metabolomics values. There is not unit associated with them.
- `LD850.ratios.quant.metabolomics.csv`: primary to secondary bile acid ratios calculated from mM.
- `embedding.LD850.csv`: [taxUmap](https://github.com/jsevo/taxumap) coordinates
- `In vitro and mouse data.xlsx`: data for *in vitro* and mouse experiment. They are organized by tabs.
- `qual_compounds.csv`: a list of all qualitative compounds to use in the heatmap
- `qual_heatmap_lookup.csv`: a template for all qualitative compounds and their classes


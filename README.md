# BeginNGS Public Compendium
Reproducible research compedium to accompany "Cost-neutral newborn population screening for 412 genetic diseases by genome sequencing with large diplotype models"
![workflow_diagram](https://github.com/TileDB-Inc/BeginNGSPub/assets/147991/51bb5965-0e31-4bba-889b-b30850098538)

The workflow diagram provides a high-level depiction of the steps involved in this analysis [workflow_diagram.pdf](https://github.com/TileDB-Inc/BeginNGSPub/files/15042044/workflow_diagram.pdf)

# Source code
Source code to run the analysis is included here as a [Jupyter notebook](beginNGS_compedium_notebook.ipynb). A [Python script](beginNGS_compedium.py) conversion is also provided.

# Data
- [Variants of interest](data/variants_of_interest_20231108.csv) (currently BeginNGS v2, 53,855 P and LP variants that map to 342 genes, 412 SCGD, and 1,603 SCGD therapeutic interventions) is normally encapsulated as a fixed resource inthe UDF, but can be implemented as an parameter. This is pre-annotated with consequence and population frequency information, but only chr-pos-ref-alt is used for the query itself
- [Blocklist](data/blocklist_20240329.csv) - entries classified as NSDCC (non-severe disease causing in childhood)
- [MOI](data/moi_20240805.txt) - mode of inheritance information

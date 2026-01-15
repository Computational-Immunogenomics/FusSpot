# Systematic identification of fusion derived neoepitopes across cancer types.

---

# Overview

The underlying study data comes from the Hartwig Medical Database, and data access was obtained through an approved public data request (DR347). To request access to the data see instructions described on the Hartwig Medical Foundation website (https://www.hartwigmedicalfoundation.nl/en/data/data-access-request/).   

To reproduce the our analysis, you need data request access to both the LINX and NEO output from WiGiTs tools, as well as the clinical data. 

---

# Folders and files

- [0_prepare_data](0_prepare_data/)
    - Collect and organize data from LINX fusions and high likelihood NEO neoepitope predictions.
    - Clean clinical data to define cohort definitions base on primary tumor type and histological markers.
    - Join cohorts, fusions, and neoepitope output together based on sampleIds.
    - Compute aggregations to identify peptide and fusion counts across cohorts.
- [1_figures](1_figures/)
    - Code to reproduce figure 2 and supplementary figures in the manuscript. 
- map.r,shortcuts.r,settings.r
    - Files to specify file locations, helper functions, and plot settings. 

![Diagram](schematic.png)

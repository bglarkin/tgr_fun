# tgr-fun
Soil fungi response to tallgrass prairie restoration in Illinois and Wisconsin

## Contents
- [Sequence Data](sequence_data.md) describes procedures to extract, transform, and load the raw OTU data. 
   - [Raw Sequence Data](otu_tables) are housed in a separate directory.
   - After ETL, [Clean Data](clean_data) are housed in a separate directory. 
- [Fungal Ecology](fungal_ecology.md) comprises analyses of alpha and beta diversity,
correlations with plant abundance, and constrained analysis of explanatory variables.
- [Edaphic](edaphic.md) describes soil chemical properties and displays differences among sites.
- Custom [Functions](functions.md) are housed separately to save space in other scripts. 
- Raw [Code](code) is housed in a separate directory to reduce root directory clutter.
- The [Resources](resources) folder contains the style sheet for figures, sidecar figure
files used by GitHub, and the [Render Reports](resources/render_reports.R) script that handles
report rendering and file path housekeeping to keep the root directory clutter-free.

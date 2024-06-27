# Computer vision reveals climate predictors of global functional trait variation in bees

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12572899.svg)](https://doi.org/10.5281/zenodo.12572899)

# Abstract
Many studies have linked climate to functional traits at local scales, yet evidence that climate variables predict global, macroecological trends in trait variation remains limited. In particular, data is limiting for functional traits that are prohibitively challenging to quantify at scale, yet many of these complex traits (e.g., coloration, pilosity) have outsize effects on organismal thermal performance. To overcome this challenge, we leveraged techniques in deep learning and trained a computer vision model to perform image segmentation tasks that allow us to quantify pilosity and coloration from images of bee specimens. We implemented this computer vision model on a large and diverse image dataset representing over 600 bee species from all bee families and subfamilies. We demonstrate that climate shapes variation in these traits at a global scale, with bee lightness increasing with maximum environmental temperatures (Thermal Melanism Hypothesis) and decreasing with annual precipitation (Gloger’s Rule). Pilosity (measured as hair coverage) was strongly correlated with lightness, with light-colored hairs covering darker integument and lightening bees’ overall surface coloration. Consequently, hair coverage increased with maximum temperatures, likely as an adaptation to reduce risk of overheating in hot, dry climates. Bees in deserts and xeric regions were significantly lighter and hairier than bees from other biomes. Together, these results provide support for major ecogeographical rules in functional trait variation, and emphasize the role of climate in shaping bee phenotypic diversity. 

# Repository Directory
### Code: Contains code for data analysis in R
pilosity_CV_analysis_May2024.R: In this document, we analyze morphometric data (ITD, head width, costal vein length, dry mass, volume, surface area) to understand inter- and intra-specific predictors of body size.

## Data: Contains the raw specimen and trait data
pilosity_lightness_data_june2024.csv: In this document, we present the raw specimen and trait data used in the analysis (pilosity_CV_analysis_May2024.R).

pilosity_lightness_data_trait_format.csv: In this document, we present the same raw data as above, but in a standardized format that lends itself to functional trait data sharing (see this [functional trait data sharing repository](https://github.com/mostwald/Functional-trait-review) for descriptors of column headers). (in progress)

specimen_list_image_figure.csv: In this document, we present specimen data for images displayed in figure X. Specimens are ordered (position column) following the image grid, reading from left to right and top to bottom.



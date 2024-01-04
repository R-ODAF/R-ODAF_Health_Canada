# Contrasts file

This folder should contain the following file:

contrasts.txt

The contrasts file is a tab-delimited text file that describe which contrasts of interest should be tested with the results() function of DESeq2.

Contrasts will provide the information necessary to make comparisons between samples. The first column must be an experimental grouping of interest (e.g., exposed) and the second column must be the baseline group against which the experimental group should be compared (e.g, vehicle_control). The names must correspond to entries in the metadata table. For example, if you have a column "dose" in metadata, then we would expect that the contrasts file contains one row for each dose group (e.g., 1000, 100, 10), while the second column might be 0 as the control for all those groups. The contents of the contrasts.txt file for this example would be:

10  0
100 0
1000    0
# Contrasts file

This folder should contain at least one tab-delimited text file that describe(s) which contrasts of interest should be tested with the results() function of DESeq2 when running the differential expression module of the workflow. You may wish to have different contrast files in order to do multiple rounds of differential expression analysis.

Specify the name of the contrasts file (ex. "contrasts.txt" in the config file).

Contrasts will provide the information necessary to make comparisons between samples. The first column must be an experimental grouping of interest (e.g., exposed) and the second column must be the baseline group against which the experimental group should be compared (e.g, vehicle_control). The names must correspond to entries in one of the columns in the metadata table; that column must be specified in the config file's "design" parameter. 

For example, if you have a column "dose" in metadata, then we would expect that the contrasts file contains one row for each dose group (e.g., 1000, 100, 10), while the second column might be 0 as the control for all those groups. 

In this example, the relevant sections of the metadata file would be:

| sample_ID | chemical | dose | group        |
|-----------|----------|------|--------------|
| s1        | BaP      | 1    | BaP_1        |
| s2        | BaP      | 1    | BaP_1        |
| s3        | BaP      | 1    | BaP_1        |
| s4        | BaP      | 10   | BaP_10       |
| s5        | BaP      | 10   | BaP_10       |
| s6        | BaP      | 10   | BaP_10       |
| s7        | DMSO     | 0    | vehicle_ctrl |
| s8        | DMSO     | 0    | vehicle_ctrl |
| s9        | DMSO     | 0    | vehicle_ctrl |

The relevant parameter in the config file would be:
design: "group"

And contents of the contrasts file would be:

BaP_1	vehicle_ctrl \\
BaP_10	vehicle_ctrl

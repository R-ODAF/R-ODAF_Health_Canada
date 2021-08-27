# Metadata

This folder should contain at least the following two files:  

1. metadata.txt  
2. contrasts.txt  

Metadata is information about your samples. One column should be a unique identifier, and should be labelled sample_id. Any additional columns should be descriptive information about the samples.

Contrasts will provide the information necessary to make comparisons between samples. The first column must be an experimental grouping of interest (e.g., exposed) and the second column must be the baseline group against which the experimental group should be compared (e.g, vehicle_control). The names must correspond to entries in the metadata table. For example, if you have a column "dose" in metadata, then we would expect that the contrasts file contains one row for each dose group (e.g., 1000, 100, 10), while the second column might be 0 as the control for all those groups.

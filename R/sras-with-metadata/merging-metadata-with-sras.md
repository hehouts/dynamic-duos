---
title: "merging metadata with sras"
output: 
  html_document: 
    keep_md: yes
date: '2022-03-02'
---

## Load libs & data

```r
library(tidyverse)
```

```
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──
```

```
## ✓ ggplot2 3.3.5     ✓ purrr   0.3.4
## ✓ tibble  3.1.6     ✓ dplyr   1.0.8
## ✓ tidyr   1.2.0     ✓ stringr 1.4.0
## ✓ readr   2.1.2     ✓ forcats 0.5.1
```

```
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## x dplyr::filter() masks stats::filter()
## x dplyr::lag()    masks stats::lag()
```

```r
library(janitor)
```

```
## 
## Attaching package: 'janitor'
```

```
## The following objects are masked from 'package:stats':
## 
##     chisq.test, fisher.test
```

```r
# "Manifest" -- call this "file manifest"
#   # has fileID, md5sum, file size, download urls,  and sample ID
file_mani <- readr::read_tsv("metadata-of-2342-samples/hmp_manifest_47683f7597.tsv")
```

```
## Rows: 2341 Columns: 5
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (4): file_id, md5, urls, sample_id
## dbl (1): size
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
# this is less likely to be useful?

#"Manifest Metadata" -- call this "sampleID manifest"
#   # has: sample_id	subject_id	subject_uuid	sample_body_site	visit_number	subject_gender	subject_race	study_full_name	project_name
samp_id_mani <- readr::read_tsv("metadata-of-2342-samples/hmp_manifest_metadata_b7ce70666.tsv")
```

```
## Rows: 2341 Columns: 9
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (7): sample_id, subject_uuid, sample_body_site, subject_gender, subject_...
## dbl (2): subject_id, visit_number
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```
## Conclusions
1. `sample_id_mani` has duplicated rows, so to join to make `full_file_manifest`, I needed to pull unique row. 
2. to join these with the ihmp metadata, the `ihmp_id`s will have to be pulled from the URL in `file_mani`

## Code Summary

```r
file_mani <- readr::read_tsv("metadata-of-2342-samples/hmp_manifest_47683f7597.tsv")
```

```
## Rows: 2341 Columns: 5
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (4): file_id, md5, urls, sample_id
## dbl (1): size
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
samp_id_mani <- readr::read_tsv("metadata-of-2342-samples/hmp_manifest_metadata_b7ce70666.tsv")
```

```
## Rows: 2341 Columns: 9
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (7): sample_id, subject_uuid, sample_body_site, subject_gender, subject_...
## dbl (2): subject_id, visit_number
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
full_file_manifest <- left_join(x = file_mani, y = unique(samp_id_mani))%>%
  separate(urls, into = c(NA, "ihmp_sample_id", NA), sep = "/raw/|.tar|.fastq.gz", remove = F) %>% 
  relocate(ihmp_sample_id, sample_id, file_id) %>% 
  arrange(ihmp_sample_id, sample_id, file_id)
```

```
## Joining, by = "sample_id"
```

```r
full_file_manifest
```

```
## # A tibble: 2,341 × 14
##    ihmp_sample_id sample_id   file_id md5     size urls  subject_id subject_uuid
##    <chr>          <chr>       <chr>   <chr>  <dbl> <chr>      <dbl> <chr>       
##  1 CSM5FZ3N_P     cad5930efe… 599427… 0081… 1.58e9 http…       3001 1419f08f554…
##  2 CSM5FZ3R_P     5994274eae… 599427… 93ee… 2.61e9 http…       3001 1419f08f554…
##  3 CSM5FZ3T_P     cad5930efe… 599427… 1331… 1.72e9 http…       3002 1419f08f554…
##  4 CSM5FZ3V_P     5994274eae… 599427… 291b… 2.85e9 http…       3001 1419f08f554…
##  5 CSM5FZ3X_P     5994274eae… 599427… 8bf2… 2.24e9 http…       3002 1419f08f554…
##  6 CSM5FZ3Z_P     5994274eae… 599427… c146… 3.82e9 http…       3002 1419f08f554…
##  7 CSM5FZ42_P     5994274eae… 599427… 6659… 5.11e9 http…       3002 1419f08f554…
##  8 CSM5FZ44_P     cad5930efe… 599427… bdde… 3.01e9 http…       3002 1419f08f554…
##  9 CSM5FZ46_P     5994274eae… 599427… 4632… 2.45e9 http…       3002 1419f08f554…
## 10 CSM5FZ4A_P     cad5930efe… 599427… 9044… 2.96e9 http…       3004 1419f08f554…
## # … with 2,331 more rows, and 6 more variables: sample_body_site <chr>,
## #   visit_number <dbl>, subject_gender <chr>, subject_race <chr>,
## #   study_full_name <chr>, project_name <chr>
```



## Methods/Notes

### why this merge is complex
shows tables arent merging perfectly on `sample_id`.


```r
nrow(file_mani)
```

```
## [1] 2341
```

```r
nrow(samp_id_mani)
```

```
## [1] 2341
```


```r
#inner join should be <= nrow, (2341), but wont be if there are duplicate rows!?
nrow(inner_join(x = file_mani, y = samp_id_mani, by = "sample_id"))
```

```
## [1] 3873
```

### properly merging file manifest & sample_id (metadata) manifest

Shows we have redundant sample IDs, no redundant file_ids in `file_mani` df.

```r
# distinct sample ids
n_distinct(file_mani$sample_id)
```

```
## [1] 1613
```

```r
n_distinct(samp_id_mani$sample_id)
```

```
## [1] 1613
```

```r
# distinct file IDs
n_distinct(file_mani$file_id)
```

```
## [1] 2341
```

Shows that sample_ids each have 2 files, from 1 dl link, that probably represent r1 & r2. 
also, metadata to match the sample ids in the ihmp metadata are probably best pulled from the url. 

```r
file_mani %>% 
  filter(sample_id == "1419f08f554e0c93f3b62fe90cf540d0")
```

```
## # A tibble: 2 × 5
##   file_id                          md5                      size urls  sample_id
##   <chr>                            <chr>                   <dbl> <chr> <chr>    
## 1 56bda9b020293b4b1d65e1eb25368a8b d8f8c1cad1fa1a62b5df2… 8.55e7 http… 1419f08f…
## 2 1419f08f554e0c93f3b62fe90cf563c5 9fa9cd301d20271365d23… 7.92e8 http… 1419f08f…
```

Shows that there are fully redundant rows in the samp_id_mani, that I suppose stems from the file_mani.

```r
samp_id_mani %>% 
  filter(sample_id == "1419f08f554e0c93f3b62fe90cf540d0")
```

```
## # A tibble: 2 × 9
##   sample_id subject_id subject_uuid sample_body_site visit_number subject_gender
##   <chr>          <dbl> <chr>        <chr>                   <dbl> <chr>         
## 1 1419f08f…       2008 1419f08f554… feces                      25 female        
## 2 1419f08f…       2008 1419f08f554… feces                      25 female        
## # … with 3 more variables: subject_race <chr>, study_full_name <chr>,
## #   project_name <chr>
```

```r
samp_id_mani %>% 
  unique()
```

```
## # A tibble: 1,613 × 9
##    sample_id               subject_id subject_uuid sample_body_site visit_number
##    <chr>                        <dbl> <chr>        <chr>                   <dbl>
##  1 1419f08f554e0c93f3b62f…       2008 1419f08f554… feces                      23
##  2 d39c1941c8f6e8b0f6ead5…       2008 1419f08f554… feces                      14
##  3 1419f08f554e0c93f3b62f…       2008 1419f08f554… feces                       6
##  4 d39c1941c8f6e8b0f6ead5…       2008 1419f08f554… feces                      15
##  5 d39c1941c8f6e8b0f6ead5…       2008 1419f08f554… feces                      21
##  6 d39c1941c8f6e8b0f6ead5…       2008 1419f08f554… feces                      12
##  7 d39c1941c8f6e8b0f6ead5…       2008 1419f08f554… feces                      18
##  8 1419f08f554e0c93f3b62f…       2008 1419f08f554… feces                      28
##  9 1419f08f554e0c93f3b62f…       2008 1419f08f554… feces                      11
## 10 1419f08f554e0c93f3b62f…       2008 1419f08f554… feces                      22
## # … with 1,603 more rows, and 4 more variables: subject_gender <chr>,
## #   subject_race <chr>, study_full_name <chr>, project_name <chr>
```

join properly, yielding the expected 2341 rows

```r
full_file_manifest <- left_join(x = file_mani, y = unique(samp_id_mani))
```

```
## Joining, by = "sample_id"
```

```r
names(full_file_manifest)
```

```
##  [1] "file_id"          "md5"              "size"             "urls"            
##  [5] "sample_id"        "subject_id"       "subject_uuid"     "sample_body_site"
##  [9] "visit_number"     "subject_gender"   "subject_race"     "study_full_name" 
## [13] "project_name"
```

```r
full_file_manifest %>%
  relocate(sample_id) %>%
  arrange(sample_id, file_id)
```

```
## # A tibble: 2,341 × 13
##    sample_id file_id md5     size urls  subject_id subject_uuid sample_body_site
##    <chr>     <chr>   <chr>  <dbl> <chr>      <dbl> <chr>        <chr>           
##  1 1419f08f… 1419f0… 6c83… 1.05e9 http…       3001 1419f08f554… feces           
##  2 1419f08f… 56bda9… 1eb2… 3.60e6 http…       3001 1419f08f554… feces           
##  3 1419f08f… 1419f0… 5334… 1.32e9 http…       3001 1419f08f554… feces           
##  4 1419f08f… 1419f0… e8cb… 9.81e8 http…       3001 1419f08f554… feces           
##  5 1419f08f… 1419f0… b5d6… 1.03e9 http…       3001 1419f08f554… feces           
##  6 1419f08f… 56bda9… ce38… 8.50e7 http…       3001 1419f08f554… feces           
##  7 1419f08f… 1419f0… f652… 2.21e8 http…       3001 1419f08f554… feces           
##  8 1419f08f… 1419f0… 2df8… 1.56e9 http…       3001 1419f08f554… feces           
##  9 1419f08f… 1419f0… 3f87… 1.49e9 http…       3001 1419f08f554… feces           
## 10 1419f08f… 56bda9… 194c… 1.34e8 http…       3001 1419f08f554… feces           
## # … with 2,331 more rows, and 5 more variables: visit_number <dbl>,
## #   subject_gender <chr>, subject_race <chr>, study_full_name <chr>,
## #   project_name <chr>
```

```r
full_file_manifest <- left_join(x = file_mani, y = unique(samp_id_mani)) %>%
  relocate(sample_id) %>%
  arrange(sample_id, file_id)
```

```
## Joining, by = "sample_id"
```

```r
full_file_manifest
```

```
## # A tibble: 2,341 × 13
##    sample_id file_id md5     size urls  subject_id subject_uuid sample_body_site
##    <chr>     <chr>   <chr>  <dbl> <chr>      <dbl> <chr>        <chr>           
##  1 1419f08f… 1419f0… 6c83… 1.05e9 http…       3001 1419f08f554… feces           
##  2 1419f08f… 56bda9… 1eb2… 3.60e6 http…       3001 1419f08f554… feces           
##  3 1419f08f… 1419f0… 5334… 1.32e9 http…       3001 1419f08f554… feces           
##  4 1419f08f… 1419f0… e8cb… 9.81e8 http…       3001 1419f08f554… feces           
##  5 1419f08f… 1419f0… b5d6… 1.03e9 http…       3001 1419f08f554… feces           
##  6 1419f08f… 56bda9… ce38… 8.50e7 http…       3001 1419f08f554… feces           
##  7 1419f08f… 1419f0… f652… 2.21e8 http…       3001 1419f08f554… feces           
##  8 1419f08f… 1419f0… 2df8… 1.56e9 http…       3001 1419f08f554… feces           
##  9 1419f08f… 1419f0… 3f87… 1.49e9 http…       3001 1419f08f554… feces           
## 10 1419f08f… 56bda9… 194c… 1.34e8 http…       3001 1419f08f554… feces           
## # … with 2,331 more rows, and 5 more variables: visit_number <dbl>,
## #   subject_gender <chr>, subject_race <chr>, study_full_name <chr>,
## #   project_name <chr>
```


Working with only the `file_mani`:


```r
full_file_manifest <- left_join(x = file_mani, y = unique(samp_id_mani)) %>%
  arrange(sample_id, file_id) %>% 
  relocate(sample_id)
```

```
## Joining, by = "sample_id"
```

```r
full_file_manifest
```

```
## # A tibble: 2,341 × 13
##    sample_id file_id md5     size urls  subject_id subject_uuid sample_body_site
##    <chr>     <chr>   <chr>  <dbl> <chr>      <dbl> <chr>        <chr>           
##  1 1419f08f… 1419f0… 6c83… 1.05e9 http…       3001 1419f08f554… feces           
##  2 1419f08f… 56bda9… 1eb2… 3.60e6 http…       3001 1419f08f554… feces           
##  3 1419f08f… 1419f0… 5334… 1.32e9 http…       3001 1419f08f554… feces           
##  4 1419f08f… 1419f0… e8cb… 9.81e8 http…       3001 1419f08f554… feces           
##  5 1419f08f… 1419f0… b5d6… 1.03e9 http…       3001 1419f08f554… feces           
##  6 1419f08f… 56bda9… ce38… 8.50e7 http…       3001 1419f08f554… feces           
##  7 1419f08f… 1419f0… f652… 2.21e8 http…       3001 1419f08f554… feces           
##  8 1419f08f… 1419f0… 2df8… 1.56e9 http…       3001 1419f08f554… feces           
##  9 1419f08f… 1419f0… 3f87… 1.49e9 http…       3001 1419f08f554… feces           
## 10 1419f08f… 56bda9… 194c… 1.34e8 http…       3001 1419f08f554… feces           
## # … with 2,331 more rows, and 5 more variables: visit_number <dbl>,
## #   subject_gender <chr>, subject_race <chr>, study_full_name <chr>,
## #   project_name <chr>
```


```r
#deets on the regex-based sep function: https://stackoverflow.com/questions/43130645/multiple-separate-arguments-in-tidyrs-separate-function


# this splits on the prefix, and one or two possible suffixes. visually inspected results, seems to work. see filtered evidence below
full_file_manifest %>%
  separate(urls, into = c(NA, "ihmp_sample_id", NA), sep = "/raw/|.tar|.fastq.gz", remove = F) %>% 
  relocate(ihmp_sample_id) %>% 
  select(ihmp_sample_id, urls) %>% 
  filter(ihmp_sample_id == "CSM79HGF_P" | ihmp_sample_id == "CSM79HPK") # notice urls extensions 
```

```
## # A tibble: 3 × 2
##   ihmp_sample_id urls                                                           
##   <chr>          <chr>                                                          
## 1 CSM79HGF_P     https://downloads.hmpdacc.org/ihmp/ibd/genome/microbiome/wgs/r…
## 2 CSM79HPK       https://downloads.hmpdacc.org/ihmp/ibd/genome/microbiome/wgs/r…
## 3 CSM79HPK       https://downloads.hmpdacc.org/ihmp/ibd/genome/microbiome/wgs/r…
```

```r
full_file_manifest <- full_file_manifest %>%
  separate(urls, into = c(NA, "ihmp_sample_id", NA), sep = "/raw/|.tar|.fastq.gz", remove = F) %>% 
  relocate(ihmp_sample_id)
full_file_manifest
```

```
## # A tibble: 2,341 × 14
##    ihmp_sample_id sample_id   file_id md5     size urls  subject_id subject_uuid
##    <chr>          <chr>       <chr>   <chr>  <dbl> <chr>      <dbl> <chr>       
##  1 CSM5MCXD       1419f08f55… 1419f0… 6c83… 1.05e9 http…       3001 1419f08f554…
##  2 CSM5MCXD       1419f08f55… 56bda9… 1eb2… 3.60e6 http…       3001 1419f08f554…
##  3 CSM5MCYS       1419f08f55… 1419f0… 5334… 1.32e9 http…       3001 1419f08f554…
##  4 CSM67U9J       1419f08f55… 1419f0… e8cb… 9.81e8 http…       3001 1419f08f554…
##  5 CSM67UA2       1419f08f55… 1419f0… b5d6… 1.03e9 http…       3001 1419f08f554…
##  6 CSM67UA2       1419f08f55… 56bda9… ce38… 8.50e7 http…       3001 1419f08f554…
##  7 CSM67UGC       1419f08f55… 1419f0… f652… 2.21e8 http…       3001 1419f08f554…
##  8 CSM79HG5       1419f08f55… 1419f0… 2df8… 1.56e9 http…       3001 1419f08f554…
##  9 CSM79HGP       1419f08f55… 1419f0… 3f87… 1.49e9 http…       3001 1419f08f554…
## 10 CSM79HGP       1419f08f55… 56bda9… 194c… 1.34e8 http…       3001 1419f08f554…
## # … with 2,331 more rows, and 6 more variables: sample_body_site <chr>,
## #   visit_number <dbl>, subject_gender <chr>, subject_race <chr>,
## #   study_full_name <chr>, project_name <chr>
```


### merging with ihmp

reminder, we are here

```r
file_mani <- readr::read_tsv("metadata-of-2342-samples/hmp_manifest_47683f7597.tsv")
```

```
## Rows: 2341 Columns: 5
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (4): file_id, md5, urls, sample_id
## dbl (1): size
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
samp_id_mani <- readr::read_tsv("metadata-of-2342-samples/hmp_manifest_metadata_b7ce70666.tsv")
```

```
## Rows: 2341 Columns: 9
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (7): sample_id, subject_uuid, sample_body_site, subject_gender, subject_...
## dbl (2): subject_id, visit_number
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
full_file_manifest <- left_join(x = file_mani, y = unique(samp_id_mani))%>%
  separate(urls, into = c(NA, "ihmp_sample_id", NA), sep = "/raw/|.tar|.fastq.gz", remove = F) %>% 
  relocate(ihmp_sample_id, sample_id, file_id) %>% 
  arrange(ihmp_sample_id, sample_id, file_id)
```

```
## Joining, by = "sample_id"
```

```r
full_file_manifest
```

```
## # A tibble: 2,341 × 14
##    ihmp_sample_id sample_id   file_id md5     size urls  subject_id subject_uuid
##    <chr>          <chr>       <chr>   <chr>  <dbl> <chr>      <dbl> <chr>       
##  1 CSM5FZ3N_P     cad5930efe… 599427… 0081… 1.58e9 http…       3001 1419f08f554…
##  2 CSM5FZ3R_P     5994274eae… 599427… 93ee… 2.61e9 http…       3001 1419f08f554…
##  3 CSM5FZ3T_P     cad5930efe… 599427… 1331… 1.72e9 http…       3002 1419f08f554…
##  4 CSM5FZ3V_P     5994274eae… 599427… 291b… 2.85e9 http…       3001 1419f08f554…
##  5 CSM5FZ3X_P     5994274eae… 599427… 8bf2… 2.24e9 http…       3002 1419f08f554…
##  6 CSM5FZ3Z_P     5994274eae… 599427… c146… 3.82e9 http…       3002 1419f08f554…
##  7 CSM5FZ42_P     5994274eae… 599427… 6659… 5.11e9 http…       3002 1419f08f554…
##  8 CSM5FZ44_P     cad5930efe… 599427… bdde… 3.01e9 http…       3002 1419f08f554…
##  9 CSM5FZ46_P     5994274eae… 599427… 4632… 2.45e9 http…       3002 1419f08f554…
## 10 CSM5FZ4A_P     cad5930efe… 599427… 9044… 2.96e9 http…       3004 1419f08f554…
## # … with 2,331 more rows, and 6 more variables: sample_body_site <chr>,
## #   visit_number <dbl>, subject_gender <chr>, subject_race <chr>,
## #   study_full_name <chr>, project_name <chr>
```

lets try merging `full_file_manifest` with `ihmp_metadata` on `ihmp_sample_id`
*lol it didnt work

```r
ihmp <- readr::read_csv("../ibdmdb_relevant_participantIDs_and_sampleIDs/paired_mgx_vrx_ihmp_metadata.csv")%>% 
  rename("ihmp_sample_id" = external_id)
```

```
## Warning: One or more parsing issues, see `problems()` for details
```

```
## Rows: 2299 Columns: 490
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: ","
## chr  (145): project, external_id, participant_id, site_sub_coll, data_type, ...
## dbl   (20): week_num, interval_days, visit_num, gssr_i_ds, wr_id, number_lan...
## lgl  (324): reads_ribosomal, delta, has_the_subject_had_a_cholecystectomy, h...
## date   (1): date_of_receipt
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
#ihmp
#
```


```r
nrow(ihmp)
```

```
## [1] 2299
```

```r
nrow(full_file_manifest)
```

```
## [1] 2341
```

```r
n_distinct(ihmp$ihmp_sample_id)
```

```
## [1] 1614
```

```r
n_distinct(full_file_manifest$ihmp_sample_id)
```

```
## [1] 1656
```



```r
df <- left_join(ihmp, full_file_manifest, "ihmp_sample_id")
#right_join(ihmp, full_file_manifest, "ihmp_sample_id")
```


```r
sra_mani <- readr::read_csv("sra-manifest.csv") 
```

```
## Rows: 378 Columns: 17
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: ","
## chr (12): Experiment Accession, Experiment Title, Organism Name, Instrument,...
## dbl  (4): Total Size, Mb, Total RUNs, Total Spots, Total Bases
## lgl  (1): Sample Title
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
sra_mani <- clean_names(sra_mani)

head(sra_mani, 2)
```

```
## # A tibble: 2 × 17
##   experiment_accession experiment_title       organism_name instrument submitter
##   <chr>                <chr>                  <chr>         <chr>      <chr>    
## 1 SRX3139489           fecal sample from HMP… human gut me… Illumina … Harvard …
## 2 SRX3139488           fecal sample from HMP… human gut me… Illumina … Harvard …
## # … with 12 more variables: study_accession <chr>, study_title <chr>,
## #   sample_accession <chr>, sample_title <lgl>, total_size_mb <dbl>,
## #   total_ru_ns <dbl>, total_spots <dbl>, total_bases <dbl>,
## #   library_name <chr>, library_strategy <chr>, library_source <chr>,
## #   library_selection <chr>
```

```r
sra_mani <- sra_mani %>% 
  separate(experiment_title, c(NA,"ihmp_sample_id",NA), sep = "Subject |_MG", remove = F) %>% 
  relocate(ihmp_sample_id)
```

```
## Warning: Expected 3 pieces. Missing pieces filled with `NA` in 78 rows [223,
## 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
## 240, 241, 242, ...].
```

```r
sra_mani
```

```
## # A tibble: 378 × 18
##    ihmp_sample_id experiment_accession experiment_title organism_name instrument
##    <chr>          <chr>                <chr>            <chr>         <chr>     
##  1 H4020C2        SRX3139489           fecal sample fr… human gut me… Illumina …
##  2 C3011C6        SRX3139488           fecal sample fr… human gut me… Illumina …
##  3 H4022C1        SRX3139487           fecal sample fr… human gut me… Illumina …
##  4 C3017C3        SRX3139486           fecal sample fr… human gut me… Illumina …
##  5 C3011C5        SRX3139485           fecal sample fr… human gut me… Illumina …
##  6 C3017C1        SRX3139484           fecal sample fr… human gut me… Illumina …
##  7 C3013C5        SRX3139483           fecal sample fr… human gut me… Illumina …
##  8 C3011C7        SRX3139482           fecal sample fr… human gut me… Illumina …
##  9 C3010C4        SRX3139481           fecal sample fr… human gut me… Illumina …
## 10 C3011C2        SRX3139480           fecal sample fr… human gut me… Illumina …
## # … with 368 more rows, and 13 more variables: submitter <chr>,
## #   study_accession <chr>, study_title <chr>, sample_accession <chr>,
## #   sample_title <lgl>, total_size_mb <dbl>, total_ru_ns <dbl>,
## #   total_spots <dbl>, total_bases <dbl>, library_name <chr>,
## #   library_strategy <chr>, library_source <chr>, library_selection <chr>
```


```r
#ihmp_sample_id, experiment_accession, study_accession, sample_accession, library_name
#names(df)

trim_df <- df %>% select(ihmp_sample_id, participant_id, data_type, week_num, interval_days, project) %>% arrange(participant_id, data_type, week_num)

trim_sra <- sra_mani %>% select(ihmp_sample_id, experiment_accession, study_accession, sample_accession, library_name) %>% arrange(ihmp_sample_id)

trim_df
```

```
## # A tibble: 3,669 × 6
##    ihmp_sample_id participant_id data_type    week_num interval_days project
##    <chr>          <chr>          <chr>           <dbl>         <dbl> <chr>  
##  1 CSM5FZ3N_P     C3001          metagenomics        0             0 G79889 
##  2 CSM5FZ3R_P     C3001          metagenomics        2            14 G79894 
##  3 CSM5YRY7_P     C3001          metagenomics        4            18 G79903 
##  4 CSM5FZ3V_P     C3001          metagenomics        6            13 G79913 
##  5 CSM5FZ4C_P     C3001          metagenomics        8            11 G79926 
##  6 CSM5MCVD_P     C3001          metagenomics       12            20 G79969 
##  7 CSM5MCVF_P     C3001          metagenomics       14            15 G80015 
##  8 CSM5MCVV_P     C3001          metagenomics       16            15 G79979 
##  9 CSM5MCWI_P     C3001          metagenomics       18            12 G80050 
## 10 CSM5MCXD       C3001          metagenomics       20            14 G110156
## # … with 3,659 more rows
```

```r
trim_sra
```

```
## # A tibble: 378 × 5
##    ihmp_sample_id experiment_acce… study_accession sample_accession library_name
##    <chr>          <chr>            <chr>           <chr>            <chr>       
##  1 C3001C1        SRX3139291       SRP115812       SRS2472163       G79889_MGX  
##  2 C3001C10       SRX3106607       SRP115812       SRS2441858       G80053_MGX  
##  3 C3001C10_MTX   SRX3106690       SRP115812       SRS2441941       G89370_MTX  
##  4 C3001C2        SRX3139294       SRP115812       SRS2472166       G79894_MGX  
##  5 C3001C3        SRX3139288       SRP115812       SRS2472160       G79903_MGX  
##  6 C3001C4        SRX3139289       SRP115812       SRS2472161       G79913_MGX  
##  7 C3001C5        SRX3106633       SRP115812       SRS2441884       G79931_MGX  
##  8 C3001C5_MTX    SRX3106735       SRP115812       SRS2441986       G89308_MTX  
##  9 C3001C7        SRX3139397       SRP115812       SRS2472234       G79969_MGX  
## 10 C3001C8        SRX3139360       SRP115812       SRS2472871       G80015_MGX  
## # … with 368 more rows
```

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("sradb")
```

```
## Bioconductor version 3.14 (BiocManager 1.30.16), R 4.1.2 (2021-11-01)
```

```
## Installing package(s) 'sradb'
```

```
## Warning in .inet_warning(msg): package 'sradb' is not available for Bioconductor version '3.14'
## 
## A version of this package for your version of R might be available elsewhere,
## see the ideas at
## https://cran.r-project.org/doc/manuals/r-patched/R-admin.html#Installing-packages
```

```
## Warning in .inet_warning(msg): Perhaps you meant 'SRAdb' ?
```

```
## Old packages: 'admisc', 'commonmark', 'desc', 'maptools', 'openssl', 'Rcpp',
##   'rlang', 'rmarkdown', 'sf'
```



```r
inner_join(trim_sra, trim_df, "ihmp_sample_id")# %>% 
```

```
## # A tibble: 0 × 10
## # … with 10 variables: ihmp_sample_id <chr>, experiment_accession <chr>,
## #   study_accession <chr>, sample_accession <chr>, library_name <chr>,
## #   participant_id <chr>, data_type <chr>, week_num <dbl>, interval_days <dbl>,
## #   project <chr>
```

```r
#   relocate(sample_accession, participant_id, data_type, week_num, interval_days, #ihmp_sample_id, project, experiment_accession, study_accession, library_name) %>% 
#  arrange("participant_id, data_type, week_num,")
#
```



### dear god, I guess I need software now


somehow, this installed  [sradb](https://www.bioconductor.org/packages/release/bioc/html/SRAdb.html) (repo [here](https://github.com/ncbi/sra-tools)) 

```r
#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#
#BiocManager::install("SRAdb")

library("SRAdb")
```

```
## Loading required package: RSQLite
```

```
## Loading required package: graph
```

```
## Loading required package: BiocGenerics
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:dplyr':
## 
##     combine, intersect, setdiff, union
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
##     union, unique, unsplit, which.max, which.min
```

```
## 
## Attaching package: 'graph'
```

```
## The following object is masked from 'package:stringr':
## 
##     boundary
```

```
## Loading required package: RCurl
```

```
## 
## Attaching package: 'RCurl'
```

```
## The following object is masked from 'package:tidyr':
## 
##     complete
```

```
## Setting options('download.file.method.GEOquery'='auto')
```

```
## Setting options('GEOquery.inmemory.gpl'=FALSE)
```







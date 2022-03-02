---
title: "merging metadata with sras"
output: 
  html_document: 
    keep_md: yes
date: '2022-03-02'
---



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
# has fileID, md5sum, file size, download urls,  and sample ID
df1<- readr::read_csv("metadata-of-2342-samples/hmp_manifest_47683f7597.tsv")
```

```
## Rows: 2341 Columns: 1
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: ","
## chr (1): file_id	md5	size	urls	sample_id
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
# has: sample_id	subject_id	subject_uuid	sample_body_site	visit_number	subject_gender	subject_race	study_full_name	project_name
df2 <- readr::read_csv("metadata-of-2342-samples/hmp_manifest_metadata_b7ce70666.tsv")
```

```
## Rows: 2341 Columns: 1
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: ","
## chr (1): sample_id	subject_id	subject_uuid	sample_body_site	visit_number	sub...
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```


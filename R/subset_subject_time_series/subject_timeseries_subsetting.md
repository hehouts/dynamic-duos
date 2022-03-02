---
title: "Subject Subsetting"
output: 
  html_document: 
    keep_md: yes
---





## Load libraries

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

```
## 
## Attaching package: 'janitor'
```

```
## The following objects are masked from 'package:stats':
## 
##     chisq.test, fisher.test
```


## load data 
This .csv was pulled from the ibdmdb webpage?

```r
setwd("/Users/hehouts/projects/dynamic-duos/R")

reference_metadata <- readr::read_csv("../metadata/ibdmdb-metadata-manifests/ihmp_metadata.csv", show_col_types = F)
```

```
## Warning: One or more parsing issues, see `problems()` for details
```

```r
#tidy names
metadata <- clean_names(reference_metadata)

head(metadata, 3)
```

```
## # A tibble: 3 × 490
##   project      external_id participant_id site_sub_coll data_type  week_num
##   <chr>        <chr>       <chr>          <chr>         <chr>         <dbl>
## 1 C3001CSC1_BP 206615      C3001          C3001CSC1     biopsy_16S        2
## 2 C3001CSC2_BP 206614      C3001          C3001CSC2     biopsy_16S        2
## 3 C3002CSC1_BP 206617      C3002          C3002CSC1     biopsy_16S        0
## # … with 484 more variables: date_of_receipt <date>, interval_days <dbl>,
## #   visit_num <dbl>, research_project <chr>, pdo_number <chr>, gssr_i_ds <dbl>,
## #   product <chr>, lcset <chr>, aggregated_lanes <chr>, wr_id <dbl>,
## #   number_lanes_in_aggregation <dbl>, reads_raw <dbl>, reads_filtered <dbl>,
## #   reads_qc_fail <dbl>, reads_human <dbl>, reads_ribosomal <dbl>,
## #   reads_viral <dbl>, delta <lgl>, interval_name <chr>,
## #   interval_sequence <dbl>, project_specific_id <dbl>, site_name <chr>, …
```

## Three subject subset

```r
subset <- metadata%>%
  filter(participant_id == "C3002"
        | participant_id == "C3003"
        | participant_id == "C3008"
        )%>%
  arrange(participant_id, week_num)
subset
```

```
## # A tibble: 151 × 490
##    project      external_id participant_id site_sub_coll data_type      week_num
##    <chr>        <chr>       <chr>          <chr>         <chr>             <dbl>
##  1 C3002CSC1_BP 206617      C3002          C3002CSC1     biopsy_16S            0
##  2 C3002CSC2_BP 206619      C3002          C3002CSC2     biopsy_16S            0
##  3 C3002CSC3_BP 206616      C3002          C3002CSC3     biopsy_16S            0
##  4 C3002CSC4_BP 206618      C3002          C3002CSC4     biopsy_16S            0
##  5 C3002CSC1_TX CSM5FZ1G    C3002          C3002CSC1     host_transcri…        0
##  6 C3002CSC2_TX CSM5FZ1F    C3002          C3002CSC2     host_transcri…        0
##  7 C3002CSC3_TX CSM5FZ1I    C3002          C3002CSC3     host_transcri…        0
##  8 C3002CSC4_TX CSM5FZ1H    C3002          C3002CSC4     host_transcri…        0
##  9 C3002C1_MBX  CSM5FZ3T    C3002          C3002C1       metabolomics          0
## 10 G79914       CSM5FZ3T_P  C3002          C3002C1       metagenomics          0
## # … with 141 more rows, and 484 more variables: date_of_receipt <date>,
## #   interval_days <dbl>, visit_num <dbl>, research_project <chr>,
## #   pdo_number <chr>, gssr_i_ds <dbl>, product <chr>, lcset <chr>,
## #   aggregated_lanes <chr>, wr_id <dbl>, number_lanes_in_aggregation <dbl>,
## #   reads_raw <dbl>, reads_filtered <dbl>, reads_qc_fail <dbl>,
## #   reads_human <dbl>, reads_ribosomal <dbl>, reads_viral <dbl>, delta <lgl>,
## #   interval_name <chr>, interval_sequence <dbl>, project_specific_id <dbl>, …
```
## metagenomics data


```r
subset_mgx <- subset%>%
  filter(data_type=="metagenomics") %>%
  arrange(participant_id, week_num)
subset_mgx
```

```
## # A tibble: 38 × 490
##    project external_id participant_id site_sub_coll data_type    week_num
##    <chr>   <chr>       <chr>          <chr>         <chr>           <dbl>
##  1 G79914  CSM5FZ3T_P  C3002          C3002C1       metagenomics        0
##  2 G79929  CSM5FZ3X_P  C3002          C3002C2       metagenomics        2
##  3 G79949  CSM5FZ3Z_P  C3002          C3002C3       metagenomics        4
##  4 G79970  CSM5FZ42_P  C3002          C3002C4       metagenomics        6
##  5 G80016  CSM5FZ44_P  C3002          C3002C5       metagenomics        8
##  6 G79994  CSM5FZ46_P  C3002          C3002C6       metagenomics       10
##  7 G80049  CSM5MCVJ_P  C3002          C3002C7       metagenomics       12
##  8 G110038 CSM5MCVL    C3002          C3002C8       metagenomics       14
##  9 G110144 CSM5MCVN    C3002          C3002C9       metagenomics       16
## 10 G110061 CSM67UBF    C3002          C3002C13      metagenomics       24
## # … with 28 more rows, and 484 more variables: date_of_receipt <date>,
## #   interval_days <dbl>, visit_num <dbl>, research_project <chr>,
## #   pdo_number <chr>, gssr_i_ds <dbl>, product <chr>, lcset <chr>,
## #   aggregated_lanes <chr>, wr_id <dbl>, number_lanes_in_aggregation <dbl>,
## #   reads_raw <dbl>, reads_filtered <dbl>, reads_qc_fail <dbl>,
## #   reads_human <dbl>, reads_ribosomal <dbl>, reads_viral <dbl>, delta <lgl>,
## #   interval_name <chr>, interval_sequence <dbl>, project_specific_id <dbl>, …
```




## vrmx data


```r
subset_vrmx <- metadata%>%
  filter(data_type=="viromics") %>%
  arrange(participant_id, week_num)
subset_vrmx
```

```
## # A tibble: 703 × 490
##    project      external_id participant_id site_sub_coll data_type week_num
##    <chr>        <chr>       <chr>          <chr>         <chr>        <dbl>
##  1 C3001C11_MVX CSM5MCXD    C3001          C3001C11      viromics        20
##  2 C3001C15_MVX CSM67UA2    C3001          C3001C15      viromics        28
##  3 C3001C20_MVX CSM79HGP    C3001          C3001C20      viromics        38
##  4 C3002C9_MVX  CSM5MCVN    C3002          C3002C9       viromics        16
##  5 C3002C14_MVX CSM67UBH    C3002          C3002C14      viromics        26
##  6 C3002C17_MVX CSM67UBN    C3002          C3002C17      viromics        32
##  7 C3002C19_MVX CSM79HJW    C3002          C3002C19      viromics        36
##  8 C3003C6_MVX  CSM5FZ4M    C3003          C3003C6       viromics         9
##  9 C3003C9_MVX  CSM5MCWQ    C3003          C3003C9       viromics        16
## 10 C3003C14_MVX CSM67UBZ    C3003          C3003C14      viromics        26
## # … with 693 more rows, and 484 more variables: date_of_receipt <date>,
## #   interval_days <dbl>, visit_num <dbl>, research_project <chr>,
## #   pdo_number <chr>, gssr_i_ds <dbl>, product <chr>, lcset <chr>,
## #   aggregated_lanes <chr>, wr_id <dbl>, number_lanes_in_aggregation <dbl>,
## #   reads_raw <dbl>, reads_filtered <dbl>, reads_qc_fail <dbl>,
## #   reads_human <dbl>, reads_ribosomal <dbl>, reads_viral <dbl>, delta <lgl>,
## #   interval_name <chr>, interval_sequence <dbl>, project_specific_id <dbl>, …
```


## One subject, mgx & vrmx combined

```r
c3002_metadata <- subset%>%
  filter(participant_id == "C3002")%>%
  filter(data_type=="viromics" | data_type=="metagenomics")%>%
  arrange(participant_id, data_type, week_num)%>%
  select(external_id, week_num, participant_id, data_type)

c3002_metadata
```

```
## # A tibble: 19 × 4
##    external_id week_num participant_id data_type   
##    <chr>          <dbl> <chr>          <chr>       
##  1 CSM5FZ3T_P         0 C3002          metagenomics
##  2 CSM5FZ3X_P         2 C3002          metagenomics
##  3 CSM5FZ3Z_P         4 C3002          metagenomics
##  4 CSM5FZ42_P         6 C3002          metagenomics
##  5 CSM5FZ44_P         8 C3002          metagenomics
##  6 CSM5FZ46_P        10 C3002          metagenomics
##  7 CSM5MCVJ_P        12 C3002          metagenomics
##  8 CSM5MCVL          14 C3002          metagenomics
##  9 CSM5MCVN          16 C3002          metagenomics
## 10 CSM67UBF          24 C3002          metagenomics
## 11 CSM67UBH          26 C3002          metagenomics
## 12 CSM67UBN          32 C3002          metagenomics
## 13 CSM67UBR          34 C3002          metagenomics
## 14 CSM79HJW          36 C3002          metagenomics
## 15 CSM79HJY          38 C3002          metagenomics
## 16 CSM5MCVN          16 C3002          viromics    
## 17 CSM67UBH          26 C3002          viromics    
## 18 CSM67UBN          32 C3002          viromics    
## 19 CSM79HJW          36 C3002          viromics
```


```r
c3002_metadata %>% select(external_id)
```

```
## # A tibble: 19 × 1
##    external_id
##    <chr>      
##  1 CSM5FZ3T_P 
##  2 CSM5FZ3X_P 
##  3 CSM5FZ3Z_P 
##  4 CSM5FZ42_P 
##  5 CSM5FZ44_P 
##  6 CSM5FZ46_P 
##  7 CSM5MCVJ_P 
##  8 CSM5MCVL   
##  9 CSM5MCVN   
## 10 CSM67UBF   
## 11 CSM67UBH   
## 12 CSM67UBN   
## 13 CSM67UBR   
## 14 CSM79HJW   
## 15 CSM79HJY   
## 16 CSM5MCVN   
## 17 CSM67UBH   
## 18 CSM67UBN   
## 19 CSM79HJW
```



```
CSM5FZ3T_P.tar			
CSM5FZ3X_P.tar			
CSM5FZ3Z_P.tar			
CSM5FZ42_P.tar			
CSM5FZ44_P.tar			
CSM5FZ46_P.tar			
CSM5MCVJ_P.tar			
CSM5MCVL.tar		
CSM5MCVN.tar		
CSM67UBF.tar
CSM67UBH.tar				
CSM67UBN.tar				
CSM67UBR.tar				
CSM79HJW.tar				
CSM79HJY.tar

CSM5MCVN.tar				
CSM67UBH.tar				
CSM67UBN.tar				
CSM79HJW.tar
```

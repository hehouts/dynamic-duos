---
title: "Metadata Summary"
output: 
  html_document: 
    keep_md: yes
date: '2022-03-16'
---

### Definitions 
MVX - viromics, Viral fraction
MGX - metagenomic, Bulk fraction


### Chunk 1 - Load Libraries and Data 

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
library(naniar)

#setwd("/Users/hehouts/projects/dynamic-duos/R")

ibdmdb <- clean_names(readr::read_csv("csvs_and_other_metadata/ihmp_metadata.csv", show_col_types = F))
ibdmdb <- ibdmdb %>% 
  filter(data_type == "metagenomics" | data_type == "viromics") %>% 
  arrange(participant_id, data_type, week_num)
```


### Chunk 2 - Determine Number of Metagenome and Viromes Samples
Conclusions:
* 2341 Metagenome & Virome samples from the iHMP IBD Study
  + 1638 Metagenomes 
  + 703 Viromes. 



```r
# 2341 total samples, mgx + mvx
ibdmdb <- ibdmdb %>% 
  filter(data_type == "metagenomics" | data_type == "viromics") %>% 
  relocate(external_id, tube_a_viromics, week_num, participant_id, data_type) %>% 
  arrange(participant_id, data_type, week_num)

# 1638 MGX samples
ibdmdb_mgx_only <- ibdmdb %>% 
  filter(data_type == "metagenomics") %>% 
  arrange(participant_id, week_num)

# 703 MVX samples
ibdmdb_mvx_only <- ibdmdb %>% 
  filter(data_type == "viromics") %>% 
  arrange(participant_id, week_num)

# hard saved csv of all mgx & vrx samples in the ibdmdb, no other 'omics. should be 2341 long. 
write_csv(ibdmdb, "csvs_and_other_metadata/ibdmdb_metadata_2341.csv")
```

### Chunk 3 - Determine Number of participants with both MVX and MGX Samples, and Determine Which Have Antibiotic Data

Conclusions:
* 130 total participants in the study
* 105 participants have both virome and metagenome samples
  + 2299 samples are from those 105 participants
* 34 participants that have viromes samples, and had antibiotics
  + 838 samples are from those 34 participants


```r
# 130 total participants
ibdmdb %>% distinct(participant_id)
```

```
## # A tibble: 130 × 1
##    participant_id
##    <chr>         
##  1 C3001         
##  2 C3002         
##  3 C3003         
##  4 C3004         
##  5 C3005         
##  6 C3006         
##  7 C3007         
##  8 C3008         
##  9 C3009         
## 10 C3010         
## # … with 120 more rows
```

```r
# 105 participant ids of ONLY participants that have vrx samples, regardless of antibiotic status, i.e. "pair" participants
pair_ids_vec <- ibdmdb %>% 
  select(data_type, participant_id) %>% 
  filter(data_type == "viromics") %>% 
  distinct(participant_id) %>% 
  pull(participant_id)

# 2299 samples from just the 105 "pair" participants
pairs_df <- ibdmdb %>% 
  filter(data_type == "metagenomics" | data_type == "viromics") %>% 
  filter(participant_id %in% pair_ids_vec) %>%
  arrange(participant_id, data_type, week_num)
#write_csv(pairs_df, "csvs_and_other_metadata/ihmp_metadata_2299_paired_mgx_vrx.csv")


### Antibiotics:

# 34 participant ids of ONLY participants that have vrx samples, and had antibiotics, i.e. "abxpair", "relevant participants"
abx_pair_ids_vec <- ibdmdb %>% 
  select(data_type, participant_id, antibiotics) %>% 
  filter(data_type == "viromics" & antibiotics=="Yes") %>%
  distinct(participant_id) %>% 
  pull(participant_id)


# 838 samples from just the 34 "abxpair" participants,
abx_pairs_df <- ibdmdb %>% 
  filter(data_type == "metagenomics" | data_type == "viromics") %>% 
  filter(participant_id %in% abx_pair_ids_vec) %>%
  arrange(participant_id, data_type, week_num)
#write_csv(abx_pairs_df, "csvs_and_other_metadata/ihmp_metadata_838_antibiotics_paired_mgx_vrx.csv")
```



### Chunk 4 - Finding existing filenames for viromes in the ibdmdb metadata

After some sleuthing, it was determined that the existing virome filenames, which start with the prefix `SM-`, correspond with the `tube_a_viromics` column

```r
ibdmdb %>% 
  select(contains("tube"))
```

```
## # A tibble: 2,341 × 14
##    tube_a_viromics serum_tube_number_1_receiv… serum_tubes_num… number_of_flora…
##    <chr>           <chr>                       <chr>                       <dbl>
##  1 SM-5QVZE        No                          No                             NA
##  2 SM-5UKLB        No                          No                             NA
##  3 SM-5YSLQ        No                          No                             NA
##  4 SM-61GG5        No                          No                             NA
##  5 SM-64NAA        No                          No                             NA
##  6 SM-6DQVE        No                          No                             NA
##  7 SM-6IEV2        No                          No                             NA
##  8 SM-6MYZ3        No                          No                             NA
##  9 SM-6QU1N        No                          No                             NA
## 10 SM-6UG78        No                          No                             NA
## # … with 2,331 more rows, and 10 more variables:
## #   number_of_tubes_collected_for_epithelial_cell_biopsies <dbl>,
## #   number_of_dna_rna_tubes_collected <dbl>, tube_a_dna_rna <chr>,
## #   tube_a_metabolomics <chr>, tube_a_storage <chr>,
## #   tube_b_fecal_calprotectin <chr>, tube_b_proteomics <chr>,
## #   stool_sample_id_tube_a_et_oh <chr>, sample_id_tube_b_no_preservative <chr>,
## #   tube_a_and_b_received_at_broad <chr>
```





### Chunk 5 - s5 subset: File IDs for viromes
I want to have the viromes and metagenomes for the following 5 participants. (2 have abx, 3 dont)

```
C3001
C3002
C3003
C3004
C3005
```


```r
ibdmdb %>% 
  filter(participant_id %in% c("C3001", "C3002", "C3003", "C3004", "C3005")) %>%
  mutate(abx_group = case_when(
      participant_id %in% abx_pair_ids_vec~T,
      !(participant_id %in% abx_pair_ids_vec)~F
      )) %>%
  select(participant_id, abx_group) %>% 
  unique()
```

```
## # A tibble: 5 × 2
##   participant_id abx_group
##   <chr>          <lgl>    
## 1 C3001          FALSE    
## 2 C3002          TRUE     
## 3 C3003          TRUE     
## 4 C3004          FALSE    
## 5 C3005          FALSE
```

```r
subset5 <- ibdmdb %>% 
  filter(participant_id %in% c("C3001", "C3002", "C3003", "C3004", "C3005")) %>%
  mutate(abx_group = case_when(
      participant_id %in% abx_pair_ids_vec~T,
      !(participant_id %in% abx_pair_ids_vec)~F
      )) %>% 
  relocate(external_id, tube_a_viromics, week_num, participant_id, abx_group, data_type) %>% 
  select(1:15)
subset5
```

```
## # A tibble: 109 × 15
##    external_id tube_a_viromics week_num participant_id abx_group data_type   
##    <chr>       <chr>              <dbl> <chr>          <lgl>     <chr>       
##  1 CSM5FZ3N_P  SM-5QVZE               0 C3001          FALSE     metagenomics
##  2 CSM5FZ3R_P  SM-5UKLB               2 C3001          FALSE     metagenomics
##  3 CSM5YRY7_P  SM-5YSLQ               4 C3001          FALSE     metagenomics
##  4 CSM5FZ3V_P  SM-61GG5               6 C3001          FALSE     metagenomics
##  5 CSM5FZ4C_P  SM-64NAA               8 C3001          FALSE     metagenomics
##  6 CSM5MCVD_P  SM-6DQVE              12 C3001          FALSE     metagenomics
##  7 CSM5MCVF_P  SM-6IEV2              14 C3001          FALSE     metagenomics
##  8 CSM5MCVV_P  SM-6MYZ3              16 C3001          FALSE     metagenomics
##  9 CSM5MCWI_P  SM-6QU1N              18 C3001          FALSE     metagenomics
## 10 CSM5MCXD    SM-6UG78              20 C3001          FALSE     metagenomics
## # … with 99 more rows, and 9 more variables: project <chr>,
## #   site_sub_coll <chr>, date_of_receipt <date>, interval_days <dbl>,
## #   visit_num <dbl>, research_project <chr>, pdo_number <chr>, gssr_i_ds <dbl>,
## #   product <chr>
```

```r
subset5_virome_filenames <- subset5 %>% 
  filter(data_type == "viromics")
```


### Chunk 6 - uids - a new ibdmdb table with download filenames and uid fro each sample


```r
ibdmdb_uids <- ibdmdb %>% 
  mutate(abx_group = case_when(
      participant_id %in% abx_pair_ids_vec~T,
      !(participant_id %in% abx_pair_ids_vec)~F
      )) %>% 
  relocate(external_id, tube_a_viromics, week_num, participant_id, abx_group, data_type
      ) %>% 
  select(1:12) %>% 
  mutate(file_prefix = case_when(
    (data_type == "viromics")~tube_a_viromics,
    (data_type == "metagenomics")~external_id
    )
  ) %>% 
  mutate(zip_type = case_when(
        (data_type == "viromics")~"VIR",
        (data_type == "metagenomics" &
           grepl('_P$', external_id) == T)~"P",
        (data_type == "metagenomics" &
           grepl('_TR$', external_id) == T)~"TR",
        (data_type == "metagenomics" &
           !grepl('_TR$', external_id) &
           !grepl('_P$', external_id) == T)~"STD"
    )) %>% 
mutate(uid = case_when(
    (data_type == "metagenomics")~external_id,
    (data_type == "viromics")~tube_a_viromics,
    )) %>% 
  mutate(file_suffix = case_when(
    (zip_type == "VIR")~".tar",
    (zip_type == "STD")~".tar",
    (zip_type == "P")~".fastq.gz",
    (zip_type == "TR" )~".tar",
    )) %>%
  relocate(file_prefix, file_suffix) %>% 
  unite("file", file_prefix:file_suffix, sep = "", remove = F) %>% 
  relocate(file, zip_type, uid)
ibdmdb_uids
```

```
## # A tibble: 2,341 × 17
##    file       zip_type uid   file_prefix file_suffix external_id tube_a_viromics
##    <chr>      <chr>    <chr> <chr>       <chr>       <chr>       <chr>          
##  1 CSM5FZ3N_… P        CSM5… CSM5FZ3N_P  .fastq.gz   CSM5FZ3N_P  SM-5QVZE       
##  2 CSM5FZ3R_… P        CSM5… CSM5FZ3R_P  .fastq.gz   CSM5FZ3R_P  SM-5UKLB       
##  3 CSM5YRY7_… P        CSM5… CSM5YRY7_P  .fastq.gz   CSM5YRY7_P  SM-5YSLQ       
##  4 CSM5FZ3V_… P        CSM5… CSM5FZ3V_P  .fastq.gz   CSM5FZ3V_P  SM-61GG5       
##  5 CSM5FZ4C_… P        CSM5… CSM5FZ4C_P  .fastq.gz   CSM5FZ4C_P  SM-64NAA       
##  6 CSM5MCVD_… P        CSM5… CSM5MCVD_P  .fastq.gz   CSM5MCVD_P  SM-6DQVE       
##  7 CSM5MCVF_… P        CSM5… CSM5MCVF_P  .fastq.gz   CSM5MCVF_P  SM-6IEV2       
##  8 CSM5MCVV_… P        CSM5… CSM5MCVV_P  .fastq.gz   CSM5MCVV_P  SM-6MYZ3       
##  9 CSM5MCWI_… P        CSM5… CSM5MCWI_P  .fastq.gz   CSM5MCWI_P  SM-6QU1N       
## 10 CSM5MCXD.… STD      CSM5… CSM5MCXD    .tar        CSM5MCXD    SM-6UG78       
## # … with 2,331 more rows, and 10 more variables: week_num <dbl>,
## #   participant_id <chr>, abx_group <lgl>, data_type <chr>, project <chr>,
## #   site_sub_coll <chr>, date_of_receipt <date>, interval_days <dbl>,
## #   visit_num <dbl>, research_project <chr>
```

```r
ibdmdb_uids_s5 <- ibdmdb_uids %>% filter(participant_id %in% c("C3001", "C3002", "C3003", "C3004", "C3005"))

ibdmdb_uids %>% select(file, file_prefix, zip_type, uid) %>% 
  write_csv("csvs_and_other_metadata/file_prefix_zip_uid.csv")

ibdmdb_uids %>% 
  filter(participant_id %in% c("C3001", "C3002", "C3003", "C3004", "C3005")) %>%
  select(file, file_prefix, zip_type, uid) %>%
  write_csv("csvs_and_other_metadata/file_prefix_zip_uid_s5.csv")
```


### Chunk 7 make the filename.txt's for all the zip subtypes


```r
# VIR
ibdmdb_uids %>% filter(zip_type == "VIR") %>% select(file) %>%  
  write_csv("../snakecrate/data/raw_data/vir/vir_filenames.txt", col_names = F)

ibdmdb_uids_s5  %>% filter(zip_type == "VIR") %>% select(file) %>%  
  write_csv("../snakecrate/data/raw_data/vir/s5_vir_filenames.txt", col_names = F)
######


# TR
ibdmdb_uids %>% filter(zip_type == "TR") %>% select(file) %>%  
  write_csv("../snakecrate/data/raw_data/tr/tr_filenames.txt", col_names = F)

ibdmdb_uids_s5  %>% filter(zip_type == "TR") %>% select(file) %>%  
  write_csv("../snakecrate/data/raw_data/tr/s5_tr_filenames.txt", col_names = F)
######


# P 
ibdmdb_uids %>% filter(zip_type == "P") %>% select(file) %>%  
  write_csv("../snakecrate/data/raw_data/p/p_filenames.txt", col_names = F)

ibdmdb_uids_s5  %>% filter(zip_type == "P") %>% select(file) %>%  
  write_csv("../snakecrate/data/raw_data/p/s5_p_filenames.txt", col_names = F)
######


# STD
ibdmdb_uids %>% filter(zip_type == "STD") %>% select(file) %>%  
  write_csv("../snakecrate/data/raw_data/std/std_filenames.txt", col_names = F)

ibdmdb_uids_s5  %>% filter(zip_type == "STD") %>% select(file) %>%  
  write_csv("../snakecrate/data/raw_data/std/s5_std_filenames.txt", col_names = F)
######
```



##Supplementary Chunks
### Chunk s1

```r
#### Here I am frantically looking for the column that its matching on 76C9Y, with external_id HSM6XRR7

# reducing the number of columns to char columns, gets down to ~300 from ~500 
df <- ibdmdb_mvx_only %>% 
  filter(external_id == "HSM6XRR7") %>% 
  select(where(is.character))


df2 <- df[!map_lgl(df, ~ all(is.na(.)))]
df2
```

```
## # A tibble: 1 × 131
##   external_id tube_a_viromics participant_id data_type project     site_sub_coll
##   <chr>       <chr>           <chr>          <chr>     <chr>       <chr>        
## 1 HSM6XRR7    SM-76C9Y        H4019          viromics  H4019C3_MVX H4019C3      
## # … with 125 more variables: research_project <chr>, interval_name <chr>,
## #   site_name <chr>, education_level <chr>, occupation <chr>,
## #   whole_blood_received_at_broad <chr>,
## #   serum_tube_number_1_received_at_csmc <chr>,
## #   serum_tubes_number_2_4_received_at_mgh <chr>, f_lora_received_at_mgh <chr>,
## #   rna_dna_received_at_broad <chr>, ecp_received_at_washington_u <chr>,
## #   diagnosis <chr>, …
```

```r
# df2 only has 131 columns, so I just clicked through them. found:

#`tube_a_viromics`!!!! 
#which matches the exact virome file name: SM-76C9Y


# looking at these "tube" columns is proabably going to be important later, hard listed below 
#names(
  ibdmdb %>% 
  select(contains("tube")) %>% 
  select(where(is.character)) %>% 
  drop_na()
```

```
## # A tibble: 449 × 11
##    tube_a_viromics serum_tube_number_1_received… serum_tubes_num… tube_a_dna_rna
##    <chr>           <chr>                         <chr>            <chr>         
##  1 SM-5QVZE        No                            No               SM-5QVZC      
##  2 SM-5UKLB        No                            No               SM-5UKL9      
##  3 SM-5YSLQ        No                            No               SM-5YSLO      
##  4 SM-61GG5        No                            No               SM-61GG3      
##  5 SM-64NAA        No                            No               SM-64NA8      
##  6 SM-6DQVE        No                            No               SM-6DQVC      
##  7 SM-6IEV2        No                            No               SM-6IEUZ      
##  8 SM-6MYZ3        No                            No               SM-6MYZ1      
##  9 SM-6QU1N        No                            No               SM-6QU1M      
## 10 SM-6UG78        No                            No               SM-6UG77      
## # … with 439 more rows, and 7 more variables: tube_a_metabolomics <chr>,
## #   tube_a_storage <chr>, tube_b_fecal_calprotectin <chr>,
## #   tube_b_proteomics <chr>, stool_sample_id_tube_a_et_oh <chr>,
## #   sample_id_tube_b_no_preservative <chr>,
## #   tube_a_and_b_received_at_broad <chr>
```

```r
#)
```
```
"serum_tube_number_1_received_at_csmc"
"serum_tubes_number_2_4_received_at_mgh"
"tube_a_dna_rna""tube_a_metabolomics"
"tube_a_storage"
"tube_a_viromics"
"tube_b_fecal_calprotectin"
"tube_b_proteomics"
"stool_sample_id_tube_a_et_oh"
"sample_id_tube_b_no_preservative"
"tube_a_and_b_received_at_broad"    
```



### Chunk s2, Download some viromes:

```r
subset5 %>% 
  mutate(
    download_URL = case_when(
      (data_type=="viromics")~"https://ibdmdb.org/tunnel/static/HMP2/Viromics/1732/"
      )
    ) %>% relocate(download_URL)
```





### Chunk s3, check that all in ibdmdb are unique external_ids

```r
#1638
ibdmdb %>% 
  filter(data_type=="metagenomics") %>% 
  distinct(external_id)
```

```
## # A tibble: 1,638 × 1
##    external_id
##    <chr>      
##  1 CSM5FZ3N_P 
##  2 CSM5FZ3R_P 
##  3 CSM5YRY7_P 
##  4 CSM5FZ3V_P 
##  5 CSM5FZ4C_P 
##  6 CSM5MCVD_P 
##  7 CSM5MCVF_P 
##  8 CSM5MCVV_P 
##  9 CSM5MCWI_P 
## 10 CSM5MCXD   
## # … with 1,628 more rows
```

```r
#703
ibdmdb %>% 
  filter(data_type=="viromics") %>% 
  distinct(external_id)
```

```
## # A tibble: 703 × 1
##    external_id
##    <chr>      
##  1 CSM5MCXD   
##  2 CSM67UA2   
##  3 CSM79HGP   
##  4 CSM5MCVN   
##  5 CSM67UBH   
##  6 CSM67UBN   
##  7 CSM79HJW   
##  8 CSM5FZ4M   
##  9 CSM5MCWQ   
## 10 CSM67UBZ   
## # … with 693 more rows
```

```r
#2341
ibdmdb
```

```
## # A tibble: 2,341 × 490
##    external_id tube_a_viromics week_num participant_id data_type    project
##    <chr>       <chr>              <dbl> <chr>          <chr>        <chr>  
##  1 CSM5FZ3N_P  SM-5QVZE               0 C3001          metagenomics G79889 
##  2 CSM5FZ3R_P  SM-5UKLB               2 C3001          metagenomics G79894 
##  3 CSM5YRY7_P  SM-5YSLQ               4 C3001          metagenomics G79903 
##  4 CSM5FZ3V_P  SM-61GG5               6 C3001          metagenomics G79913 
##  5 CSM5FZ4C_P  SM-64NAA               8 C3001          metagenomics G79926 
##  6 CSM5MCVD_P  SM-6DQVE              12 C3001          metagenomics G79969 
##  7 CSM5MCVF_P  SM-6IEV2              14 C3001          metagenomics G80015 
##  8 CSM5MCVV_P  SM-6MYZ3              16 C3001          metagenomics G79979 
##  9 CSM5MCWI_P  SM-6QU1N              18 C3001          metagenomics G80050 
## 10 CSM5MCXD    SM-6UG78              20 C3001          metagenomics G110156
## # … with 2,331 more rows, and 484 more variables: site_sub_coll <chr>,
## #   date_of_receipt <date>, interval_days <dbl>, visit_num <dbl>,
## #   research_project <chr>, pdo_number <chr>, gssr_i_ds <dbl>, product <chr>,
## #   lcset <chr>, aggregated_lanes <chr>, wr_id <dbl>,
## #   number_lanes_in_aggregation <dbl>, reads_raw <dbl>, reads_filtered <dbl>,
## #   reads_qc_fail <dbl>, reads_human <dbl>, reads_ribosomal <dbl>,
## #   reads_viral <dbl>, delta <lgl>, interval_name <chr>, …
```

```r
#1656
ibdmdb %>% 
  distinct(external_id)
```

```
## # A tibble: 1,656 × 1
##    external_id
##    <chr>      
##  1 CSM5FZ3N_P 
##  2 CSM5FZ3R_P 
##  3 CSM5YRY7_P 
##  4 CSM5FZ3V_P 
##  5 CSM5FZ4C_P 
##  6 CSM5MCVD_P 
##  7 CSM5MCVF_P 
##  8 CSM5MCVV_P 
##  9 CSM5MCWI_P 
## 10 CSM5MCXD   
## # … with 1,646 more rows
```

```r
1638-703
```

```
## [1] 935
```

```r
#935

2341-935
```

```
## [1] 1406
```

```r
#1406
```

### bash sandbox


```bash
while IFS="," read -r file prefix zip_type uid
do
  echo "file is $file"
  echo "file prefix is $prefix"
  echo "zip_type is $zip_type"
  echo "uid is $uid"
done < <(tail -n +2 csvs_and_other_metadata/file_prefix_zip_uid.csv)
```

```
## file is CSM5FZ3N_P.fastq.gz
## file prefix is CSM5FZ3N_P
## zip_type is P
## uid is CSM5FZ3N_P
## file is CSM5FZ3R_P.fastq.gz
## file prefix is CSM5FZ3R_P
## zip_type is P
## uid is CSM5FZ3R_P
## file is CSM5YRY7_P.fastq.gz
## file prefix is CSM5YRY7_P
## zip_type is P
## uid is CSM5YRY7_P
## file is CSM5FZ3V_P.fastq.gz
## file prefix is CSM5FZ3V_P
## zip_type is P
## uid is CSM5FZ3V_P
## file is CSM5FZ4C_P.fastq.gz
## file prefix is CSM5FZ4C_P
## zip_type is P
## uid is CSM5FZ4C_P
## file is CSM5MCVD_P.fastq.gz
## file prefix is CSM5MCVD_P
## zip_type is P
## uid is CSM5MCVD_P
## file is CSM5MCVF_P.fastq.gz
## file prefix is CSM5MCVF_P
## zip_type is P
## uid is CSM5MCVF_P
## file is CSM5MCVV_P.fastq.gz
## file prefix is CSM5MCVV_P
## zip_type is P
## uid is CSM5MCVV_P
## file is CSM5MCWI_P.fastq.gz
## file prefix is CSM5MCWI_P
## zip_type is P
## uid is CSM5MCWI_P
## file is CSM5MCXD.tar
## file prefix is CSM5MCXD
## zip_type is STD
## uid is CSM5MCXD
## file is CSM5MCYS.tar
## file prefix is CSM5MCYS
## zip_type is STD
## uid is CSM5MCYS
## file is CSM67U9J.tar
## file prefix is CSM67U9J
## zip_type is STD
## uid is CSM67U9J
## file is CSM67UA2.tar
## file prefix is CSM67UA2
## zip_type is STD
## uid is CSM67UA2
## file is CSM67UGC.tar
## file prefix is CSM67UGC
## zip_type is STD
## uid is CSM67UGC
## file is CSM79HG5.tar
## file prefix is CSM79HG5
## zip_type is STD
## uid is CSM79HG5
## file is CSM79HGP.tar
## file prefix is CSM79HGP
## zip_type is STD
## uid is CSM79HGP
## file is SM-6UG79.tar
## file prefix is SM-6UG79
## zip_type is VIR
## uid is SM-6UG79
## file is SM-71WY3.tar
## file prefix is SM-71WY3
## zip_type is VIR
## uid is SM-71WY3
## file is SM-7CRX4.tar
## file prefix is SM-7CRX4
## zip_type is VIR
## uid is SM-7CRX4
## file is CSM5FZ3T_P.fastq.gz
## file prefix is CSM5FZ3T_P
## zip_type is P
## uid is CSM5FZ3T_P
## file is CSM5FZ3X_P.fastq.gz
## file prefix is CSM5FZ3X_P
## zip_type is P
## uid is CSM5FZ3X_P
## file is CSM5FZ3Z_P.fastq.gz
## file prefix is CSM5FZ3Z_P
## zip_type is P
## uid is CSM5FZ3Z_P
## file is CSM5FZ42_P.fastq.gz
## file prefix is CSM5FZ42_P
## zip_type is P
## uid is CSM5FZ42_P
## file is CSM5FZ44_P.fastq.gz
## file prefix is CSM5FZ44_P
## zip_type is P
## uid is CSM5FZ44_P
## file is CSM5FZ46_P.fastq.gz
## file prefix is CSM5FZ46_P
## zip_type is P
## uid is CSM5FZ46_P
## file is CSM5MCVJ_P.fastq.gz
## file prefix is CSM5MCVJ_P
## zip_type is P
## uid is CSM5MCVJ_P
## file is CSM5MCVL.tar
## file prefix is CSM5MCVL
## zip_type is STD
## uid is CSM5MCVL
## file is CSM5MCVN.tar
## file prefix is CSM5MCVN
## zip_type is STD
## uid is CSM5MCVN
## file is CSM67UBF.tar
## file prefix is CSM67UBF
## zip_type is STD
## uid is CSM67UBF
## file is CSM67UBH.tar
## file prefix is CSM67UBH
## zip_type is STD
## uid is CSM67UBH
## file is CSM67UBN.tar
## file prefix is CSM67UBN
## zip_type is STD
## uid is CSM67UBN
## file is CSM67UBR.tar
## file prefix is CSM67UBR
## zip_type is STD
## uid is CSM67UBR
## file is CSM79HJW.tar
## file prefix is CSM79HJW
## zip_type is STD
## uid is CSM79HJW
## file is CSM79HJY.tar
## file prefix is CSM79HJY
## zip_type is STD
## uid is CSM79HJY
## file is SM-6X9WV.tar
## file prefix is SM-6X9WV
## zip_type is VIR
## uid is SM-6X9WV
## file is SM-76CAU.tar
## file prefix is SM-76CAU
## zip_type is VIR
## uid is SM-76CAU
## file is SM-7CP3H.tar
## file prefix is SM-7CP3H
## zip_type is VIR
## uid is SM-7CP3H
## file is SM-7EWTT.tar
## file prefix is SM-7EWTT
## zip_type is VIR
## uid is SM-7EWTT
## file is CSM5FZ4E_P.fastq.gz
## file prefix is CSM5FZ4E_P
## zip_type is P
## uid is CSM5FZ4E_P
## file is CSM5FZ4G_P.fastq.gz
## file prefix is CSM5FZ4G_P
## zip_type is P
## uid is CSM5FZ4G_P
## file is CSM5FZ4K_P.fastq.gz
## file prefix is CSM5FZ4K_P
## zip_type is P
## uid is CSM5FZ4K_P
## file is CSM5FZ4M.tar
## file prefix is CSM5FZ4M
## zip_type is STD
## uid is CSM5FZ4M
## file is CSM5MCWM_P.fastq.gz
## file prefix is CSM5MCWM_P
## zip_type is P
## uid is CSM5MCWM_P
## file is CSM5MCWQ.tar
## file prefix is CSM5MCWQ
## zip_type is STD
## uid is CSM5MCWQ
## file is CSM67UBX.tar
## file prefix is CSM67UBX
## zip_type is STD
## uid is CSM67UBX
## file is CSM67UBZ.tar
## file prefix is CSM67UBZ
## zip_type is STD
## uid is CSM67UBZ
## file is CSM67UC6.tar
## file prefix is CSM67UC6
## zip_type is STD
## uid is CSM67UC6
## file is CSM79HLM.tar
## file prefix is CSM79HLM
## zip_type is STD
## uid is CSM79HLM
## file is SM-6OSOR.tar
## file prefix is SM-6OSOR
## zip_type is VIR
## uid is SM-6OSOR
## file is SM-6Y2V3.tar
## file prefix is SM-6Y2V3
## zip_type is VIR
## uid is SM-6Y2V3
## file is SM-77VQ4.tar
## file prefix is SM-77VQ4
## zip_type is VIR
## uid is SM-77VQ4
## file is SM-7EWUL.tar
## file prefix is SM-7EWUL
## zip_type is VIR
## uid is SM-7EWUL
## file is CSM5FZ4A_P.fastq.gz
## file prefix is CSM5FZ4A_P
## zip_type is P
## uid is CSM5FZ4A_P
## file is CSM5MCU8_P.fastq.gz
## file prefix is CSM5MCU8_P
## zip_type is P
## uid is CSM5MCU8_P
## file is CSM5MCUA_P.fastq.gz
## file prefix is CSM5MCUA_P
## zip_type is P
## uid is CSM5MCUA_P
## file is CSM5MCUC_P.fastq.gz
## file prefix is CSM5MCUC_P
## zip_type is P
## uid is CSM5MCUC_P
## file is CSM5MCUE_P.fastq.gz
## file prefix is CSM5MCUE_P
## zip_type is P
## uid is CSM5MCUE_P
## file is CSM5MCXH.tar
## file prefix is CSM5MCXH
## zip_type is STD
## uid is CSM5MCXH
## file is CSM5MCXJ.tar
## file prefix is CSM5MCXJ
## zip_type is STD
## uid is CSM5MCXJ
## file is CSM5MCXL.tar
## file prefix is CSM5MCXL
## zip_type is STD
## uid is CSM5MCXL
## file is CSM5MCXN.tar
## file prefix is CSM5MCXN
## zip_type is STD
## uid is CSM5MCXN
## file is CSM5MCXP.tar
## file prefix is CSM5MCXP
## zip_type is STD
## uid is CSM5MCXP
## file is CSM5MCXR.tar
## file prefix is CSM5MCXR
## zip_type is STD
## uid is CSM5MCXR
## file is CSM67UDF.tar
## file prefix is CSM67UDF
## zip_type is STD
## uid is CSM67UDF
## file is CSM67UDJ.tar
## file prefix is CSM67UDJ
## zip_type is STD
## uid is CSM67UDJ
## file is CSM67UDN.tar
## file prefix is CSM67UDN
## zip_type is STD
## uid is CSM67UDN
## file is CSM67UDR_TR.tar
## file prefix is CSM67UDR_TR
## zip_type is TR
## uid is CSM67UDR_TR
## file is CSM67UDR.tar
## file prefix is CSM67UDR
## zip_type is STD
## uid is CSM67UDR
## file is CSM67UDY.tar
## file prefix is CSM67UDY
## zip_type is STD
## uid is CSM67UDY
## file is CSM79HLA_TR.tar
## file prefix is CSM79HLA_TR
## zip_type is TR
## uid is CSM79HLA_TR
## file is CSM79HLA.tar
## file prefix is CSM79HLA
## zip_type is STD
## uid is CSM79HLA
## file is CSM79HLG.tar
## file prefix is CSM79HLG
## zip_type is STD
## uid is CSM79HLG
## file is CSM79HLE.tar
## file prefix is CSM79HLE
## zip_type is STD
## uid is CSM79HLE
## file is CSM79HLC.tar
## file prefix is CSM79HLC
## zip_type is STD
## uid is CSM79HLC
## file is CSM79HLI.tar
## file prefix is CSM79HLI
## zip_type is STD
## uid is CSM79HLI
## file is CSM79HLK.tar
## file prefix is CSM79HLK
## zip_type is STD
## uid is CSM79HLK
## file is SM-6WJN6.tar
## file prefix is SM-6WJN6
## zip_type is VIR
## uid is SM-6WJN6
## file is SM-6X9X4.tar
## file prefix is SM-6X9X4
## zip_type is VIR
## uid is SM-6X9X4
## file is SM-6YAZO.tar
## file prefix is SM-6YAZO
## zip_type is VIR
## uid is SM-6YAZO
## file is SM-6ZKSM.tar
## file prefix is SM-6ZKSM
## zip_type is VIR
## uid is SM-6ZKSM
## file is SM-71WXM.tar
## file prefix is SM-71WXM
## zip_type is VIR
## uid is SM-71WXM
## file is SM-73JY4.tar
## file prefix is SM-73JY4
## zip_type is VIR
## uid is SM-73JY4
## file is SM-76EOJ.tar
## file prefix is SM-76EOJ
## zip_type is VIR
## uid is SM-76EOJ
## file is SM-791BX.tar
## file prefix is SM-791BX
## zip_type is VIR
## uid is SM-791BX
## file is SM-7BF2J.tar
## file prefix is SM-7BF2J
## zip_type is VIR
## uid is SM-7BF2J
## file is SM-7CRJ8.tar
## file prefix is SM-7CRJ8
## zip_type is VIR
## uid is SM-7CRJ8
## file is SM-7EWUH.tar
## file prefix is SM-7EWUH
## zip_type is VIR
## uid is SM-7EWUH
## file is SM-7GYJO.tar
## file prefix is SM-7GYJO
## zip_type is VIR
## uid is SM-7GYJO
## file is SM-7IL1E.tar
## file prefix is SM-7IL1E
## zip_type is VIR
## uid is SM-7IL1E
## file is SM-7KPUR.tar
## file prefix is SM-7KPUR
## zip_type is VIR
## uid is SM-7KPUR
## file is SM-7MOQJ.tar
## file prefix is SM-7MOQJ
## zip_type is VIR
## uid is SM-7MOQJ
## file is CSM5MCUQ_P.fastq.gz
## file prefix is CSM5MCUQ_P
## zip_type is P
## uid is CSM5MCUQ_P
## file is CSM5MCUS_P.fastq.gz
## file prefix is CSM5MCUS_P
## zip_type is P
## uid is CSM5MCUS_P
## file is CSM5MCUW_P.fastq.gz
## file prefix is CSM5MCUW_P
## zip_type is P
## uid is CSM5MCUW_P
## file is CSM5MCUY_P.fastq.gz
## file prefix is CSM5MCUY_P
## zip_type is P
## uid is CSM5MCUY_P
## file is CSM5MCY4.tar
## file prefix is CSM5MCY4
## zip_type is STD
## uid is CSM5MCY4
## file is CSM5MCY8.tar
## file prefix is CSM5MCY8
## zip_type is STD
## uid is CSM5MCY8
## file is CSM67UE3.tar
## file prefix is CSM67UE3
## zip_type is STD
## uid is CSM67UE3
## file is CSM67UE7.tar
## file prefix is CSM67UE7
## zip_type is STD
## uid is CSM67UE7
## file is CSM67UEA.tar
## file prefix is CSM67UEA
## zip_type is STD
## uid is CSM67UEA
## file is CSM67UEM.tar
## file prefix is CSM67UEM
## zip_type is STD
## uid is CSM67UEM
## file is CSM67UEI.tar
## file prefix is CSM67UEI
## zip_type is STD
## uid is CSM67UEI
## file is CSM79HO1.tar
## file prefix is CSM79HO1
## zip_type is STD
## uid is CSM79HO1
## file is SM-6WOCE.tar
## file prefix is SM-6WOCE
## zip_type is VIR
## uid is SM-6WOCE
## file is SM-6ZEVI.tar
## file prefix is SM-6ZEVI
## zip_type is VIR
## uid is SM-6ZEVI
## file is SM-7AA2E.tar
## file prefix is SM-7AA2E
## zip_type is VIR
## uid is SM-7AA2E
## file is SM-7BP5L.tar
## file prefix is SM-7BP5L
## zip_type is VIR
## uid is SM-7BP5L
## file is SM-7I1G8.tar
## file prefix is SM-7I1G8
## zip_type is VIR
## uid is SM-7I1G8
## file is SM-7K25Z.tar
## file prefix is SM-7K25Z
## zip_type is VIR
## uid is SM-7K25Z
## file is CSM5MCTZ_P.fastq.gz
## file prefix is CSM5MCTZ_P
## zip_type is P
## uid is CSM5MCTZ_P
## file is CSM5MCUG_P.fastq.gz
## file prefix is CSM5MCUG_P
## zip_type is P
## uid is CSM5MCUG_P
## file is CSM5MCUK_P.fastq.gz
## file prefix is CSM5MCUK_P
## zip_type is P
## uid is CSM5MCUK_P
## file is CSM5MCUO.tar
## file prefix is CSM5MCUO
## zip_type is STD
## uid is CSM5MCUO
## file is CSM5MCX3.tar
## file prefix is CSM5MCX3
## zip_type is STD
## uid is CSM5MCX3
## file is CSM67UFV.tar
## file prefix is CSM67UFV
## zip_type is STD
## uid is CSM67UFV
## file is CSM67UFZ.tar
## file prefix is CSM67UFZ
## zip_type is STD
## uid is CSM67UFZ
## file is CSM67UG8.tar
## file prefix is CSM67UG8
## zip_type is STD
## uid is CSM67UG8
## file is CSM79HMN.tar
## file prefix is CSM79HMN
## zip_type is STD
## uid is CSM79HMN
## file is CSM79HMP.tar
## file prefix is CSM79HMP
## zip_type is STD
## uid is CSM79HMP
## file is CSM79HMT.tar
## file prefix is CSM79HMT
## zip_type is STD
## uid is CSM79HMT
## file is SM-6V2CZ.tar
## file prefix is SM-6V2CZ
## zip_type is VIR
## uid is SM-6V2CZ
## file is SM-6ZKSZ.tar
## file prefix is SM-6ZKSZ
## zip_type is VIR
## uid is SM-6ZKSZ
## file is SM-7BF2O.tar
## file prefix is SM-7BF2O
## zip_type is VIR
## uid is SM-7BF2O
## file is SM-7ORRO.tar
## file prefix is SM-7ORRO
## zip_type is VIR
## uid is SM-7ORRO
## file is CSM5MCVB_P.fastq.gz
## file prefix is CSM5MCVB_P
## zip_type is P
## uid is CSM5MCVB_P
## file is CSM5MCV1_P.fastq.gz
## file prefix is CSM5MCV1_P
## zip_type is P
## uid is CSM5MCV1_P
## file is CSM5MCV5_P.fastq.gz
## file prefix is CSM5MCV5_P
## zip_type is P
## uid is CSM5MCV5_P
## file is CSM5MCU4_P.fastq.gz
## file prefix is CSM5MCU4_P
## zip_type is P
## uid is CSM5MCU4_P
## file is CSM5MCVZ_P.fastq.gz
## file prefix is CSM5MCVZ_P
## zip_type is P
## uid is CSM5MCVZ_P
## file is CSM5MCW4_P.fastq.gz
## file prefix is CSM5MCW4_P
## zip_type is P
## uid is CSM5MCW4_P
## file is CSM5MCW6.tar
## file prefix is CSM5MCW6
## zip_type is STD
## uid is CSM5MCW6
## file is CSM5MCYW.tar
## file prefix is CSM5MCYW
## zip_type is STD
## uid is CSM5MCYW
## file is CSM5MCZ3.tar
## file prefix is CSM5MCZ3
## zip_type is STD
## uid is CSM5MCZ3
## file is CSM5MCZ5.tar
## file prefix is CSM5MCZ5
## zip_type is STD
## uid is CSM5MCZ5
## file is CSM5MCZ7.tar
## file prefix is CSM5MCZ7
## zip_type is STD
## uid is CSM5MCZ7
## file is CSM79HKT.tar
## file prefix is CSM79HKT
## zip_type is STD
## uid is CSM79HKT
## file is CSM79HKV.tar
## file prefix is CSM79HKV
## zip_type is STD
## uid is CSM79HKV
## file is CSM79HQ9.tar
## file prefix is CSM79HQ9
## zip_type is STD
## uid is CSM79HQ9
## file is CSM79HQB.tar
## file prefix is CSM79HQB
## zip_type is STD
## uid is CSM79HQB
## file is CSM79HQF.tar
## file prefix is CSM79HQF
## zip_type is STD
## uid is CSM79HQF
## file is SM-6XSVD.tar
## file prefix is SM-6XSVD
## zip_type is VIR
## uid is SM-6XSVD
## file is SM-74Y95.tar
## file prefix is SM-74Y95
## zip_type is VIR
## uid is SM-74Y95
## file is SM-7K267.tar
## file prefix is SM-7K267
## zip_type is VIR
## uid is SM-7K267
## file is SM-9LB6K.tar
## file prefix is SM-9LB6K
## zip_type is VIR
## uid is SM-9LB6K
## file is CSM5MCWA_P.fastq.gz
## file prefix is CSM5MCWA_P
## zip_type is P
## uid is CSM5MCWA_P
## file is CSM5MCWC.tar
## file prefix is CSM5MCWC
## zip_type is STD
## uid is CSM5MCWC
## file is CSM5MCWE.tar
## file prefix is CSM5MCWE
## zip_type is STD
## uid is CSM5MCWE
## file is CSM5MCWG.tar
## file prefix is CSM5MCWG
## zip_type is STD
## uid is CSM5MCWG
## file is CSM67UAA.tar
## file prefix is CSM67UAA
## zip_type is STD
## uid is CSM67UAA
## file is CSM67UAG.tar
## file prefix is CSM67UAG
## zip_type is STD
## uid is CSM67UAG
## file is CSM79HI3.tar
## file prefix is CSM79HI3
## zip_type is STD
## uid is CSM79HI3
## file is CSM79HI7.tar
## file prefix is CSM79HI7
## zip_type is STD
## uid is CSM79HI7
## file is CSM79HIB.tar
## file prefix is CSM79HIB
## zip_type is STD
## uid is CSM79HIB
## file is CSM7KOJE.tar
## file prefix is CSM7KOJE
## zip_type is STD
## uid is CSM7KOJE
## file is CSM7KOJG.tar
## file prefix is CSM7KOJG
## zip_type is STD
## uid is CSM7KOJG
## file is CSM7KOJO.tar
## file prefix is CSM7KOJO
## zip_type is STD
## uid is CSM7KOJO
## file is SM-6X9WZ.tar
## file prefix is SM-6X9WZ
## zip_type is VIR
## uid is SM-6X9WZ
## file is SM-6YAZS.tar
## file prefix is SM-6YAZS
## zip_type is VIR
## uid is SM-6YAZS
## file is SM-76CAM.tar
## file prefix is SM-76CAM
## zip_type is VIR
## uid is SM-76CAM
## file is SM-7EWUE.tar
## file prefix is SM-7EWUE
## zip_type is VIR
## uid is SM-7EWUE
## file is SM-7IL1P.tar
## file prefix is SM-7IL1P
## zip_type is VIR
## uid is SM-7IL1P
## file is SM-7MCUD.tar
## file prefix is SM-7MCUD
## zip_type is VIR
## uid is SM-7MCUD
## file is SM-9HC8A.tar
## file prefix is SM-9HC8A
## zip_type is VIR
## uid is SM-9HC8A
## file is SM-9RRZO.tar
## file prefix is SM-9RRZO
## zip_type is VIR
## uid is SM-9RRZO
## file is CSM5MCWK_P.fastq.gz
## file prefix is CSM5MCWK_P
## zip_type is P
## uid is CSM5MCWK_P
## file is CSM5MCXT.tar
## file prefix is CSM5MCXT
## zip_type is STD
## uid is CSM5MCXT
## file is CSM5MCXV.tar
## file prefix is CSM5MCXV
## zip_type is STD
## uid is CSM5MCXV
## file is CSM5MCXX_P.fastq.gz
## file prefix is CSM5MCXX_P
## zip_type is P
## uid is CSM5MCXX_P
## file is CSM5MCXZ_P.fastq.gz
## file prefix is CSM5MCXZ_P
## zip_type is P
## uid is CSM5MCXZ_P
## file is CSM5MCY2.tar
## file prefix is CSM5MCY2
## zip_type is STD
## uid is CSM5MCY2
## file is CSM67UCK.tar
## file prefix is CSM67UCK
## zip_type is STD
## uid is CSM67UCK
## file is CSM79HK9.tar
## file prefix is CSM79HK9
## zip_type is STD
## uid is CSM79HK9
## file is CSM79HKB.tar
## file prefix is CSM79HKB
## zip_type is STD
## uid is CSM79HKB
## file is CSM79HOF.tar
## file prefix is CSM79HOF
## zip_type is STD
## uid is CSM79HOF
## file is CSM79HOH.tar
## file prefix is CSM79HOH
## zip_type is STD
## uid is CSM79HOH
## file is CSM7KOL4.tar
## file prefix is CSM7KOL4
## zip_type is STD
## uid is CSM7KOL4
## file is CSM7KOLA.tar
## file prefix is CSM7KOLA
## zip_type is STD
## uid is CSM7KOLA
## file is CSM7KOLE.tar
## file prefix is CSM7KOLE
## zip_type is STD
## uid is CSM7KOLE
## file is SM-6WU1R.tar
## file prefix is SM-6WU1R
## zip_type is VIR
## uid is SM-6WU1R
## file is SM-72PHI.tar
## file prefix is SM-72PHI
## zip_type is VIR
## uid is SM-72PHI
## file is SM-7IL1U.tar
## file prefix is SM-7IL1U
## zip_type is VIR
## uid is SM-7IL1U
## file is SM-9IT18.tar
## file prefix is SM-9IT18
## zip_type is VIR
## uid is SM-9IT18
## file is SM-9QMOQ.tar
## file prefix is SM-9QMOQ
## zip_type is VIR
## uid is SM-9QMOQ
## file is SM-9VWCS.tar
## file prefix is SM-9VWCS
## zip_type is VIR
## uid is SM-9VWCS
## file is CSM5MCXB_P.fastq.gz
## file prefix is CSM5MCXB_P
## zip_type is P
## uid is CSM5MCXB_P
## file is CSM5MCYI_P.fastq.gz
## file prefix is CSM5MCYI_P
## zip_type is P
## uid is CSM5MCYI_P
## file is CSM5MCYM_P.fastq.gz
## file prefix is CSM5MCYM_P
## zip_type is P
## uid is CSM5MCYM_P
## file is CSM5MCYO_P.fastq.gz
## file prefix is CSM5MCYO_P
## zip_type is P
## uid is CSM5MCYO_P
## file is CSM5MCYQ_P.fastq.gz
## file prefix is CSM5MCYQ_P
## zip_type is P
## uid is CSM5MCYQ_P
## file is CSM67UEP_P.fastq.gz
## file prefix is CSM67UEP_P
## zip_type is P
## uid is CSM67UEP_P
## file is CSM67UET_P.fastq.gz
## file prefix is CSM67UET_P
## zip_type is P
## uid is CSM67UET_P
## file is CSM67UEW_TR.tar
## file prefix is CSM67UEW_TR
## zip_type is TR
## uid is CSM67UEW_TR
## file is CSM67UEW_P.fastq.gz
## file prefix is CSM67UEW_P
## zip_type is P
## uid is CSM67UEW_P
## file is CSM67UEW.tar
## file prefix is CSM67UEW
## zip_type is STD
## uid is CSM67UEW
## file is CSM67UF1_P.fastq.gz
## file prefix is CSM67UF1_P
## zip_type is P
## uid is CSM67UF1_P
## file is CSM67UF1.tar
## file prefix is CSM67UF1
## zip_type is STD
## uid is CSM67UF1
## file is CSM67UF5.tar
## file prefix is CSM67UF5
## zip_type is STD
## uid is CSM67UF5
## file is CSM79HM1.tar
## file prefix is CSM79HM1
## zip_type is STD
## uid is CSM79HM1
## file is CSM79HQT_P.fastq.gz
## file prefix is CSM79HQT_P
## zip_type is P
## uid is CSM79HQT_P
## file is CSM79HM5_P.fastq.gz
## file prefix is CSM79HM5_P
## zip_type is P
## uid is CSM79HM5_P
## file is CSM79HM7.tar
## file prefix is CSM79HM7
## zip_type is STD
## uid is CSM79HM7
## file is CSM79HM9_P.fastq.gz
## file prefix is CSM79HM9_P
## zip_type is P
## uid is CSM79HM9_P
## file is CSM7KOMP.tar
## file prefix is CSM7KOMP
## zip_type is STD
## uid is CSM7KOMP
## file is CSM7KOMR_P.fastq.gz
## file prefix is CSM7KOMR_P
## zip_type is P
## uid is CSM7KOMR_P
## file is CSM7KOMT.tar
## file prefix is CSM7KOMT
## zip_type is STD
## uid is CSM7KOMT
## file is CSM7KOMV_P.fastq.gz
## file prefix is CSM7KOMV_P
## zip_type is P
## uid is CSM7KOMV_P
## file is SM-7CM7N.tar
## file prefix is SM-7CM7N
## zip_type is VIR
## uid is SM-7CM7N
## file is SM-7FVSR.tar
## file prefix is SM-7FVSR
## zip_type is VIR
## uid is SM-7FVSR
## file is SM-9PPV2.tar
## file prefix is SM-9PPV2
## zip_type is VIR
## uid is SM-9PPV2
## file is SM-9VE4Z.tar
## file prefix is SM-9VE4Z
## zip_type is VIR
## uid is SM-9VE4Z
## file is CSM5MCXF_P.fastq.gz
## file prefix is CSM5MCXF_P
## zip_type is P
## uid is CSM5MCXF_P
## file is CSM5MCZB.tar
## file prefix is CSM5MCZB
## zip_type is STD
## uid is CSM5MCZB
## file is CSM5MCZD.tar
## file prefix is CSM5MCZD
## zip_type is STD
## uid is CSM5MCZD
## file is CSM5MCZF.tar
## file prefix is CSM5MCZF
## zip_type is STD
## uid is CSM5MCZF
## file is CSM67U9B.tar
## file prefix is CSM67U9B
## zip_type is STD
## uid is CSM67U9B
## file is CSM67U9D.tar
## file prefix is CSM67U9D
## zip_type is STD
## uid is CSM67U9D
## file is CSM67UGO.tar
## file prefix is CSM67UGO
## zip_type is STD
## uid is CSM67UGO
## file is CSM79HOJ.tar
## file prefix is CSM79HOJ
## zip_type is STD
## uid is CSM79HOJ
## file is CSM79HOL.tar
## file prefix is CSM79HOL
## zip_type is STD
## uid is CSM79HOL
## file is CSM79HOT.tar
## file prefix is CSM79HOT
## zip_type is STD
## uid is CSM79HOT
## file is CSM7KOMX.tar
## file prefix is CSM7KOMX
## zip_type is STD
## uid is CSM7KOMX
## file is CSM7KOMZ.tar
## file prefix is CSM7KOMZ
## zip_type is STD
## uid is CSM7KOMZ
## file is CSM7KON2.tar
## file prefix is CSM7KON2
## zip_type is STD
## uid is CSM7KON2
## file is CSM7KON8.tar
## file prefix is CSM7KON8
## zip_type is STD
## uid is CSM7KON8
## file is SM-6ZKSV.tar
## file prefix is SM-6ZKSV
## zip_type is VIR
## uid is SM-6ZKSV
## file is SM-6ZT43.tar
## file prefix is SM-6ZT43
## zip_type is VIR
## uid is SM-6ZT43
## file is SM-72PHU.tar
## file prefix is SM-72PHU
## zip_type is VIR
## uid is SM-72PHU
## file is SM-77FXK.tar
## file prefix is SM-77FXK
## zip_type is VIR
## uid is SM-77FXK
## file is SM-7DN3W.tar
## file prefix is SM-7DN3W
## zip_type is VIR
## uid is SM-7DN3W
## file is SM-7LDIT.tar
## file prefix is SM-7LDIT
## zip_type is VIR
## uid is SM-7LDIT
## file is SM-9SIJ3.tar
## file prefix is SM-9SIJ3
## zip_type is VIR
## uid is SM-9SIJ3
## file is SM-9ZA5W.tar
## file prefix is SM-9ZA5W
## zip_type is VIR
## uid is SM-9ZA5W
## file is CSM5MCYU_P.fastq.gz
## file prefix is CSM5MCYU_P
## zip_type is P
## uid is CSM5MCYU_P
## file is CSM67U9H_P.fastq.gz
## file prefix is CSM67U9H_P
## zip_type is P
## uid is CSM67U9H_P
## file is CSM67U9H.tar
## file prefix is CSM67U9H
## zip_type is STD
## uid is CSM67U9H
## file is CSM67U9N.tar
## file prefix is CSM67U9N
## zip_type is STD
## uid is CSM67U9N
## file is CSM67U9P_P.fastq.gz
## file prefix is CSM67U9P_P
## zip_type is P
## uid is CSM67U9P_P
## file is CSM67U9R_P.fastq.gz
## file prefix is CSM67U9R_P
## zip_type is P
## uid is CSM67U9R_P
## file is CSM79HGD_P.fastq.gz
## file prefix is CSM79HGD_P
## zip_type is P
## uid is CSM79HGD_P
## file is CSM79HGF_P.fastq.gz
## file prefix is CSM79HGF_P
## zip_type is P
## uid is CSM79HGF_P
## file is CSM79HGF.tar
## file prefix is CSM79HGF
## zip_type is STD
## uid is CSM79HGF
## file is CSM79HGH_P.fastq.gz
## file prefix is CSM79HGH_P
## zip_type is P
## uid is CSM79HGH_P
## file is CSM79HGJ_P.fastq.gz
## file prefix is CSM79HGJ_P
## zip_type is P
## uid is CSM79HGJ_P
## file is CSM79HGL_P.fastq.gz
## file prefix is CSM79HGL_P
## zip_type is P
## uid is CSM79HGL_P
## file is CSM79HGN_P.fastq.gz
## file prefix is CSM79HGN_P
## zip_type is P
## uid is CSM79HGN_P
## file is CSM79HPK.tar
## file prefix is CSM79HPK
## zip_type is STD
## uid is CSM79HPK
## file is CSM79HPM_P.fastq.gz
## file prefix is CSM79HPM_P
## zip_type is P
## uid is CSM79HPM_P
## file is CSM79HPS.tar
## file prefix is CSM79HPS
## zip_type is STD
## uid is CSM79HPS
## file is CSM79HPQ_P.fastq.gz
## file prefix is CSM79HPQ_P
## zip_type is P
## uid is CSM79HPQ_P
## file is CSM79HPO.tar
## file prefix is CSM79HPO
## zip_type is STD
## uid is CSM79HPO
## file is CSM79HPU.tar
## file prefix is CSM79HPU
## zip_type is STD
## uid is CSM79HPU
## file is CSM7KONS_P.fastq.gz
## file prefix is CSM7KONS_P
## zip_type is P
## uid is CSM7KONS_P
## file is CSM7KONU.tar
## file prefix is CSM7KONU
## zip_type is STD
## uid is CSM7KONU
## file is CSM7KONW_P.fastq.gz
## file prefix is CSM7KONW_P
## zip_type is P
## uid is CSM7KONW_P
## file is SM-6ZUHG.tar
## file prefix is SM-6ZUHG
## zip_type is VIR
## uid is SM-6ZUHG
## file is SM-7DN3O.tar
## file prefix is SM-7DN3O
## zip_type is VIR
## uid is SM-7DN3O
## file is SM-7LDJ6.tar
## file prefix is SM-7LDJ6
## zip_type is VIR
## uid is SM-7LDJ6
## file is SM-9JY9Z.tar
## file prefix is SM-9JY9Z
## zip_type is VIR
## uid is SM-9JY9Z
## file is SM-9SNKS.tar
## file prefix is SM-9SNKS
## zip_type is VIR
## uid is SM-9SNKS
## file is CSM67U9T_P.fastq.gz
## file prefix is CSM67U9T_P
## zip_type is P
## uid is CSM67U9T_P
## file is CSM67UAK.tar
## file prefix is CSM67UAK
## zip_type is STD
## uid is CSM67UAK
## file is CSM67UAM.tar
## file prefix is CSM67UAM
## zip_type is STD
## uid is CSM67UAM
## file is CSM67UAO.tar
## file prefix is CSM67UAO
## zip_type is STD
## uid is CSM67UAO
## file is CSM67UAQ.tar
## file prefix is CSM67UAQ
## zip_type is STD
## uid is CSM67UAQ
## file is CSM67UAS.tar
## file prefix is CSM67UAS
## zip_type is STD
## uid is CSM67UAS
## file is CSM79HID.tar
## file prefix is CSM79HID
## zip_type is STD
## uid is CSM79HID
## file is CSM79HIF.tar
## file prefix is CSM79HIF
## zip_type is STD
## uid is CSM79HIF
## file is CSM79HIH.tar
## file prefix is CSM79HIH
## zip_type is STD
## uid is CSM79HIH
## file is CSM79HIJ.tar
## file prefix is CSM79HIJ
## zip_type is STD
## uid is CSM79HIJ
## file is CSM79HIL.tar
## file prefix is CSM79HIL
## zip_type is STD
## uid is CSM79HIL
## file is CSM79HIN.tar
## file prefix is CSM79HIN
## zip_type is STD
## uid is CSM79HIN
## file is CSM7KOJQ.tar
## file prefix is CSM7KOJQ
## zip_type is STD
## uid is CSM7KOJQ
## file is CSM7KOJS.tar
## file prefix is CSM7KOJS
## zip_type is STD
## uid is CSM7KOJS
## file is CSM7KOJU.tar
## file prefix is CSM7KOJU
## zip_type is STD
## uid is CSM7KOJU
## file is CSM7KOJW.tar
## file prefix is CSM7KOJW
## zip_type is STD
## uid is CSM7KOJW
## file is CSM7KOJY.tar
## file prefix is CSM7KOJY
## zip_type is STD
## uid is CSM7KOJY
## file is CSM7KOK1.tar
## file prefix is CSM7KOK1
## zip_type is STD
## uid is CSM7KOK1
## file is CSM7KOPG.tar
## file prefix is CSM7KOPG
## zip_type is STD
## uid is CSM7KOPG
## file is CSM7KOPI.tar
## file prefix is CSM7KOPI
## zip_type is STD
## uid is CSM7KOPI
## file is CSM7KOPK.tar
## file prefix is CSM7KOPK
## zip_type is STD
## uid is CSM7KOPK
## file is CSM7KOPM.tar
## file prefix is CSM7KOPM
## zip_type is STD
## uid is CSM7KOPM
## file is CSM7KOPO.tar
## file prefix is CSM7KOPO
## zip_type is STD
## uid is CSM7KOPO
## file is SM-74Y8S.tar
## file prefix is SM-74Y8S
## zip_type is VIR
## uid is SM-74Y8S
## file is SM-76CAQ.tar
## file prefix is SM-76CAQ
## zip_type is VIR
## uid is SM-76CAQ
## file is SM-7EDB2.tar
## file prefix is SM-7EDB2
## zip_type is VIR
## uid is SM-7EDB2
## file is SM-7M8SX.tar
## file prefix is SM-7M8SX
## zip_type is VIR
## uid is SM-7M8SX
## file is SM-9WOBU.tar
## file prefix is SM-9WOBU
## zip_type is VIR
## uid is SM-9WOBU
## file is CSM67U9V_P.fastq.gz
## file prefix is CSM67U9V_P
## zip_type is P
## uid is CSM67U9V_P
## file is CSM67UAU.tar
## file prefix is CSM67UAU
## zip_type is STD
## uid is CSM67UAU
## file is CSM67UAW.tar
## file prefix is CSM67UAW
## zip_type is STD
## uid is CSM67UAW
## file is CSM67UAY.tar
## file prefix is CSM67UAY
## zip_type is STD
## uid is CSM67UAY
## file is CSM67UB1.tar
## file prefix is CSM67UB1
## zip_type is STD
## uid is CSM67UB1
## file is CSM67UB3.tar
## file prefix is CSM67UB3
## zip_type is STD
## uid is CSM67UB3
## file is CSM79HIR.tar
## file prefix is CSM79HIR
## zip_type is STD
## uid is CSM79HIR
## file is CSM79HIT.tar
## file prefix is CSM79HIT
## zip_type is STD
## uid is CSM79HIT
## file is CSM79HIV.tar
## file prefix is CSM79HIV
## zip_type is STD
## uid is CSM79HIV
## file is CSM79HIX.tar
## file prefix is CSM79HIX
## zip_type is STD
## uid is CSM79HIX
## file is CSM79HIZ.tar
## file prefix is CSM79HIZ
## zip_type is STD
## uid is CSM79HIZ
## file is CSM7KOK3.tar
## file prefix is CSM7KOK3
## zip_type is STD
## uid is CSM7KOK3
## file is CSM7KOK5.tar
## file prefix is CSM7KOK5
## zip_type is STD
## uid is CSM7KOK5
## file is CSM7KOK7.tar
## file prefix is CSM7KOK7
## zip_type is STD
## uid is CSM7KOK7
## file is CSM7KOKB.tar
## file prefix is CSM7KOKB
## zip_type is STD
## uid is CSM7KOKB
## file is CSM7KOKD.tar
## file prefix is CSM7KOKD
## zip_type is STD
## uid is CSM7KOKD
## file is CSM7KOPS.tar
## file prefix is CSM7KOPS
## zip_type is STD
## uid is CSM7KOPS
## file is CSM7KOPU.tar
## file prefix is CSM7KOPU
## zip_type is STD
## uid is CSM7KOPU
## file is CSM7KOPW.tar
## file prefix is CSM7KOPW
## zip_type is STD
## uid is CSM7KOPW
## file is CSM7KOQ1.tar
## file prefix is CSM7KOQ1
## zip_type is STD
## uid is CSM7KOQ1
## file is SM-73JXY.tar
## file prefix is SM-73JXY
## zip_type is VIR
## uid is SM-73JXY
## file is SM-76EON.tar
## file prefix is SM-76EON
## zip_type is VIR
## uid is SM-76EON
## file is SM-791BO.tar
## file prefix is SM-791BO
## zip_type is VIR
## uid is SM-791BO
## file is SM-7CRX7.tar
## file prefix is SM-7CRX7
## zip_type is VIR
## uid is SM-7CRX7
## file is SM-7F43V.tar
## file prefix is SM-7F43V
## zip_type is VIR
## uid is SM-7F43V
## file is SM-7HAI3.tar
## file prefix is SM-7HAI3
## zip_type is VIR
## uid is SM-7HAI3
## file is SM-7IWFI.tar
## file prefix is SM-7IWFI
## zip_type is VIR
## uid is SM-7IWFI
## file is SM-7L29M.tar
## file prefix is SM-7L29M
## zip_type is VIR
## uid is SM-7L29M
## file is SM-7MCUT.tar
## file prefix is SM-7MCUT
## zip_type is VIR
## uid is SM-7MCUT
## file is SM-9PJ1H.tar
## file prefix is SM-9PJ1H
## zip_type is VIR
## uid is SM-9PJ1H
## file is SM-9RYIR.tar
## file prefix is SM-9RYIR
## zip_type is VIR
## uid is SM-9RYIR
## file is SM-9UW5R.tar
## file prefix is SM-9UW5R
## zip_type is VIR
## uid is SM-9UW5R
## file is SM-9WO25.tar
## file prefix is SM-9WO25
## zip_type is VIR
## uid is SM-9WO25
## file is SM-9Y7E4.tar
## file prefix is SM-9Y7E4
## zip_type is VIR
## uid is SM-9Y7E4
## file is SM-A18Z6.tar
## file prefix is SM-A18Z6
## zip_type is VIR
## uid is SM-A18Z6
## file is SM-A61QH.tar
## file prefix is SM-A61QH
## zip_type is VIR
## uid is SM-A61QH
## file is CSM67U9X_P.fastq.gz
## file prefix is CSM67U9X_P
## zip_type is P
## uid is CSM67U9X_P
## file is CSM67UB5_P.fastq.gz
## file prefix is CSM67UB5_P
## zip_type is P
## uid is CSM67UB5_P
## file is CSM67UB7_P.fastq.gz
## file prefix is CSM67UB7_P
## zip_type is P
## uid is CSM67UB7_P
## file is CSM67UB9_P.fastq.gz
## file prefix is CSM67UB9_P
## zip_type is P
## uid is CSM67UB9_P
## file is CSM67UB9.tar
## file prefix is CSM67UB9
## zip_type is STD
## uid is CSM67UB9
## file is CSM67UBB.tar
## file prefix is CSM67UBB
## zip_type is STD
## uid is CSM67UBB
## file is CSM79HJ2_P.fastq.gz
## file prefix is CSM79HJ2_P
## zip_type is P
## uid is CSM79HJ2_P
## file is CSM79HJ4_P.fastq.gz
## file prefix is CSM79HJ4_P
## zip_type is P
## uid is CSM79HJ4_P
## file is CSM79HJ6_P.fastq.gz
## file prefix is CSM79HJ6_P
## zip_type is P
## uid is CSM79HJ6_P
## file is CSM79HJ8_P.fastq.gz
## file prefix is CSM79HJ8_P
## zip_type is P
## uid is CSM79HJ8_P
## file is CSM79HJA.tar
## file prefix is CSM79HJA
## zip_type is STD
## uid is CSM79HJA
## file is CSM79HJC_P.fastq.gz
## file prefix is CSM79HJC_P
## zip_type is P
## uid is CSM79HJC_P
## file is CSM7KOKF.tar
## file prefix is CSM7KOKF
## zip_type is STD
## uid is CSM7KOKF
## file is CSM7KOKH_P.fastq.gz
## file prefix is CSM7KOKH_P
## zip_type is P
## uid is CSM7KOKH_P
## file is CSM7KOKJ.tar
## file prefix is CSM7KOKJ
## zip_type is STD
## uid is CSM7KOKJ
## file is CSM7KOKL_P.fastq.gz
## file prefix is CSM7KOKL_P
## zip_type is P
## uid is CSM7KOKL_P
## file is CSM7KOKN.tar
## file prefix is CSM7KOKN
## zip_type is STD
## uid is CSM7KOKN
## file is CSM7KOKP_P.fastq.gz
## file prefix is CSM7KOKP_P
## zip_type is P
## uid is CSM7KOKP_P
## file is CSM7KOQX.tar
## file prefix is CSM7KOQX
## zip_type is STD
## uid is CSM7KOQX
## file is CSM7KOQZ_P.fastq.gz
## file prefix is CSM7KOQZ_P
## zip_type is P
## uid is CSM7KOQZ_P
## file is CSM7KOR2.tar
## file prefix is CSM7KOR2
## zip_type is STD
## uid is CSM7KOR2
## file is CSM7KOR4_P.fastq.gz
## file prefix is CSM7KOR4_P
## zip_type is P
## uid is CSM7KOR4_P
## file is CSM7KOR8_P.fastq.gz
## file prefix is CSM7KOR8_P
## zip_type is P
## uid is CSM7KOR8_P
## file is SM-791C2.tar
## file prefix is SM-791C2
## zip_type is VIR
## uid is SM-791C2
## file is SM-7CRJC.tar
## file prefix is SM-7CRJC
## zip_type is VIR
## uid is SM-7CRJC
## file is SM-7KPVK.tar
## file prefix is SM-7KPVK
## zip_type is VIR
## uid is SM-7KPVK
## file is SM-9JUDG.tar
## file prefix is SM-9JUDG
## zip_type is VIR
## uid is SM-9JUDG
## file is SM-9OS8F.tar
## file prefix is SM-9OS8F
## zip_type is VIR
## uid is SM-9OS8F
## file is SM-9Y7D8.tar
## file prefix is SM-9Y7D8
## zip_type is VIR
## uid is SM-9Y7D8
## file is CSM67UCU_P.fastq.gz
## file prefix is CSM67UCU_P
## zip_type is P
## uid is CSM67UCU_P
## file is CSM67UAI_P.fastq.gz
## file prefix is CSM67UAI_P
## zip_type is P
## uid is CSM67UAI_P
## file is CSM79HHO.tar
## file prefix is CSM79HHO
## zip_type is STD
## uid is CSM79HHO
## file is CSM79HHM.tar
## file prefix is CSM79HHM
## zip_type is STD
## uid is CSM79HHM
## file is CSM79HHU.tar
## file prefix is CSM79HHU
## zip_type is STD
## uid is CSM79HHU
## file is CSM79HN2.tar
## file prefix is CSM79HN2
## zip_type is STD
## uid is CSM79HN2
## file is CSM79HN6.tar
## file prefix is CSM79HN6
## zip_type is STD
## uid is CSM79HN6
## file is CSM7KOLK.tar
## file prefix is CSM7KOLK
## zip_type is STD
## uid is CSM7KOLK
## file is CSM7KOLM.tar
## file prefix is CSM7KOLM
## zip_type is STD
## uid is CSM7KOLM
## file is CSM7KOSV.tar
## file prefix is CSM7KOSV
## zip_type is STD
## uid is CSM7KOSV
## file is CSM7KOSX.tar
## file prefix is CSM7KOSX
## zip_type is STD
## uid is CSM7KOSX
## file is SM-7D34P.tar
## file prefix is SM-7D34P
## zip_type is VIR
## uid is SM-7D34P
## file is SM-7F444.tar
## file prefix is SM-7F444
## zip_type is VIR
## uid is SM-7F444
## file is SM-7HXKA.tar
## file prefix is SM-7HXKA
## zip_type is VIR
## uid is SM-7HXKA
## file is SM-7NSOY.tar
## file prefix is SM-7NSOY
## zip_type is VIR
## uid is SM-7NSOY
## file is SM-9RD3E.tar
## file prefix is SM-9RD3E
## zip_type is VIR
## uid is SM-9RD3E
## file is SM-9ZE8U.tar
## file prefix is SM-9ZE8U
## zip_type is VIR
## uid is SM-9ZE8U
## file is SM-A1HUO.tar
## file prefix is SM-A1HUO
## zip_type is VIR
## uid is SM-A1HUO
## file is CSM67UH7.tar
## file prefix is CSM67UH7
## zip_type is STD
## uid is CSM67UH7
## file is CSM79HGR_P.fastq.gz
## file prefix is CSM79HGR_P
## zip_type is P
## uid is CSM79HGR_P
## file is CSM79HGV_P.fastq.gz
## file prefix is CSM79HGV_P
## zip_type is P
## uid is CSM79HGV_P
## file is CSM79HGX.tar
## file prefix is CSM79HGX
## zip_type is STD
## uid is CSM79HGX
## file is CSM79HGZ.tar
## file prefix is CSM79HGZ
## zip_type is STD
## uid is CSM79HGZ
## file is CSM79HOV.tar
## file prefix is CSM79HOV
## zip_type is STD
## uid is CSM79HOV
## file is CSM79HOX.tar
## file prefix is CSM79HOX
## zip_type is STD
## uid is CSM79HOX
## file is CSM79HOZ.tar
## file prefix is CSM79HOZ
## zip_type is STD
## uid is CSM79HOZ
## file is CSM79HP2.tar
## file prefix is CSM79HP2
## zip_type is STD
## uid is CSM79HP2
## file is CSM79HP4.tar
## file prefix is CSM79HP4
## zip_type is STD
## uid is CSM79HP4
## file is CSM79HP6.tar
## file prefix is CSM79HP6
## zip_type is STD
## uid is CSM79HP6
## file is CSM7KOOH.tar
## file prefix is CSM7KOOH
## zip_type is STD
## uid is CSM7KOOH
## file is CSM7KOOJ.tar
## file prefix is CSM7KOOJ
## zip_type is STD
## uid is CSM7KOOJ
## file is CSM7KOOL.tar
## file prefix is CSM7KOOL
## zip_type is STD
## uid is CSM7KOOL
## file is CSM7KOON.tar
## file prefix is CSM7KOON
## zip_type is STD
## uid is CSM7KOON
## file is CSM7KOOP.tar
## file prefix is CSM7KOOP
## zip_type is STD
## uid is CSM7KOOP
## file is CSM7KOOR.tar
## file prefix is CSM7KOOR
## zip_type is STD
## uid is CSM7KOOR
## file is CSMA8M9R.tar
## file prefix is CSMA8M9R
## zip_type is STD
## uid is CSMA8M9R
## file is CSMAAEUA.tar
## file prefix is CSMAAEUA
## zip_type is STD
## uid is CSMAAEUA
## file is CSMAE44D.tar
## file prefix is CSMAE44D
## zip_type is STD
## uid is CSMAE44D
## file is CSMAG78W.tar
## file prefix is CSMAG78W
## zip_type is STD
## uid is CSMAG78W
## file is CSMAHYLR.tar
## file prefix is CSMAHYLR
## zip_type is STD
## uid is CSMAHYLR
## file is SM-7H4HE.tar
## file prefix is SM-7H4HE
## zip_type is VIR
## uid is SM-7H4HE
## file is SM-7IWFB.tar
## file prefix is SM-7IWFB
## zip_type is VIR
## uid is SM-7IWFB
## file is SM-7ORRS.tar
## file prefix is SM-7ORRS
## zip_type is VIR
## uid is SM-7ORRS
## file is SM-9KOPT.tar
## file prefix is SM-9KOPT
## zip_type is VIR
## uid is SM-9KOPT
## file is SM-9SD8W.tar
## file prefix is SM-9SD8W
## zip_type is VIR
## uid is SM-9SD8W
## file is SM-A3J7G.tar
## file prefix is SM-A3J7G
## zip_type is VIR
## uid is SM-A3J7G
## file is SM-AG791.tar
## file prefix is SM-AG791
## zip_type is VIR
## uid is SM-AG791
## file is CSM79HH2_P.fastq.gz
## file prefix is CSM79HH2_P
## zip_type is P
## uid is CSM79HH2_P
## file is CSM79HH4.tar
## file prefix is CSM79HH4
## zip_type is STD
## uid is CSM79HH4
## file is CSM79HH8.tar
## file prefix is CSM79HH8
## zip_type is STD
## uid is CSM79HH8
## file is CSM79HHA.tar
## file prefix is CSM79HHA
## zip_type is STD
## uid is CSM79HHA
## file is CSM79HPA_TR.tar
## file prefix is CSM79HPA_TR
## zip_type is TR
## uid is CSM79HPA_TR
## file is CSM79HPA.tar
## file prefix is CSM79HPA
## zip_type is STD
## uid is CSM79HPA
## file is CSM79HPC.tar
## file prefix is CSM79HPC
## zip_type is STD
## uid is CSM79HPC
## file is CSM7KONA.tar
## file prefix is CSM7KONA
## zip_type is STD
## uid is CSM7KONA
## file is CSM7KONK.tar
## file prefix is CSM7KONK
## zip_type is STD
## uid is CSM7KONK
## file is CSM7KOTA.tar
## file prefix is CSM7KOTA
## zip_type is STD
## uid is CSM7KOTA
## file is CSM7KOTC.tar
## file prefix is CSM7KOTC
## zip_type is STD
## uid is CSM7KOTC
## file is CSM7KOTK.tar
## file prefix is CSM7KOTK
## zip_type is STD
## uid is CSM7KOTK
## file is SM-7EDAY.tar
## file prefix is SM-7EDAY
## zip_type is VIR
## uid is SM-7EDAY
## file is SM-7GYK1.tar
## file prefix is SM-7GYK1
## zip_type is VIR
## uid is SM-7GYK1
## file is SM-7M8RY.tar
## file prefix is SM-7M8RY
## zip_type is VIR
## uid is SM-7M8RY
## file is SM-7ORRI.tar
## file prefix is SM-7ORRI
## zip_type is VIR
## uid is SM-7ORRI
## file is SM-9RRZK.tar
## file prefix is SM-9RRZK
## zip_type is VIR
## uid is SM-9RRZK
## file is SM-A12L1.tar
## file prefix is SM-A12L1
## zip_type is VIR
## uid is SM-A12L1
## file is SM-A5QM2.tar
## file prefix is SM-A5QM2
## zip_type is VIR
## uid is SM-A5QM2
## file is SM-AFSIL.tar
## file prefix is SM-AFSIL
## zip_type is VIR
## uid is SM-AFSIL
## file is CSM79HG7_P.fastq.gz
## file prefix is CSM79HG7_P
## zip_type is P
## uid is CSM79HG7_P
## file is CSM79HHW_P.fastq.gz
## file prefix is CSM79HHW_P
## zip_type is P
## uid is CSM79HHW_P
## file is CSM79HHW.tar
## file prefix is CSM79HHW
## zip_type is STD
## uid is CSM79HHW
## file is CSM79HJM.tar
## file prefix is CSM79HJM
## zip_type is STD
## uid is CSM79HJM
## file is CSM79HJO.tar
## file prefix is CSM79HJO
## zip_type is STD
## uid is CSM79HJO
## file is CSM79HJQ.tar
## file prefix is CSM79HJQ
## zip_type is STD
## uid is CSM79HJQ
## file is CSM79HJS.tar
## file prefix is CSM79HJS
## zip_type is STD
## uid is CSM79HJS
## file is CSM79HJU.tar
## file prefix is CSM79HJU
## zip_type is STD
## uid is CSM79HJU
## file is CSM79HQV.tar
## file prefix is CSM79HQV
## zip_type is STD
## uid is CSM79HQV
## file is CSM79HQX.tar
## file prefix is CSM79HQX
## zip_type is STD
## uid is CSM79HQX
## file is CSM79HQZ.tar
## file prefix is CSM79HQZ
## zip_type is STD
## uid is CSM79HQZ
## file is CSM79HR2.tar
## file prefix is CSM79HR2
## zip_type is STD
## uid is CSM79HR2
## file is CSM79HR4.tar
## file prefix is CSM79HR4
## zip_type is STD
## uid is CSM79HR4
## file is CSM79HR6.tar
## file prefix is CSM79HR6
## zip_type is STD
## uid is CSM79HR6
## file is CSM7KOOT.tar
## file prefix is CSM7KOOT
## zip_type is STD
## uid is CSM7KOOT
## file is CSM7KOOV.tar
## file prefix is CSM7KOOV
## zip_type is STD
## uid is CSM7KOOV
## file is CSM7KOOX.tar
## file prefix is CSM7KOOX
## zip_type is STD
## uid is CSM7KOOX
## file is CSM7KOOZ.tar
## file prefix is CSM7KOOZ
## zip_type is STD
## uid is CSM7KOOZ
## file is CSM7KOP2.tar
## file prefix is CSM7KOP2
## zip_type is STD
## uid is CSM7KOP2
## file is CSMA88CB.tar
## file prefix is CSMA88CB
## zip_type is STD
## uid is CSMA88CB
## file is CSMA9J65.tar
## file prefix is CSMA9J65
## zip_type is STD
## uid is CSMA9J65
## file is CSMACTZP.tar
## file prefix is CSMACTZP
## zip_type is STD
## uid is CSMACTZP
## file is CSMAF72L.tar
## file prefix is CSMAF72L
## zip_type is STD
## uid is CSMAF72L
## file is CSMAH393.tar
## file prefix is CSMAH393
## zip_type is STD
## uid is CSMAH393
## file is CSMAIG7X.tar
## file prefix is CSMAIG7X
## zip_type is STD
## uid is CSMAIG7X
## file is SM-7EWUT.tar
## file prefix is SM-7EWUT
## zip_type is VIR
## uid is SM-7EWUT
## file is SM-7FK4R.tar
## file prefix is SM-7FK4R
## zip_type is VIR
## uid is SM-7FK4R
## file is SM-7K1V8.tar
## file prefix is SM-7K1V8
## zip_type is VIR
## uid is SM-7K1V8
## file is SM-9KONQ.tar
## file prefix is SM-9KONQ
## zip_type is VIR
## uid is SM-9KONQ
## file is SM-9QMNN.tar
## file prefix is SM-9QMNN
## zip_type is VIR
## uid is SM-9QMNN
## file is SM-9W3B5.tar
## file prefix is SM-9W3B5
## zip_type is VIR
## uid is SM-9W3B5
## file is SM-A9J9P.tar
## file prefix is SM-A9J9P
## zip_type is VIR
## uid is SM-A9J9P
## file is SM-AIG8N.tar
## file prefix is SM-AIG8N
## zip_type is VIR
## uid is SM-AIG8N
## file is CSM79HKX.tar
## file prefix is CSM79HKX
## zip_type is STD
## uid is CSM79HKX
## file is CSM79HKZ.tar
## file prefix is CSM79HKZ
## zip_type is STD
## uid is CSM79HKZ
## file is CSM79HL4.tar
## file prefix is CSM79HL4
## zip_type is STD
## uid is CSM79HL4
## file is CSM79HL6.tar
## file prefix is CSM79HL6
## zip_type is STD
## uid is CSM79HL6
## file is CSM7KOKR.tar
## file prefix is CSM7KOKR
## zip_type is STD
## uid is CSM7KOKR
## file is CSM7KOKT.tar
## file prefix is CSM7KOKT
## zip_type is STD
## uid is CSM7KOKT
## file is CSM7KOKZ.tar
## file prefix is CSM7KOKZ
## zip_type is STD
## uid is CSM7KOKZ
## file is CSM7KOL2.tar
## file prefix is CSM7KOL2
## zip_type is STD
## uid is CSM7KOL2
## file is CSM7KORC.tar
## file prefix is CSM7KORC
## zip_type is STD
## uid is CSM7KORC
## file is CSM7KORK.tar
## file prefix is CSM7KORK
## zip_type is STD
## uid is CSM7KORK
## file is CSM7KORI.tar
## file prefix is CSM7KORI
## zip_type is STD
## uid is CSM7KORI
## file is CSM7KORG.tar
## file prefix is CSM7KORG
## zip_type is STD
## uid is CSM7KORG
## file is SM-7LDIY.tar
## file prefix is SM-7LDIY
## zip_type is VIR
## uid is SM-7LDIY
## file is SM-9J5JL.tar
## file prefix is SM-9J5JL
## zip_type is VIR
## uid is SM-9J5JL
## file is SM-9T55N.tar
## file prefix is SM-9T55N
## zip_type is VIR
## uid is SM-9T55N
## file is SM-9Y4A9.tar
## file prefix is SM-9Y4A9
## zip_type is VIR
## uid is SM-9Y4A9
## file is SM-A8XLN.tar
## file prefix is SM-A8XLN
## zip_type is VIR
## uid is SM-A8XLN
## file is SM-AEYQ9.tar
## file prefix is SM-AEYQ9
## zip_type is VIR
## uid is SM-AEYQ9
## file is CSM79HNE.tar
## file prefix is CSM79HNE
## zip_type is STD
## uid is CSM79HNE
## file is CSM79HNG_P.fastq.gz
## file prefix is CSM79HNG_P
## zip_type is P
## uid is CSM79HNG_P
## file is CSM79HNI.tar
## file prefix is CSM79HNI
## zip_type is STD
## uid is CSM79HNI
## file is CSM79HNK.tar
## file prefix is CSM79HNK
## zip_type is STD
## uid is CSM79HNK
## file is CSM79HNM.tar
## file prefix is CSM79HNM
## zip_type is STD
## uid is CSM79HNM
## file is CSM7KOLY.tar
## file prefix is CSM7KOLY
## zip_type is STD
## uid is CSM7KOLY
## file is CSM7KOS7.tar
## file prefix is CSM7KOS7
## zip_type is STD
## uid is CSM7KOS7
## file is CSM7KOSH.tar
## file prefix is CSM7KOSH
## zip_type is STD
## uid is CSM7KOSH
## file is CSM9X1Y5.tar
## file prefix is CSM9X1Y5
## zip_type is STD
## uid is CSM9X1Y5
## file is SM-9SIJ7.tar
## file prefix is SM-9SIJ7
## zip_type is VIR
## uid is SM-9SIJ7
## file is SM-A1ZVM.tar
## file prefix is SM-A1ZVM
## zip_type is VIR
## uid is SM-A1ZVM
## file is SM-AGUNO.tar
## file prefix is SM-AGUNO
## zip_type is VIR
## uid is SM-AGUNO
## file is CSM79HJI_P.fastq.gz
## file prefix is CSM79HJI_P
## zip_type is P
## uid is CSM79HJI_P
## file is CSM79HNO.tar
## file prefix is CSM79HNO
## zip_type is STD
## uid is CSM79HNO
## file is CSM79HNU.tar
## file prefix is CSM79HNU
## zip_type is STD
## uid is CSM79HNU
## file is CSM79HNW.tar
## file prefix is CSM79HNW
## zip_type is STD
## uid is CSM79HNW
## file is CSM7KOMB.tar
## file prefix is CSM7KOMB
## zip_type is STD
## uid is CSM7KOMB
## file is CSM7KOMH.tar
## file prefix is CSM7KOMH
## zip_type is STD
## uid is CSM7KOMH
## file is CSM7KOSL.tar
## file prefix is CSM7KOSL
## zip_type is STD
## uid is CSM7KOSL
## file is CSM7KOSP.tar
## file prefix is CSM7KOSP
## zip_type is STD
## uid is CSM7KOSP
## file is CSM7KOST.tar
## file prefix is CSM7KOST
## zip_type is STD
## uid is CSM7KOST
## file is CSM7KOSJ.tar
## file prefix is CSM7KOSJ
## zip_type is STD
## uid is CSM7KOSJ
## file is SM-7K1WM.tar
## file prefix is SM-7K1WM
## zip_type is VIR
## uid is SM-7K1WM
## file is SM-7LDJA.tar
## file prefix is SM-7LDJA
## zip_type is VIR
## uid is SM-7LDJA
## file is SM-7R3AL.tar
## file prefix is SM-7R3AL
## zip_type is VIR
## uid is SM-7R3AL
## file is SM-9RD3I.tar
## file prefix is SM-9RD3I
## zip_type is VIR
## uid is SM-9RD3I
## file is SM-A1ZX8.tar
## file prefix is SM-A1ZX8
## zip_type is VIR
## uid is SM-A1ZX8
## file is SM-A77UX.tar
## file prefix is SM-A77UX
## zip_type is VIR
## uid is SM-A77UX
## file is CSM79HNY_P.fastq.gz
## file prefix is CSM79HNY_P
## zip_type is P
## uid is CSM79HNY_P
## file is CSM79HR8.tar
## file prefix is CSM79HR8
## zip_type is STD
## uid is CSM79HR8
## file is CSM79HRC.tar
## file prefix is CSM79HRC
## zip_type is STD
## uid is CSM79HRC
## file is CSM79HRE.tar
## file prefix is CSM79HRE
## zip_type is STD
## uid is CSM79HRE
## file is CSM79HRG.tar
## file prefix is CSM79HRG
## zip_type is STD
## uid is CSM79HRG
## file is CSM7KOO5.tar
## file prefix is CSM7KOO5
## zip_type is STD
## uid is CSM7KOO5
## file is CSM7KOO9.tar
## file prefix is CSM7KOO9
## zip_type is STD
## uid is CSM7KOO9
## file is CSM7KOOD.tar
## file prefix is CSM7KOOD
## zip_type is STD
## uid is CSM7KOOD
## file is CSM7KOOF.tar
## file prefix is CSM7KOOF
## zip_type is STD
## uid is CSM7KOOF
## file is CSM9X1Z4.tar
## file prefix is CSM9X1Z4
## zip_type is STD
## uid is CSM9X1Z4
## file is CSM9X1ZC.tar
## file prefix is CSM9X1ZC
## zip_type is STD
## uid is CSM9X1ZC
## file is SM-9HV41.tar
## file prefix is SM-9HV41
## zip_type is VIR
## uid is SM-9HV41
## file is SM-9JYAT.tar
## file prefix is SM-9JYAT
## zip_type is VIR
## uid is SM-9JYAT
## file is SM-9O9QZ.tar
## file prefix is SM-9O9QZ
## zip_type is VIR
## uid is SM-9O9QZ
## file is SM-9SVR3.tar
## file prefix is SM-9SVR3
## zip_type is VIR
## uid is SM-9SVR3
## file is SM-9YTMB.tar
## file prefix is SM-9YTMB
## zip_type is VIR
## uid is SM-9YTMB
## file is SM-A55Y2.tar
## file prefix is SM-A55Y2
## zip_type is VIR
## uid is SM-A55Y2
## file is SM-AFSIQ.tar
## file prefix is SM-AFSIQ
## zip_type is VIR
## uid is SM-AFSIQ
## file is SM-ALGUT.tar
## file prefix is SM-ALGUT
## zip_type is VIR
## uid is SM-ALGUT
## file is SM-AVPGY.tar
## file prefix is SM-AVPGY
## zip_type is VIR
## uid is SM-AVPGY
## file is CSM7KOP6.tar
## file prefix is CSM7KOP6
## zip_type is STD
## uid is CSM7KOP6
## file is CSM7KOP8.tar
## file prefix is CSM7KOP8
## zip_type is STD
## uid is CSM7KOP8
## file is CSM7KOPE.tar
## file prefix is CSM7KOPE
## zip_type is STD
## uid is CSM7KOPE
## file is CSM7KOU7_P.fastq.gz
## file prefix is CSM7KOU7_P
## zip_type is P
## uid is CSM7KOU7_P
## file is CSM7KOU9.tar
## file prefix is CSM7KOU9
## zip_type is STD
## uid is CSM7KOU9
## file is CSM7KOUB.tar
## file prefix is CSM7KOUB
## zip_type is STD
## uid is CSM7KOUB
## file is CSM9X1ZO.tar
## file prefix is CSM9X1ZO
## zip_type is STD
## uid is CSM9X1ZO
## file is CSM9X1ZQ.tar
## file prefix is CSM9X1ZQ
## zip_type is STD
## uid is CSM9X1ZQ
## file is CSM9X1ZY.tar
## file prefix is CSM9X1ZY
## zip_type is STD
## uid is CSM9X1ZY
## file is CSM9X22G.tar
## file prefix is CSM9X22G
## zip_type is STD
## uid is CSM9X22G
## file is CSM9X22I.tar
## file prefix is CSM9X22I
## zip_type is STD
## uid is CSM9X22I
## file is CSM9X22K.tar
## file prefix is CSM9X22K
## zip_type is STD
## uid is CSM9X22K
## file is SM-9ZEPB.tar
## file prefix is SM-9ZEPB
## zip_type is VIR
## uid is SM-9ZEPB
## file is SM-A7J1W.tar
## file prefix is SM-A7J1W
## zip_type is VIR
## uid is SM-A7J1W
## file is SM-AY1XU.tar
## file prefix is SM-AY1XU
## zip_type is VIR
## uid is SM-AY1XU
## file is CSM7KOQP_P.fastq.gz
## file prefix is CSM7KOQP_P
## zip_type is P
## uid is CSM7KOQP_P
## file is CSM79HQR_P.fastq.gz
## file prefix is CSM79HQR_P
## zip_type is P
## uid is CSM79HQR_P
## file is CSM7KORO.tar
## file prefix is CSM7KORO
## zip_type is STD
## uid is CSM7KORO
## file is CSM7KORM.tar
## file prefix is CSM7KORM
## zip_type is STD
## uid is CSM7KORM
## file is CSM7KORU.tar
## file prefix is CSM7KORU
## zip_type is STD
## uid is CSM7KORU
## file is CSM7KORS.tar
## file prefix is CSM7KORS
## zip_type is STD
## uid is CSM7KORS
## file is CSM9X1XU.tar
## file prefix is CSM9X1XU
## zip_type is STD
## uid is CSM9X1XU
## file is CSM9X1Y3.tar
## file prefix is CSM9X1Y3
## zip_type is STD
## uid is CSM9X1Y3
## file is CSM9X21J.tar
## file prefix is CSM9X21J
## zip_type is STD
## uid is CSM9X21J
## file is CSM9X21L.tar
## file prefix is CSM9X21L
## zip_type is STD
## uid is CSM9X21L
## file is CSM9X21N.tar
## file prefix is CSM9X21N
## zip_type is STD
## uid is CSM9X21N
## file is SM-9WOBW.tar
## file prefix is SM-9WOBW
## zip_type is VIR
## uid is SM-9WOBW
## file is SM-9Y7E9.tar
## file prefix is SM-9Y7E9
## zip_type is VIR
## uid is SM-9Y7E9
## file is SM-A3743.tar
## file prefix is SM-A3743
## zip_type is VIR
## uid is SM-A3743
## file is SM-AEB6V.tar
## file prefix is SM-AEB6V
## zip_type is VIR
## uid is SM-AEB6V
## file is SM-ANU6U.tar
## file prefix is SM-ANU6U
## zip_type is VIR
## uid is SM-ANU6U
## file is SM-AVADM.tar
## file prefix is SM-AVADM
## zip_type is VIR
## uid is SM-AVADM
## file is CSM7KOTO.tar
## file prefix is CSM7KOTO
## zip_type is STD
## uid is CSM7KOTO
## file is CSM7KOTQ.tar
## file prefix is CSM7KOTQ
## zip_type is STD
## uid is CSM7KOTQ
## file is CSM7KOTS.tar
## file prefix is CSM7KOTS
## zip_type is STD
## uid is CSM7KOTS
## file is CSM7KOTU.tar
## file prefix is CSM7KOTU
## zip_type is STD
## uid is CSM7KOTU
## file is CSM9X1YV.tar
## file prefix is CSM9X1YV
## zip_type is STD
## uid is CSM9X1YV
## file is CSM9X21R.tar
## file prefix is CSM9X21R
## zip_type is STD
## uid is CSM9X21R
## file is CSM9X21T.tar
## file prefix is CSM9X21T
## zip_type is STD
## uid is CSM9X21T
## file is CSM9X222.tar
## file prefix is CSM9X222
## zip_type is STD
## uid is CSM9X222
## file is CSM9X233.tar
## file prefix is CSM9X233
## zip_type is STD
## uid is CSM9X233
## file is CSM9X235.tar
## file prefix is CSM9X235
## zip_type is STD
## uid is CSM9X235
## file is CSM9X237.tar
## file prefix is CSM9X237
## zip_type is STD
## uid is CSM9X237
## file is CSM9X23B.tar
## file prefix is CSM9X23B
## zip_type is STD
## uid is CSM9X23B
## file is SM-A8MIU.tar
## file prefix is SM-A8MIU
## zip_type is VIR
## uid is SM-A8MIU
## file is SM-AC5QV.tar
## file prefix is SM-AC5QV
## zip_type is VIR
## uid is SM-AC5QV
## file is SM-AFSIF.tar
## file prefix is SM-AFSIF
## zip_type is VIR
## uid is SM-AFSIF
## file is SM-AJDM7.tar
## file prefix is SM-AJDM7
## zip_type is VIR
## uid is SM-AJDM7
## file is SM-BXXS9.tar
## file prefix is SM-BXXS9
## zip_type is VIR
## uid is SM-BXXS9
## file is SM-C1WB3.tar
## file prefix is SM-C1WB3
## zip_type is VIR
## uid is SM-C1WB3
## file is SM-CH32O.tar
## file prefix is SM-CH32O
## zip_type is VIR
## uid is SM-CH32O
## file is CSM7KOQ5_P.fastq.gz
## file prefix is CSM7KOQ5_P
## zip_type is P
## uid is CSM7KOQ5_P
## file is CSM7KOUJ_P.fastq.gz
## file prefix is CSM7KOUJ_P
## zip_type is P
## uid is CSM7KOUJ_P
## file is CSM7KOUL.tar
## file prefix is CSM7KOUL
## zip_type is STD
## uid is CSM7KOUL
## file is CSM7KOUN.tar
## file prefix is CSM7KOUN
## zip_type is STD
## uid is CSM7KOUN
## file is CSM9X219.tar
## file prefix is CSM9X219
## zip_type is STD
## uid is CSM9X219
## file is CSM9X213.tar
## file prefix is CSM9X213
## zip_type is STD
## uid is CSM9X213
## file is CSM9X215.tar
## file prefix is CSM9X215
## zip_type is STD
## uid is CSM9X215
## file is CSM9X211.tar
## file prefix is CSM9X211
## zip_type is STD
## uid is CSM9X211
## file is CSM9X22S.tar
## file prefix is CSM9X22S
## zip_type is STD
## uid is CSM9X22S
## file is CSM9X22U.tar
## file prefix is CSM9X22U
## zip_type is STD
## uid is CSM9X22U
## file is CSM9X23H.tar
## file prefix is CSM9X23H
## zip_type is STD
## uid is CSM9X23H
## file is CSM9X23N.tar
## file prefix is CSM9X23N
## zip_type is STD
## uid is CSM9X23N
## file is SM-A47S2.tar
## file prefix is SM-A47S2
## zip_type is VIR
## uid is SM-A47S2
## file is SM-A9F74.tar
## file prefix is SM-A9F74
## zip_type is VIR
## uid is SM-A9F74
## file is SM-AF6NS.tar
## file prefix is SM-AF6NS
## zip_type is VIR
## uid is SM-AF6NS
## file is SM-AGUNG.tar
## file prefix is SM-AGUNG
## zip_type is VIR
## uid is SM-AGUNG
## file is SM-B2OQT.tar
## file prefix is SM-B2OQT
## zip_type is VIR
## uid is SM-B2OQT
## file is SM-B5F9D.tar
## file prefix is SM-B5F9D
## zip_type is VIR
## uid is SM-B5F9D
## file is SM-CF3XF.tar
## file prefix is SM-CF3XF
## zip_type is VIR
## uid is SM-CF3XF
## file is SM-CL4HS.tar
## file prefix is SM-CL4HS
## zip_type is VIR
## uid is SM-CL4HS
## file is ESM5MEDZ_P.fastq.gz
## file prefix is ESM5MEDZ_P
## zip_type is P
## uid is ESM5MEDZ_P
## file is ESM5MEE2_P.fastq.gz
## file prefix is ESM5MEE2_P
## zip_type is P
## uid is ESM5MEE2_P
## file is HSM5MEE5_P.fastq.gz
## file prefix is HSM5MEE5_P
## zip_type is P
## uid is HSM5MEE5_P
## file is ESM5MEE6_P.fastq.gz
## file prefix is ESM5MEE6_P
## zip_type is P
## uid is ESM5MEE6_P
## file is ESM5ME9D_P.fastq.gz
## file prefix is ESM5ME9D_P
## zip_type is P
## uid is ESM5ME9D_P
## file is ESM5ME9G_P.fastq.gz
## file prefix is ESM5ME9G_P
## zip_type is P
## uid is ESM5ME9G_P
## file is ESM5MEBA_P.fastq.gz
## file prefix is ESM5MEBA_P
## zip_type is P
## uid is ESM5MEBA_P
## file is ESM5MEBE.tar
## file prefix is ESM5MEBE
## zip_type is STD
## uid is ESM5MEBE
## file is ESM5GEXY.tar
## file prefix is ESM5GEXY
## zip_type is STD
## uid is ESM5GEXY
## file is ESM5MEBG.tar
## file prefix is ESM5MEBG
## zip_type is STD
## uid is ESM5MEBG
## file is ESM5MECQ.tar
## file prefix is ESM5MECQ
## zip_type is STD
## uid is ESM5MECQ
## file is ESM5MEBI.tar
## file prefix is ESM5MEBI
## zip_type is STD
## uid is ESM5MEBI
## file is ESM5MECL.tar
## file prefix is ESM5MECL
## zip_type is STD
## uid is ESM5MECL
## file is SM-6ZKS7.tar
## file prefix is SM-6ZKS7
## zip_type is VIR
## uid is SM-6ZKS7
## file is SM-76ENX.tar
## file prefix is SM-76ENX
## zip_type is VIR
## uid is SM-76ENX
## file is SM-7DGV6.tar
## file prefix is SM-7DGV6
## zip_type is VIR
## uid is SM-7DGV6
## file is ESM5GEYU_P.fastq.gz
## file prefix is ESM5GEYU_P
## zip_type is P
## uid is ESM5GEYU_P
## file is ESM5GEZ1_P.fastq.gz
## file prefix is ESM5GEZ1_P
## zip_type is P
## uid is ESM5GEZ1_P
## file is ESM5GEZ3_P.fastq.gz
## file prefix is ESM5GEZ3_P
## zip_type is P
## uid is ESM5GEZ3_P
## file is ESM5MEBP_P.fastq.gz
## file prefix is ESM5MEBP_P
## zip_type is P
## uid is ESM5MEBP_P
## file is ESM5ME9H_P.fastq.gz
## file prefix is ESM5ME9H_P
## zip_type is P
## uid is ESM5ME9H_P
## file is ESM5GEYX_P.fastq.gz
## file prefix is ESM5GEYX_P
## zip_type is P
## uid is ESM5GEYX_P
## file is ESM5MEEJ_P.fastq.gz
## file prefix is ESM5MEEJ_P
## zip_type is P
## uid is ESM5MEEJ_P
## file is ESM5GEZ4_P.fastq.gz
## file prefix is ESM5GEZ4_P
## zip_type is P
## uid is ESM5GEZ4_P
## file is ESM5GEZ6_P.fastq.gz
## file prefix is ESM5GEZ6_P
## zip_type is P
## uid is ESM5GEZ6_P
## file is ESM5GEZA_P.fastq.gz
## file prefix is ESM5GEZA_P
## zip_type is P
## uid is ESM5GEZA_P
## file is ESM5MEBS.tar
## file prefix is ESM5MEBS
## zip_type is STD
## uid is ESM5MEBS
## file is ESM5MEBU.tar
## file prefix is ESM5MEBU
## zip_type is STD
## uid is ESM5MEBU
## file is ESM5MEC3.tar
## file prefix is ESM5MEC3
## zip_type is STD
## uid is ESM5MEC3
## file is ESM5MEC5.tar
## file prefix is ESM5MEC5
## zip_type is STD
## uid is ESM5MEC5
## file is ESM5MEDU.tar
## file prefix is ESM5MEDU
## zip_type is STD
## uid is ESM5MEDU
## file is ESM5MEC9.tar
## file prefix is ESM5MEC9
## zip_type is STD
## uid is ESM5MEC9
## file is ESM718SY.tar
## file prefix is ESM718SY
## zip_type is STD
## uid is ESM718SY
## file is ESM5ME9U.tar
## file prefix is ESM5ME9U
## zip_type is STD
## uid is ESM5ME9U
## file is SM-6V2DC.tar
## file prefix is SM-6V2DC
## zip_type is VIR
## uid is SM-6V2DC
## file is SM-6XJTO.tar
## file prefix is SM-6XJTO
## zip_type is VIR
## uid is SM-6XJTO
## file is SM-7LDIM.tar
## file prefix is SM-7LDIM
## zip_type is VIR
## uid is SM-7LDIM
## file is ESM5GEYY_P.fastq.gz
## file prefix is ESM5GEYY_P
## zip_type is P
## uid is ESM5GEYY_P
## file is ESM5GEYW_P.fastq.gz
## file prefix is ESM5GEYW_P
## zip_type is P
## uid is ESM5GEYW_P
## file is ESM5MEAB_P.fastq.gz
## file prefix is ESM5MEAB_P
## zip_type is P
## uid is ESM5MEAB_P
## file is ESM5MEA7_P.fastq.gz
## file prefix is ESM5MEA7_P
## zip_type is P
## uid is ESM5MEA7_P
## file is ESM5MEA9_P.fastq.gz
## file prefix is ESM5MEA9_P
## zip_type is P
## uid is ESM5MEA9_P
## file is ESM5MEDN.tar
## file prefix is ESM5MEDN
## zip_type is STD
## uid is ESM5MEDN
## file is ESM5MEDP_P.fastq.gz
## file prefix is ESM5MEDP_P
## zip_type is P
## uid is ESM5MEDP_P
## file is ESM5MEDK.tar
## file prefix is ESM5MEDK
## zip_type is STD
## uid is ESM5MEDK
## file is ESM5MEDD.tar
## file prefix is ESM5MEDD
## zip_type is STD
## uid is ESM5MEDD
## file is ESM5MEDF.tar
## file prefix is ESM5MEDF
## zip_type is STD
## uid is ESM5MEDF
## file is ESM718UH.tar
## file prefix is ESM718UH
## zip_type is STD
## uid is ESM718UH
## file is ESM7F5AK.tar
## file prefix is ESM7F5AK
## zip_type is STD
## uid is ESM7F5AK
## file is ESM7F5AM.tar
## file prefix is ESM7F5AM
## zip_type is STD
## uid is ESM7F5AM
## file is ESM7F5CB.tar
## file prefix is ESM7F5CB
## zip_type is STD
## uid is ESM7F5CB
## file is ESM7F5CD.tar
## file prefix is ESM7F5CD
## zip_type is STD
## uid is ESM7F5CD
## file is ESM7F5CF.tar
## file prefix is ESM7F5CF
## zip_type is STD
## uid is ESM7F5CF
## file is SM-6X9WB.tar
## file prefix is SM-6X9WB
## zip_type is VIR
## uid is SM-6X9WB
## file is SM-74Y6X.tar
## file prefix is SM-74Y6X
## zip_type is VIR
## uid is SM-74Y6X
## file is SM-9HC7V.tar
## file prefix is SM-9HC7V
## zip_type is VIR
## uid is SM-9HC7V
## file is ESM5MEB9_P.fastq.gz
## file prefix is ESM5MEB9_P
## zip_type is P
## uid is ESM5MEB9_P
## file is ESM5MEB7.tar
## file prefix is ESM5MEB7
## zip_type is STD
## uid is ESM5MEB7
## file is ESM5MED2.tar
## file prefix is ESM5MED2
## zip_type is STD
## uid is ESM5MED2
## file is ESM718V8.tar
## file prefix is ESM718V8
## zip_type is STD
## uid is ESM718V8
## file is ESM718V4.tar
## file prefix is ESM718V4
## zip_type is STD
## uid is ESM718V4
## file is ESM718TK.tar
## file prefix is ESM718TK
## zip_type is STD
## uid is ESM718TK
## file is ESM718TM.tar
## file prefix is ESM718TM
## zip_type is STD
## uid is ESM718TM
## file is ESM718TF.tar
## file prefix is ESM718TF
## zip_type is STD
## uid is ESM718TF
## file is ESM718T9.tar
## file prefix is ESM718T9
## zip_type is STD
## uid is ESM718T9
## file is ESM718T7.tar
## file prefix is ESM718T7
## zip_type is STD
## uid is ESM718T7
## file is ESM7F5C5.tar
## file prefix is ESM7F5C5
## zip_type is STD
## uid is ESM7F5C5
## file is ESM7F5C7.tar
## file prefix is ESM7F5C7
## zip_type is STD
## uid is ESM7F5C7
## file is ESM9IEP1.tar
## file prefix is ESM9IEP1
## zip_type is STD
## uid is ESM9IEP1
## file is SM-6ZUH7.tar
## file prefix is SM-6ZUH7
## zip_type is VIR
## uid is SM-6ZUH7
## file is SM-72XH2.tar
## file prefix is SM-72XH2
## zip_type is VIR
## uid is SM-72XH2
## file is SM-77M5L.tar
## file prefix is SM-77M5L
## zip_type is VIR
## uid is SM-77M5L
## file is SM-7E7H6.tar
## file prefix is SM-7E7H6
## zip_type is VIR
## uid is SM-7E7H6
## file is SM-7K1V4.tar
## file prefix is SM-7K1V4
## zip_type is VIR
## uid is SM-7K1V4
## file is SM-7PO8P.tar
## file prefix is SM-7PO8P
## zip_type is VIR
## uid is SM-7PO8P
## file is SM-9JG91.tar
## file prefix is SM-9JG91
## zip_type is VIR
## uid is SM-9JG91
## file is ESM7F5AE_P.fastq.gz
## file prefix is ESM7F5AE_P
## zip_type is P
## uid is ESM7F5AE_P
## file is ESM5MEC7_P.fastq.gz
## file prefix is ESM5MEC7_P
## zip_type is P
## uid is ESM5MEC7_P
## file is ESM718U9_P.fastq.gz
## file prefix is ESM718U9_P
## zip_type is P
## uid is ESM718U9_P
## file is HSM5FZBQ_P.fastq.gz
## file prefix is HSM5FZBQ_P
## zip_type is P
## uid is HSM5FZBQ_P
## file is HSM5FZBR_P.fastq.gz
## file prefix is HSM5FZBR_P
## zip_type is P
## uid is HSM5FZBR_P
## file is HSM5FZBP_P.fastq.gz
## file prefix is HSM5FZBP_P
## zip_type is P
## uid is HSM5FZBP_P
## file is HSM5MD7W_P.fastq.gz
## file prefix is HSM5MD7W_P
## zip_type is P
## uid is HSM5MD7W_P
## file is HSM5MD5H_P.fastq.gz
## file prefix is HSM5MD5H_P
## zip_type is P
## uid is HSM5MD5H_P
## file is HSM5MD5K_P.fastq.gz
## file prefix is HSM5MD5K_P
## zip_type is P
## uid is HSM5MD5K_P
## file is HSM5FZC2_P.fastq.gz
## file prefix is HSM5FZC2_P
## zip_type is P
## uid is HSM5FZC2_P
## file is HSM5FZBZ.tar
## file prefix is HSM5FZBZ
## zip_type is STD
## uid is HSM5FZBZ
## file is HSM5MD7J.tar
## file prefix is HSM5MD7J
## zip_type is STD
## uid is HSM5MD7J
## file is HSM5MD79.tar
## file prefix is HSM5MD79
## zip_type is STD
## uid is HSM5MD79
## file is HSM67VDZ.tar
## file prefix is HSM67VDZ
## zip_type is STD
## uid is HSM67VDZ
## file is HSM67VE4.tar
## file prefix is HSM67VE4
## zip_type is STD
## uid is HSM67VE4
## file is SM-6T7W7.tar
## file prefix is SM-6T7W7
## zip_type is VIR
## uid is SM-6T7W7
## file is SM-6ZEU9.tar
## file prefix is SM-6ZEU9
## zip_type is VIR
## uid is SM-6ZEU9
## file is SM-7DGUP.tar
## file prefix is SM-7DGUP
## zip_type is VIR
## uid is SM-7DGUP
## file is HSM5MD4U_P.fastq.gz
## file prefix is HSM5MD4U_P
## zip_type is P
## uid is HSM5MD4U_P
## file is HSM5MD4W_P.fastq.gz
## file prefix is HSM5MD4W_P
## zip_type is P
## uid is HSM5MD4W_P
## file is HSM5MD4Y.tar
## file prefix is HSM5MD4Y
## zip_type is STD
## uid is HSM5MD4Y
## file is HSM5MD53.tar
## file prefix is HSM5MD53
## zip_type is STD
## uid is HSM5MD53
## file is HSM5MD5P.tar
## file prefix is HSM5MD5P
## zip_type is STD
## uid is HSM5MD5P
## file is HSM6XRSG.tar
## file prefix is HSM6XRSG
## zip_type is STD
## uid is HSM6XRSG
## file is HSM6XRSI.tar
## file prefix is HSM6XRSI
## zip_type is STD
## uid is HSM6XRSI
## file is HSM67VDP.tar
## file prefix is HSM67VDP
## zip_type is STD
## uid is HSM67VDP
## file is HSM7CYZT.tar
## file prefix is HSM7CYZT
## zip_type is STD
## uid is HSM7CYZT
## file is HSM7CYZV.tar
## file prefix is HSM7CYZV
## zip_type is STD
## uid is HSM7CYZV
## file is HSM7CZ14.tar
## file prefix is HSM7CZ14
## zip_type is STD
## uid is HSM7CZ14
## file is SM-6ZEUH.tar
## file prefix is SM-6ZEUH
## zip_type is VIR
## uid is SM-6ZEUH
## file is SM-77M6B.tar
## file prefix is SM-77M6B
## zip_type is VIR
## uid is SM-77M6B
## file is SM-7AA27.tar
## file prefix is SM-7AA27
## zip_type is VIR
## uid is SM-7AA27
## file is SM-7HXJV.tar
## file prefix is SM-7HXJV
## zip_type is VIR
## uid is SM-7HXJV
## file is SM-9ISZZ.tar
## file prefix is SM-9ISZZ
## zip_type is VIR
## uid is SM-9ISZZ
## file is HSM5MD4B_P.fastq.gz
## file prefix is HSM5MD4B_P
## zip_type is P
## uid is HSM5MD4B_P
## file is HSM5MD4A_P.fastq.gz
## file prefix is HSM5MD4A_P
## zip_type is P
## uid is HSM5MD4A_P
## file is HSM5MD49_P.fastq.gz
## file prefix is HSM5MD49_P
## zip_type is P
## uid is HSM5MD49_P
## file is HSM5MD48.tar
## file prefix is HSM5MD48
## zip_type is STD
## uid is HSM5MD48
## file is HSM5MD7K.tar
## file prefix is HSM5MD7K
## zip_type is STD
## uid is HSM5MD7K
## file is HSM5MD7M.tar
## file prefix is HSM5MD7M
## zip_type is STD
## uid is HSM5MD7M
## file is HSM5MD7O.tar
## file prefix is HSM5MD7O
## zip_type is STD
## uid is HSM5MD7O
## file is HSM5MD7Q.tar
## file prefix is HSM5MD7Q
## zip_type is STD
## uid is HSM5MD7Q
## file is HSM5MD7S.tar
## file prefix is HSM5MD7S
## zip_type is STD
## uid is HSM5MD7S
## file is HSM5MD7U.tar
## file prefix is HSM5MD7U
## zip_type is STD
## uid is HSM5MD7U
## file is HSM67VEC.tar
## file prefix is HSM67VEC
## zip_type is STD
## uid is HSM67VEC
## file is HSM67VEE.tar
## file prefix is HSM67VEE
## zip_type is STD
## uid is HSM67VEE
## file is HSM67VEG.tar
## file prefix is HSM67VEG
## zip_type is STD
## uid is HSM67VEG
## file is HSM67VEI.tar
## file prefix is HSM67VEI
## zip_type is STD
## uid is HSM67VEI
## file is HSM67VEK.tar
## file prefix is HSM67VEK
## zip_type is STD
## uid is HSM67VEK
## file is HSM67VEM.tar
## file prefix is HSM67VEM
## zip_type is STD
## uid is HSM67VEM
## file is HSM67VEM_TR.tar
## file prefix is HSM67VEM_TR
## zip_type is TR
## uid is HSM67VEM_TR
## file is HSM7CYX2.tar
## file prefix is HSM7CYX2
## zip_type is STD
## uid is HSM7CYX2
## file is HSM7CYX4.tar
## file prefix is HSM7CYX4
## zip_type is STD
## uid is HSM7CYX4
## file is HSM7CYX6.tar
## file prefix is HSM7CYX6
## zip_type is STD
## uid is HSM7CYX6
## file is HSM7CYX8.tar
## file prefix is HSM7CYX8
## zip_type is STD
## uid is HSM7CYX8
## file is HSM7CYXA.tar
## file prefix is HSM7CYXA
## zip_type is STD
## uid is HSM7CYXA
## file is HSM7CYXC.tar
## file prefix is HSM7CYXC
## zip_type is STD
## uid is HSM7CYXC
## file is SM-6WJMT.tar
## file prefix is SM-6WJMT
## zip_type is VIR
## uid is SM-6WJMT
## file is SM-6XSLB.tar
## file prefix is SM-6XSLB
## zip_type is VIR
## uid is SM-6XSLB
## file is SM-6ZT3I.tar
## file prefix is SM-6ZT3I
## zip_type is VIR
## uid is SM-6ZT3I
## file is SM-72PGG.tar
## file prefix is SM-72PGG
## zip_type is VIR
## uid is SM-72PGG
## file is SM-75ZGU.tar
## file prefix is SM-75ZGU
## zip_type is VIR
## uid is SM-75ZGU
## file is SM-787IA.tar
## file prefix is SM-787IA
## zip_type is VIR
## uid is SM-787IA
## file is SM-7AMIE.tar
## file prefix is SM-7AMIE
## zip_type is VIR
## uid is SM-7AMIE
## file is SM-7BP4W.tar
## file prefix is SM-7BP4W
## zip_type is VIR
## uid is SM-7BP4W
## file is SM-7FK3L.tar
## file prefix is SM-7FK3L
## zip_type is VIR
## uid is SM-7FK3L
## file is SM-7GYIS.tar
## file prefix is SM-7GYIS
## zip_type is VIR
## uid is SM-7GYIS
## file is SM-7I1FT.tar
## file prefix is SM-7I1FT
## zip_type is VIR
## uid is SM-7I1FT
## file is SM-7K1V2.tar
## file prefix is SM-7K1V2
## zip_type is VIR
## uid is SM-7K1V2
## file is SM-7LDIG.tar
## file prefix is SM-7LDIG
## zip_type is VIR
## uid is SM-7LDIG
## file is SM-7NH5C.tar
## file prefix is SM-7NH5C
## zip_type is VIR
## uid is SM-7NH5C
## file is SM-7R3BC.tar
## file prefix is SM-7R3BC
## zip_type is VIR
## uid is SM-7R3BC
## file is SM-9KON3.tar
## file prefix is SM-9KON3
## zip_type is VIR
## uid is SM-9KON3
## file is HSM5MD7Z_P.fastq.gz
## file prefix is HSM5MD7Z_P
## zip_type is P
## uid is HSM5MD7Z_P
## file is HSM5MD4P_P.fastq.gz
## file prefix is HSM5MD4P_P
## zip_type is P
## uid is HSM5MD4P_P
## file is HSM5MD4O.tar
## file prefix is HSM5MD4O
## zip_type is STD
## uid is HSM5MD4O
## file is HSM5MD4N.tar
## file prefix is HSM5MD4N
## zip_type is STD
## uid is HSM5MD4N
## file is HSM6XRSN.tar
## file prefix is HSM6XRSN
## zip_type is STD
## uid is HSM6XRSN
## file is HSM6XRST.tar
## file prefix is HSM6XRST
## zip_type is STD
## uid is HSM6XRST
## file is HSM67VHQ.tar
## file prefix is HSM67VHQ
## zip_type is STD
## uid is HSM67VHQ
## file is HSM67VHS.tar
## file prefix is HSM67VHS
## zip_type is STD
## uid is HSM67VHS
## file is HSM67VHW.tar
## file prefix is HSM67VHW
## zip_type is STD
## uid is HSM67VHW
## file is HSM67VI1.tar
## file prefix is HSM67VI1
## zip_type is STD
## uid is HSM67VI1
## file is HSM7CYXQ.tar
## file prefix is HSM7CYXQ
## zip_type is STD
## uid is HSM7CYXQ
## file is HSM7CYXS.tar
## file prefix is HSM7CYXS
## zip_type is STD
## uid is HSM7CYXS
## file is SM-6Y2UQ.tar
## file prefix is SM-6Y2UQ
## zip_type is VIR
## uid is SM-6Y2UQ
## file is SM-72XGU.tar
## file prefix is SM-72XGU
## zip_type is VIR
## uid is SM-72XGU
## file is SM-7AA1X.tar
## file prefix is SM-7AA1X
## zip_type is VIR
## uid is SM-7AA1X
## file is SM-7GYI7.tar
## file prefix is SM-7GYI7
## zip_type is VIR
## uid is SM-7GYI7
## file is HSM5MD8A_P.fastq.gz
## file prefix is HSM5MD8A_P
## zip_type is P
## uid is HSM5MD8A_P
## file is HSM5MD57_P.fastq.gz
## file prefix is HSM5MD57_P
## zip_type is P
## uid is HSM5MD57_P
## file is HSM5MD59_P.fastq.gz
## file prefix is HSM5MD59_P
## zip_type is P
## uid is HSM5MD59_P
## file is HSM5MD5B_P.fastq.gz
## file prefix is HSM5MD5B_P
## zip_type is P
## uid is HSM5MD5B_P
## file is HSM5MD5B.tar
## file prefix is HSM5MD5B
## zip_type is STD
## uid is HSM5MD5B
## file is HSM5MD5D_P.fastq.gz
## file prefix is HSM5MD5D_P
## zip_type is P
## uid is HSM5MD5D_P
## file is HSM5MD5D.tar
## file prefix is HSM5MD5D
## zip_type is STD
## uid is HSM5MD5D
## file is HSM5MD5F_P.fastq.gz
## file prefix is HSM5MD5F_P
## zip_type is P
## uid is HSM5MD5F_P
## file is HSM6XRSX_P.fastq.gz
## file prefix is HSM6XRSX_P
## zip_type is P
## uid is HSM6XRSX_P
## file is HSM6XRSX.tar
## file prefix is HSM6XRSX
## zip_type is STD
## uid is HSM6XRSX
## file is HSM6XRT6_P.fastq.gz
## file prefix is HSM6XRT6_P
## zip_type is P
## uid is HSM6XRT6_P
## file is HSM6XRSZ_P.fastq.gz
## file prefix is HSM6XRSZ_P
## zip_type is P
## uid is HSM6XRSZ_P
## file is HSM6XRT4_P.fastq.gz
## file prefix is HSM6XRT4_P
## zip_type is P
## uid is HSM6XRT4_P
## file is HSM6XRT2_P.fastq.gz
## file prefix is HSM6XRT2_P
## zip_type is P
## uid is HSM6XRT2_P
## file is HSM67VFX_P.fastq.gz
## file prefix is HSM67VFX_P
## zip_type is P
## uid is HSM67VFX_P
## file is HSM67VFX.tar
## file prefix is HSM67VFX
## zip_type is STD
## uid is HSM67VFX
## file is HSM67VFZ.tar
## file prefix is HSM67VFZ
## zip_type is STD
## uid is HSM67VFZ
## file is HSM67VG2_P.fastq.gz
## file prefix is HSM67VG2_P
## zip_type is P
## uid is HSM67VG2_P
## file is HSM67VG6_P.fastq.gz
## file prefix is HSM67VG6_P
## zip_type is P
## uid is HSM67VG6_P
## file is HSM67VG8.tar
## file prefix is HSM67VG8
## zip_type is STD
## uid is HSM67VG8
## file is CSM7CZ2F_P.fastq.gz
## file prefix is CSM7CZ2F_P
## zip_type is P
## uid is CSM7CZ2F_P
## file is HSM7CZ2H.tar
## file prefix is HSM7CZ2H
## zip_type is STD
## uid is HSM7CZ2H
## file is HSM7CZ2J_P.fastq.gz
## file prefix is HSM7CZ2J_P
## zip_type is P
## uid is HSM7CZ2J_P
## file is HSM7CZ2L_P.fastq.gz
## file prefix is HSM7CZ2L_P
## zip_type is P
## uid is HSM7CZ2L_P
## file is SM-6UG73.tar
## file prefix is SM-6UG73
## zip_type is VIR
## uid is SM-6UG73
## file is SM-6X9W7.tar
## file prefix is SM-6X9W7
## zip_type is VIR
## uid is SM-6X9W7
## file is SM-6ZJXW.tar
## file prefix is SM-6ZJXW
## zip_type is VIR
## uid is SM-6ZJXW
## file is SM-7EDAD.tar
## file prefix is SM-7EDAD
## zip_type is VIR
## uid is SM-7EDAD
## file is SM-7L29I.tar
## file prefix is SM-7L29I
## zip_type is VIR
## uid is SM-7L29I
## file is HSM5MD82_P.fastq.gz
## file prefix is HSM5MD82_P
## zip_type is P
## uid is HSM5MD82_P
## file is HSM5MD8P_P.fastq.gz
## file prefix is HSM5MD8P_P
## zip_type is P
## uid is HSM5MD8P_P
## file is HSM5MD8N_P.fastq.gz
## file prefix is HSM5MD8N_P
## zip_type is P
## uid is HSM5MD8N_P
## file is HSM5MD8L_P.fastq.gz
## file prefix is HSM5MD8L_P
## zip_type is P
## uid is HSM5MD8L_P
## file is HSM5MD8B_P.fastq.gz
## file prefix is HSM5MD8B_P
## zip_type is P
## uid is HSM5MD8B_P
## file is HSM5MD8D_P.fastq.gz
## file prefix is HSM5MD8D_P
## zip_type is P
## uid is HSM5MD8D_P
## file is HSM6XRT8.tar
## file prefix is HSM6XRT8
## zip_type is STD
## uid is HSM6XRT8
## file is HSM6XRTA_P.fastq.gz
## file prefix is HSM6XRTA_P
## zip_type is P
## uid is HSM6XRTA_P
## file is HSM6XRTC_P.fastq.gz
## file prefix is HSM6XRTC_P
## zip_type is P
## uid is HSM6XRTC_P
## file is HSM6XRTE_P.fastq.gz
## file prefix is HSM6XRTE_P
## zip_type is P
## uid is HSM6XRTE_P
## file is HSM6XRTG_P.fastq.gz
## file prefix is HSM6XRTG_P
## zip_type is P
## uid is HSM6XRTG_P
## file is HSM6XRTG.tar
## file prefix is HSM6XRTG
## zip_type is STD
## uid is HSM6XRTG
## file is HSM67VGA_P.fastq.gz
## file prefix is HSM67VGA_P
## zip_type is P
## uid is HSM67VGA_P
## file is HSM67VGA.tar
## file prefix is HSM67VGA
## zip_type is STD
## uid is HSM67VGA
## file is HSM67VGC.tar
## file prefix is HSM67VGC
## zip_type is STD
## uid is HSM67VGC
## file is HSM67VGE_P.fastq.gz
## file prefix is HSM67VGE_P
## zip_type is P
## uid is HSM67VGE_P
## file is HSM67VGG.tar
## file prefix is HSM67VGG
## zip_type is STD
## uid is HSM67VGG
## file is HSM67VGI_P.fastq.gz
## file prefix is HSM67VGI_P
## zip_type is P
## uid is HSM67VGI_P
## file is HSM67VGK.tar
## file prefix is HSM67VGK
## zip_type is STD
## uid is HSM67VGK
## file is HSM67VGK_TR.tar
## file prefix is HSM67VGK_TR
## zip_type is TR
## uid is HSM67VGK_TR
## file is HSM7CYYR_P.fastq.gz
## file prefix is HSM7CYYR_P
## zip_type is P
## uid is HSM7CYYR_P
## file is HSM7CYYV_P.fastq.gz
## file prefix is HSM7CYYV_P
## zip_type is P
## uid is HSM7CYYV_P
## file is SM-6ZKRY.tar
## file prefix is SM-6ZKRY
## zip_type is VIR
## uid is SM-6ZKRY
## file is SM-7CRWM.tar
## file prefix is SM-7CRWM
## zip_type is VIR
## uid is SM-7CRWM
## file is SM-7F1RM.tar
## file prefix is SM-7F1RM
## zip_type is VIR
## uid is SM-7F1RM
## file is SM-7IWFK.tar
## file prefix is SM-7IWFK
## zip_type is VIR
## uid is SM-7IWFK
## file is SM-7M8RE.tar
## file prefix is SM-7M8RE
## zip_type is VIR
## uid is SM-7M8RE
## file is HSM5MD87_P.fastq.gz
## file prefix is HSM5MD87_P
## zip_type is P
## uid is HSM5MD87_P
## file is HSM5MD47_P.fastq.gz
## file prefix is HSM5MD47_P
## zip_type is P
## uid is HSM5MD47_P
## file is HSM5MD44_P.fastq.gz
## file prefix is HSM5MD44_P
## zip_type is P
## uid is HSM5MD44_P
## file is HSM5MD43_P.fastq.gz
## file prefix is HSM5MD43_P
## zip_type is P
## uid is HSM5MD43_P
## file is HSM5MD43.tar
## file prefix is HSM5MD43
## zip_type is STD
## uid is HSM5MD43
## file is HSM5MD3Y.tar
## file prefix is HSM5MD3Y
## zip_type is STD
## uid is HSM5MD3Y
## file is HSM5MD41.tar
## file prefix is HSM5MD41
## zip_type is STD
## uid is HSM5MD41
## file is HSM6XRTM.tar
## file prefix is HSM6XRTM
## zip_type is STD
## uid is HSM6XRTM
## file is HSM6XRTQ.tar
## file prefix is HSM6XRTQ
## zip_type is STD
## uid is HSM6XRTQ
## file is HSM6XRTO.tar
## file prefix is HSM6XRTO
## zip_type is STD
## uid is HSM6XRTO
## file is HSM6XRTS.tar
## file prefix is HSM6XRTS
## zip_type is STD
## uid is HSM6XRTS
## file is HSM67VGY.tar
## file prefix is HSM67VGY
## zip_type is STD
## uid is HSM67VGY
## file is HSM67VH1.tar
## file prefix is HSM67VH1
## zip_type is STD
## uid is HSM67VH1
## file is SM-6XJUD.tar
## file prefix is SM-6XJUD
## zip_type is VIR
## uid is SM-6XJUD
## file is SM-6ZT3E.tar
## file prefix is SM-6ZT3E
## zip_type is VIR
## uid is SM-6ZT3E
## file is SM-73JXF.tar
## file prefix is SM-73JXF
## zip_type is VIR
## uid is SM-73JXF
## file is SM-7CP2W.tar
## file prefix is SM-7CP2W
## zip_type is VIR
## uid is SM-7CP2W
## file is SM-7FHYC.tar
## file prefix is SM-7FHYC
## zip_type is VIR
## uid is SM-7FHYC
## file is SM-7PO8X.tar
## file prefix is SM-7PO8X
## zip_type is VIR
## uid is SM-7PO8X
## file is SM-9T55I.tar
## file prefix is SM-9T55I
## zip_type is VIR
## uid is SM-9T55I
## file is HSM5MD8J_P.fastq.gz
## file prefix is HSM5MD8J_P
## zip_type is P
## uid is HSM5MD8J_P
## file is HSM5MD8F_P.fastq.gz
## file prefix is HSM5MD8F_P
## zip_type is P
## uid is HSM5MD8F_P
## file is HSM5MD8H_P.fastq.gz
## file prefix is HSM5MD8H_P
## zip_type is P
## uid is HSM5MD8H_P
## file is HSM5MD6E_P.fastq.gz
## file prefix is HSM5MD6E_P
## zip_type is P
## uid is HSM5MD6E_P
## file is HSM5MD6I_P.fastq.gz
## file prefix is HSM5MD6I_P
## zip_type is P
## uid is HSM5MD6I_P
## file is HSM5MD6I.tar
## file prefix is HSM5MD6I
## zip_type is STD
## uid is HSM5MD6I
## file is HSM5MD6K_P.fastq.gz
## file prefix is HSM5MD6K_P
## zip_type is P
## uid is HSM5MD6K_P
## file is HSM5MD6K.tar
## file prefix is HSM5MD6K
## zip_type is STD
## uid is HSM5MD6K
## file is HSM5MD6M.tar
## file prefix is HSM5MD6M
## zip_type is STD
## uid is HSM5MD6M
## file is HSM6XRRJ.tar
## file prefix is HSM6XRRJ
## zip_type is STD
## uid is HSM6XRRJ
## file is HSM6XRVA.tar
## file prefix is HSM6XRVA
## zip_type is STD
## uid is HSM6XRVA
## file is HSM6XRVC.tar
## file prefix is HSM6XRVC
## zip_type is STD
## uid is HSM6XRVC
## file is HSM6XRVC_TR.tar
## file prefix is HSM6XRVC_TR
## zip_type is TR
## uid is HSM6XRVC_TR
## file is HSM6XRVK.tar
## file prefix is HSM6XRVK
## zip_type is STD
## uid is HSM6XRVK
## file is HSM7J4QZ.tar
## file prefix is HSM7J4QZ
## zip_type is STD
## uid is HSM7J4QZ
## file is HSM7J4R2.tar
## file prefix is HSM7J4R2
## zip_type is STD
## uid is HSM7J4R2
## file is SM-6ZEUX.tar
## file prefix is SM-6ZEUX
## zip_type is VIR
## uid is SM-6ZEUX
## file is SM-6ZUH3.tar
## file prefix is SM-6ZUH3
## zip_type is VIR
## uid is SM-6ZUH3
## file is SM-7AMIB.tar
## file prefix is SM-7AMIB
## zip_type is VIR
## uid is SM-7AMIB
## file is SM-7I1FP.tar
## file prefix is SM-7I1FP
## zip_type is VIR
## uid is SM-7I1FP
## file is SM-7NH6B.tar
## file prefix is SM-7NH6B
## zip_type is VIR
## uid is SM-7NH6B
## file is SM-9R78E.tar
## file prefix is SM-9R78E
## zip_type is VIR
## uid is SM-9R78E
## file is HSM5MD3L_P.fastq.gz
## file prefix is HSM5MD3L_P
## zip_type is P
## uid is HSM5MD3L_P
## file is HSM5MD66.tar
## file prefix is HSM5MD66
## zip_type is STD
## uid is HSM5MD66
## file is HSM5MD6A.tar
## file prefix is HSM5MD6A
## zip_type is STD
## uid is HSM5MD6A
## file is HSM5MD6A_TR.tar
## file prefix is HSM5MD6A_TR
## zip_type is TR
## uid is HSM5MD6A_TR
## file is HSM5MD6C.tar
## file prefix is HSM5MD6C
## zip_type is STD
## uid is HSM5MD6C
## file is HSM6XRRV.tar
## file prefix is HSM6XRRV
## zip_type is STD
## uid is HSM6XRRV
## file is HSM6XRS2.tar
## file prefix is HSM6XRS2
## zip_type is STD
## uid is HSM6XRS2
## file is HSM6XRVM.tar
## file prefix is HSM6XRVM
## zip_type is STD
## uid is HSM6XRVM
## file is HSM6XRVO.tar
## file prefix is HSM6XRVO
## zip_type is STD
## uid is HSM6XRVO
## file is HSM6XRVU.tar
## file prefix is HSM6XRVU
## zip_type is STD
## uid is HSM6XRVU
## file is HSM6XRVW.tar
## file prefix is HSM6XRVW
## zip_type is STD
## uid is HSM6XRVW
## file is HSM7J4QT.tar
## file prefix is HSM7J4QT
## zip_type is STD
## uid is HSM7J4QT
## file is SM-6XSLF.tar
## file prefix is SM-6XSLF
## zip_type is VIR
## uid is SM-6XSLF
## file is SM-6ZT2W.tar
## file prefix is SM-6ZT2W
## zip_type is VIR
## uid is SM-6ZT2W
## file is SM-7AMII.tar
## file prefix is SM-7AMII
## zip_type is VIR
## uid is SM-7AMII
## file is SM-7GYIK.tar
## file prefix is SM-7GYIK
## zip_type is VIR
## uid is SM-7GYIK
## file is SM-7IKZG.tar
## file prefix is SM-7IKZG
## zip_type is VIR
## uid is SM-7IKZG
## file is SM-9KOMQ.tar
## file prefix is SM-9KOMQ
## zip_type is VIR
## uid is SM-9KOMQ
## file is SM-9ZJIJ.tar
## file prefix is SM-9ZJIJ
## zip_type is VIR
## uid is SM-9ZJIJ
## file is HSM5MD5X_P.fastq.gz
## file prefix is HSM5MD5X_P
## zip_type is P
## uid is HSM5MD5X_P
## file is HSM5MD62.tar
## file prefix is HSM5MD62
## zip_type is STD
## uid is HSM5MD62
## file is HSM5MD6Y.tar
## file prefix is HSM5MD6Y
## zip_type is STD
## uid is HSM5MD6Y
## file is HSM5MD71.tar
## file prefix is HSM5MD71
## zip_type is STD
## uid is HSM5MD71
## file is HSM5MD73.tar
## file prefix is HSM5MD73
## zip_type is STD
## uid is HSM5MD73
## file is HSM5MD75.tar
## file prefix is HSM5MD75
## zip_type is STD
## uid is HSM5MD75
## file is HSM6XRS4.tar
## file prefix is HSM6XRS4
## zip_type is STD
## uid is HSM6XRS4
## file is HSM6XRS6.tar
## file prefix is HSM6XRS6
## zip_type is STD
## uid is HSM6XRS6
## file is HSM6XRS8.tar
## file prefix is HSM6XRS8
## zip_type is STD
## uid is HSM6XRS8
## file is HSM6XRSE.tar
## file prefix is HSM6XRSE
## zip_type is STD
## uid is HSM6XRSE
## file is HSM7CYZ5.tar
## file prefix is HSM7CYZ5
## zip_type is STD
## uid is HSM7CYZ5
## file is HSM7CYZ7.tar
## file prefix is HSM7CYZ7
## zip_type is STD
## uid is HSM7CYZ7
## file is HSM7CYZ9.tar
## file prefix is HSM7CYZ9
## zip_type is STD
## uid is HSM7CYZ9
## file is HSM7CYZB.tar
## file prefix is HSM7CYZB
## zip_type is STD
## uid is HSM7CYZB
## file is HSM7CYZD.tar
## file prefix is HSM7CYZD
## zip_type is STD
## uid is HSM7CYZD
## file is HSM7CYZF.tar
## file prefix is HSM7CYZF
## zip_type is STD
## uid is HSM7CYZF
## file is HSM7J4QB.tar
## file prefix is HSM7J4QB
## zip_type is STD
## uid is HSM7J4QB
## file is HSM7J4QD.tar
## file prefix is HSM7J4QD
## zip_type is STD
## uid is HSM7J4QD
## file is HSM7J4QF.tar
## file prefix is HSM7J4QF
## zip_type is STD
## uid is HSM7J4QF
## file is HSM7J4QH.tar
## file prefix is HSM7J4QH
## zip_type is STD
## uid is HSM7J4QH
## file is HSM7J4QJ.tar
## file prefix is HSM7J4QJ
## zip_type is STD
## uid is HSM7J4QJ
## file is HSM7J4QL.tar
## file prefix is HSM7J4QL
## zip_type is STD
## uid is HSM7J4QL
## file is SM-6X9W3.tar
## file prefix is SM-6X9W3
## zip_type is VIR
## uid is SM-6X9W3
## file is SM-6YBJJ.tar
## file prefix is SM-6YBJJ
## zip_type is VIR
## uid is SM-6YBJJ
## file is SM-727GR.tar
## file prefix is SM-727GR
## zip_type is VIR
## uid is SM-727GR
## file is SM-76CA2.tar
## file prefix is SM-76CA2
## zip_type is VIR
## uid is SM-76CA2
## file is SM-7IKZM.tar
## file prefix is SM-7IKZM
## zip_type is VIR
## uid is SM-7IKZM
## file is SM-9P5G8.tar
## file prefix is SM-9P5G8
## zip_type is VIR
## uid is SM-9P5G8
## file is SM-9UCZE.tar
## file prefix is SM-9UCZE
## zip_type is VIR
## uid is SM-9UCZE
## file is HSM5MD5Z_P.fastq.gz
## file prefix is HSM5MD5Z_P
## zip_type is P
## uid is HSM5MD5Z_P
## file is HSM5MD6O_P.fastq.gz
## file prefix is HSM5MD6O_P
## zip_type is P
## uid is HSM5MD6O_P
## file is HSM5MD6Q.tar
## file prefix is HSM5MD6Q
## zip_type is STD
## uid is HSM5MD6Q
## file is HSM5MD6S_P.fastq.gz
## file prefix is HSM5MD6S_P
## zip_type is P
## uid is HSM5MD6S_P
## file is HSM5MD6U_P.fastq.gz
## file prefix is HSM5MD6U_P
## zip_type is P
## uid is HSM5MD6U_P
## file is HSM5MD6W.tar
## file prefix is HSM5MD6W
## zip_type is STD
## uid is HSM5MD6W
## file is HSM67VF3.tar
## file prefix is HSM67VF3
## zip_type is STD
## uid is HSM67VF3
## file is HSM7CYXE.tar
## file prefix is HSM7CYXE
## zip_type is STD
## uid is HSM7CYXE
## file is HSM7CYXG.tar
## file prefix is HSM7CYXG
## zip_type is STD
## uid is HSM7CYXG
## file is HSM7CYXI.tar
## file prefix is HSM7CYXI
## zip_type is STD
## uid is HSM7CYXI
## file is HSM7CYXO.tar
## file prefix is HSM7CYXO
## zip_type is STD
## uid is HSM7CYXO
## file is HSM7J4I3.tar
## file prefix is HSM7J4I3
## zip_type is STD
## uid is HSM7J4I3
## file is HSM7J4KO.tar
## file prefix is HSM7J4KO
## zip_type is STD
## uid is HSM7J4KO
## file is HSM7J4KQ.tar
## file prefix is HSM7J4KQ
## zip_type is STD
## uid is HSM7J4KQ
## file is SM-791BC.tar
## file prefix is SM-791BC
## zip_type is VIR
## uid is SM-791BC
## file is SM-7ED9S.tar
## file prefix is SM-7ED9S
## zip_type is VIR
## uid is SM-7ED9S
## file is SM-9P5FI.tar
## file prefix is SM-9P5FI
## zip_type is VIR
## uid is SM-9P5FI
## file is SM-A2HYL.tar
## file prefix is SM-A2HYL
## zip_type is VIR
## uid is SM-A2HYL
## file is HSM6XRQB_P.fastq.gz
## file prefix is HSM6XRQB_P
## zip_type is P
## uid is HSM6XRQB_P
## file is HSM6XRQB.tar
## file prefix is HSM6XRQB
## zip_type is STD
## uid is HSM6XRQB
## file is HSM67VHK.tar
## file prefix is HSM67VHK
## zip_type is STD
## uid is HSM67VHK
## file is HSM6XRQI.tar
## file prefix is HSM6XRQI
## zip_type is STD
## uid is HSM6XRQI
## file is HSM6XRQK.tar
## file prefix is HSM6XRQK
## zip_type is STD
## uid is HSM6XRQK
## file is HSM6XRQM.tar
## file prefix is HSM6XRQM
## zip_type is STD
## uid is HSM6XRQM
## file is HSM6XRQO.tar
## file prefix is HSM6XRQO
## zip_type is STD
## uid is HSM6XRQO
## file is HSM67VF9.tar
## file prefix is HSM67VF9
## zip_type is STD
## uid is HSM67VF9
## file is HSM67VFD.tar
## file prefix is HSM67VFD
## zip_type is STD
## uid is HSM67VFD
## file is HSM67VFF.tar
## file prefix is HSM67VFF
## zip_type is STD
## uid is HSM67VFF
## file is HSM67VFH.tar
## file prefix is HSM67VFH
## zip_type is STD
## uid is HSM67VFH
## file is HSM67VFJ.tar
## file prefix is HSM67VFJ
## zip_type is STD
## uid is HSM67VFJ
## file is HSM7CYY3.tar
## file prefix is HSM7CYY3
## zip_type is STD
## uid is HSM7CYY3
## file is HSM7CYY5.tar
## file prefix is HSM7CYY5
## zip_type is STD
## uid is HSM7CYY5
## file is HSM7CYY7.tar
## file prefix is HSM7CYY7
## zip_type is STD
## uid is HSM7CYY7
## file is HSM7CYYD.tar
## file prefix is HSM7CYYD
## zip_type is STD
## uid is HSM7CYYD
## file is HSM7CYY9.tar
## file prefix is HSM7CYY9
## zip_type is STD
## uid is HSM7CYY9
## file is HSM7CYYB.tar
## file prefix is HSM7CYYB
## zip_type is STD
## uid is HSM7CYYB
## file is SM-6ZKY4.tar
## file prefix is SM-6ZKY4
## zip_type is VIR
## uid is SM-6ZKY4
## file is SM-73JXJ.tar
## file prefix is SM-73JXJ
## zip_type is VIR
## uid is SM-73JXJ
## file is SM-76CAA.tar
## file prefix is SM-76CAA
## zip_type is VIR
## uid is SM-76CAA
## file is SM-7AK5W.tar
## file prefix is SM-7AK5W
## zip_type is VIR
## uid is SM-7AK5W
## file is SM-7DN2Z.tar
## file prefix is SM-7DN2Z
## zip_type is VIR
## uid is SM-7DN2Z
## file is SM-7HXJZ.tar
## file prefix is SM-7HXJZ
## zip_type is VIR
## uid is SM-7HXJZ
## file is SM-7NSOJ.tar
## file prefix is SM-7NSOJ
## zip_type is VIR
## uid is SM-7NSOJ
## file is SM-9KOPH.tar
## file prefix is SM-9KOPH
## zip_type is VIR
## uid is SM-9KOPH
## file is SM-9QMNI.tar
## file prefix is SM-9QMNI
## zip_type is VIR
## uid is SM-9QMNI
## file is SM-9RYIJ.tar
## file prefix is SM-9RYIJ
## zip_type is VIR
## uid is SM-9RYIJ
## file is SM-9U261.tar
## file prefix is SM-9U261
## zip_type is VIR
## uid is SM-9U261
## file is HSM6XRQC_P.fastq.gz
## file prefix is HSM6XRQC_P
## zip_type is P
## uid is HSM6XRQC_P
## file is HSM6XRQS.tar
## file prefix is HSM6XRQS
## zip_type is STD
## uid is HSM6XRQS
## file is HSM6XRQU.tar
## file prefix is HSM6XRQU
## zip_type is STD
## uid is HSM6XRQU
## file is HSM6XRQW.tar
## file prefix is HSM6XRQW
## zip_type is STD
## uid is HSM6XRQW
## file is HSM6XRQY.tar
## file prefix is HSM6XRQY
## zip_type is STD
## uid is HSM6XRQY
## file is HSM67VFR.tar
## file prefix is HSM67VFR
## zip_type is STD
## uid is HSM67VFR
## file is HSM7CYYF.tar
## file prefix is HSM7CYYF
## zip_type is STD
## uid is HSM7CYYF
## file is HSM7CYYH.tar
## file prefix is HSM7CYYH
## zip_type is STD
## uid is HSM7CYYH
## file is HSM7CYYP.tar
## file prefix is HSM7CYYP
## zip_type is STD
## uid is HSM7CYYP
## file is HSM7J4LD.tar
## file prefix is HSM7J4LD
## zip_type is STD
## uid is HSM7J4LD
## file is HSM7J4LF.tar
## file prefix is HSM7J4LF
## zip_type is STD
## uid is HSM7J4LF
## file is HSM7J4LH.tar
## file prefix is HSM7J4LH
## zip_type is STD
## uid is HSM7J4LH
## file is HSM7J4LN.tar
## file prefix is HSM7J4LN
## zip_type is STD
## uid is HSM7J4LN
## file is SM-76ENQ.tar
## file prefix is SM-76ENQ
## zip_type is VIR
## uid is SM-76ENQ
## file is SM-791B8.tar
## file prefix is SM-791B8
## zip_type is VIR
## uid is SM-791B8
## file is SM-7GYIC.tar
## file prefix is SM-7GYIC
## zip_type is VIR
## uid is SM-7GYIC
## file is SM-9RTV9.tar
## file prefix is SM-9RTV9
## zip_type is VIR
## uid is SM-9RTV9
## file is SM-9WOB7.tar
## file prefix is SM-9WOB7
## zip_type is VIR
## uid is SM-9WOB7
## file is SM-A373Q.tar
## file prefix is SM-A373Q
## zip_type is VIR
## uid is SM-A373Q
## file is HSM6XRQE_P.fastq.gz
## file prefix is HSM6XRQE_P
## zip_type is P
## uid is HSM6XRQE_P
## file is HSM6XRR5.tar
## file prefix is HSM6XRR5
## zip_type is STD
## uid is HSM6XRR5
## file is HSM6XRR7.tar
## file prefix is HSM6XRR7
## zip_type is STD
## uid is HSM6XRR7
## file is HSM6XRR9.tar
## file prefix is HSM6XRR9
## zip_type is STD
## uid is HSM6XRR9
## file is HSM6XRRB.tar
## file prefix is HSM6XRRB
## zip_type is STD
## uid is HSM6XRRB
## file is HSM6XRRD.tar
## file prefix is HSM6XRRD
## zip_type is STD
## uid is HSM6XRRD
## file is HSM6XRUX.tar
## file prefix is HSM6XRUX
## zip_type is STD
## uid is HSM6XRUX
## file is HSM6XRUZ.tar
## file prefix is HSM6XRUZ
## zip_type is STD
## uid is HSM6XRUZ
## file is HSM6XRV2.tar
## file prefix is HSM6XRV2
## zip_type is STD
## uid is HSM6XRV2
## file is HSM6XRV4.tar
## file prefix is HSM6XRV4
## zip_type is STD
## uid is HSM6XRV4
## file is HSM6XRV6.tar
## file prefix is HSM6XRV6
## zip_type is STD
## uid is HSM6XRV6
## file is HSM6XRV8.tar
## file prefix is HSM6XRV8
## zip_type is STD
## uid is HSM6XRV8
## file is HSM7J4PC.tar
## file prefix is HSM7J4PC
## zip_type is STD
## uid is HSM7J4PC
## file is HSM7J4PE.tar
## file prefix is HSM7J4PE
## zip_type is STD
## uid is HSM7J4PE
## file is HSM7J4PG.tar
## file prefix is HSM7J4PG
## zip_type is STD
## uid is HSM7J4PG
## file is HSM7J4PI.tar
## file prefix is HSM7J4PI
## zip_type is STD
## uid is HSM7J4PI
## file is HSM7J4PK.tar
## file prefix is HSM7J4PK
## zip_type is STD
## uid is HSM7J4PK
## file is HSM7J4PM.tar
## file prefix is HSM7J4PM
## zip_type is STD
## uid is HSM7J4PM
## file is HSM7J4M4.tar
## file prefix is HSM7J4M4
## zip_type is STD
## uid is HSM7J4M4
## file is HSM7J4M6.tar
## file prefix is HSM7J4M6
## zip_type is STD
## uid is HSM7J4M6
## file is HSM7J4M8.tar
## file prefix is HSM7J4M8
## zip_type is STD
## uid is HSM7J4M8
## file is HSM7J4MA.tar
## file prefix is HSM7J4MA
## zip_type is STD
## uid is HSM7J4MA
## file is HSM7J4MC.tar
## file prefix is HSM7J4MC
## zip_type is STD
## uid is HSM7J4MC
## file is HSM7J4ME.tar
## file prefix is HSM7J4ME
## zip_type is STD
## uid is HSM7J4ME
## file is SM-73BNC.tar
## file prefix is SM-73BNC
## zip_type is VIR
## uid is SM-73BNC
## file is SM-76C9Y.tar
## file prefix is SM-76C9Y
## zip_type is VIR
## uid is SM-76C9Y
## file is SM-7BP52.tar
## file prefix is SM-7BP52
## zip_type is VIR
## uid is SM-7BP52
## file is SM-7IKZL.tar
## file prefix is SM-7IKZL
## zip_type is VIR
## uid is SM-7IKZL
## file is SM-7KPU3.tar
## file prefix is SM-7KPU3
## zip_type is VIR
## uid is SM-7KPU3
## file is SM-7MCUM.tar
## file prefix is SM-7MCUM
## zip_type is VIR
## uid is SM-7MCUM
## file is SM-9NK6G.tar
## file prefix is SM-9NK6G
## zip_type is VIR
## uid is SM-9NK6G
## file is SM-9U26B.tar
## file prefix is SM-9U26B
## zip_type is VIR
## uid is SM-9U26B
## file is SM-9Y7DS.tar
## file prefix is SM-9Y7DS
## zip_type is VIR
## uid is SM-9Y7DS
## file is SM-A1ZWQ.tar
## file prefix is SM-A1ZWQ
## zip_type is VIR
## uid is SM-A1ZWQ
## file is SM-A4GBD.tar
## file prefix is SM-A4GBD
## zip_type is VIR
## uid is SM-A4GBD
## file is HSM6XRR3.tar
## file prefix is HSM6XRR3
## zip_type is STD
## uid is HSM6XRR3
## file is HSM67VCX_P.fastq.gz
## file prefix is HSM67VCX_P
## zip_type is P
## uid is HSM67VCX_P
## file is HSM67VCZ.tar
## file prefix is HSM67VCZ
## zip_type is STD
## uid is HSM67VCZ
## file is HSM67VD2.tar
## file prefix is HSM67VD2
## zip_type is STD
## uid is HSM67VD2
## file is HSM67VD4.tar
## file prefix is HSM67VD4
## zip_type is STD
## uid is HSM67VD4
## file is HSM67VD6.tar
## file prefix is HSM67VD6
## zip_type is STD
## uid is HSM67VD6
## file is HSM67VHB.tar
## file prefix is HSM67VHB
## zip_type is STD
## uid is HSM67VHB
## file is HSM67VHD.tar
## file prefix is HSM67VHD
## zip_type is STD
## uid is HSM67VHD
## file is HSM67VHF.tar
## file prefix is HSM67VHF
## zip_type is STD
## uid is HSM67VHF
## file is HSM67VHH.tar
## file prefix is HSM67VHH
## zip_type is STD
## uid is HSM67VHH
## file is HSM67VHJ.tar
## file prefix is HSM67VHJ
## zip_type is STD
## uid is HSM67VHJ
## file is HSM6XRUV.tar
## file prefix is HSM6XRUV
## zip_type is STD
## uid is HSM6XRUV
## file is HSM7J4PO.tar
## file prefix is HSM7J4PO
## zip_type is STD
## uid is HSM7J4PO
## file is HSM7J4PQ.tar
## file prefix is HSM7J4PQ
## zip_type is STD
## uid is HSM7J4PQ
## file is HSM7J4PS.tar
## file prefix is HSM7J4PS
## zip_type is STD
## uid is HSM7J4PS
## file is HSM7J4PU.tar
## file prefix is HSM7J4PU
## zip_type is STD
## uid is HSM7J4PU
## file is HSM7J4PW.tar
## file prefix is HSM7J4PW
## zip_type is STD
## uid is HSM7J4PW
## file is HSM7J4PY.tar
## file prefix is HSM7J4PY
## zip_type is STD
## uid is HSM7J4PY
## file is HSM7J4IO.tar
## file prefix is HSM7J4IO
## zip_type is STD
## uid is HSM7J4IO
## file is HSM7J4IP.tar
## file prefix is HSM7J4IP
## zip_type is STD
## uid is HSM7J4IP
## file is HSM7J4IQ.tar
## file prefix is HSM7J4IQ
## zip_type is STD
## uid is HSM7J4IQ
## file is HSM7J4IR.tar
## file prefix is HSM7J4IR
## zip_type is STD
## uid is HSM7J4IR
## file is HSM7J4IS.tar
## file prefix is HSM7J4IS
## zip_type is STD
## uid is HSM7J4IS
## file is SM-77FY2.tar
## file prefix is SM-77FY2
## zip_type is VIR
## uid is SM-77FY2
## file is SM-7AA1P.tar
## file prefix is SM-7AA1P
## zip_type is VIR
## uid is SM-7AA1P
## file is SM-7C83T.tar
## file prefix is SM-7C83T
## zip_type is VIR
## uid is SM-7C83T
## file is SM-7DGV2.tar
## file prefix is SM-7DGV2
## zip_type is VIR
## uid is SM-7DGV2
## file is SM-7EWSV.tar
## file prefix is SM-7EWSV
## zip_type is VIR
## uid is SM-7EWSV
## file is SM-7FY2T.tar
## file prefix is SM-7FY2T
## zip_type is VIR
## uid is SM-7FY2T
## file is SM-7I1GZ.tar
## file prefix is SM-7I1GZ
## zip_type is VIR
## uid is SM-7I1GZ
## file is SM-7LSBG.tar
## file prefix is SM-7LSBG
## zip_type is VIR
## uid is SM-7LSBG
## file is SM-7RYAM.tar
## file prefix is SM-7RYAM
## zip_type is VIR
## uid is SM-7RYAM
## file is SM-9JM8R.tar
## file prefix is SM-9JM8R
## zip_type is VIR
## uid is SM-9JM8R
## file is SM-9KOMU.tar
## file prefix is SM-9KOMU
## zip_type is VIR
## uid is SM-9KOMU
## file is SM-9RD3A.tar
## file prefix is SM-9RD3A
## zip_type is VIR
## uid is SM-9RD3A
## file is SM-9UKB5.tar
## file prefix is SM-9UKB5
## zip_type is VIR
## uid is SM-9UKB5
## file is SM-9WO1G.tar
## file prefix is SM-9WO1G
## zip_type is VIR
## uid is SM-9WO1G
## file is SM-9XKX1.tar
## file prefix is SM-9XKX1
## zip_type is VIR
## uid is SM-9XKX1
## file is SM-9ZK9I.tar
## file prefix is SM-9ZK9I
## zip_type is VIR
## uid is SM-9ZK9I
## file is SM-A1ZWF.tar
## file prefix is SM-A1ZWF
## zip_type is VIR
## uid is SM-A1ZWF
## file is SM-A4O3D.tar
## file prefix is SM-A4O3D
## zip_type is VIR
## uid is SM-A4O3D
## file is HSM67VDX_P.fastq.gz
## file prefix is HSM67VDX_P
## zip_type is P
## uid is HSM67VDX_P
## file is HSM67VEO.tar
## file prefix is HSM67VEO
## zip_type is STD
## uid is HSM67VEO
## file is HSM67VEQ.tar
## file prefix is HSM67VEQ
## zip_type is STD
## uid is HSM67VEQ
## file is HSM67VES.tar
## file prefix is HSM67VES
## zip_type is STD
## uid is HSM67VES
## file is HSM67VEU.tar
## file prefix is HSM67VEU
## zip_type is STD
## uid is HSM67VEU
## file is HSM67VEW.tar
## file prefix is HSM67VEW
## zip_type is STD
## uid is HSM67VEW
## file is HSM7CYZJ.tar
## file prefix is HSM7CYZJ
## zip_type is STD
## uid is HSM7CYZJ
## file is HSM7CYZL.tar
## file prefix is HSM7CYZL
## zip_type is STD
## uid is HSM7CYZL
## file is HSM7CYZR.tar
## file prefix is HSM7CYZR
## zip_type is STD
## uid is HSM7CYZR
## file is HSM7J4RE.tar
## file prefix is HSM7J4RE
## zip_type is STD
## uid is HSM7J4RE
## file is HSM7J4G1.tar
## file prefix is HSM7J4G1
## zip_type is STD
## uid is HSM7J4G1
## file is HSM7J4G8.tar
## file prefix is HSM7J4G8
## zip_type is STD
## uid is HSM7J4G8
## file is HSM7J4IU.tar
## file prefix is HSM7J4IU
## zip_type is STD
## uid is HSM7J4IU
## file is HSM7J4IW.tar
## file prefix is HSM7J4IW
## zip_type is STD
## uid is HSM7J4IW
## file is SM-791AN.tar
## file prefix is SM-791AN
## zip_type is VIR
## uid is SM-791AN
## file is SM-7CP35.tar
## file prefix is SM-7CP35
## zip_type is VIR
## uid is SM-7CP35
## file is SM-7ED9W.tar
## file prefix is SM-7ED9W
## zip_type is VIR
## uid is SM-7ED9W
## file is SM-7IWEQ.tar
## file prefix is SM-7IWEQ
## zip_type is VIR
## uid is SM-7IWEQ
## file is SM-9PJ1G.tar
## file prefix is SM-9PJ1G
## zip_type is VIR
## uid is SM-9PJ1G
## file is HSM67VDR_P.fastq.gz
## file prefix is HSM67VDR_P
## zip_type is P
## uid is HSM67VDR_P
## file is HSM6XRUL.tar
## file prefix is HSM6XRUL
## zip_type is STD
## uid is HSM6XRUL
## file is HSM6XRUN.tar
## file prefix is HSM6XRUN
## zip_type is STD
## uid is HSM6XRUN
## file is HSM6XRUR.tar
## file prefix is HSM6XRUR
## zip_type is STD
## uid is HSM6XRUR
## file is HSM6XRQ8.tar
## file prefix is HSM6XRQ8
## zip_type is STD
## uid is HSM6XRQ8
## file is HSM7CZ16.tar
## file prefix is HSM7CZ16
## zip_type is STD
## uid is HSM7CZ16
## file is HSM7CZ18.tar
## file prefix is HSM7CZ18
## zip_type is STD
## uid is HSM7CZ18
## file is HSM7CZ1A.tar
## file prefix is HSM7CZ1A
## zip_type is STD
## uid is HSM7CZ1A
## file is HSM7CZ1C.tar
## file prefix is HSM7CZ1C
## zip_type is STD
## uid is HSM7CZ1C
## file is HSM7CZ1E.tar
## file prefix is HSM7CZ1E
## zip_type is STD
## uid is HSM7CZ1E
## file is HSM7CZ1G.tar
## file prefix is HSM7CZ1G
## zip_type is STD
## uid is HSM7CZ1G
## file is HSM7J4HA.tar
## file prefix is HSM7J4HA
## zip_type is STD
## uid is HSM7J4HA
## file is HSM7J4HC.tar
## file prefix is HSM7J4HC
## zip_type is STD
## uid is HSM7J4HC
## file is HSM7J4HE.tar
## file prefix is HSM7J4HE
## zip_type is STD
## uid is HSM7J4HE
## file is HSM7J4HG.tar
## file prefix is HSM7J4HG
## zip_type is STD
## uid is HSM7J4HG
## file is HSM7J4HI.tar
## file prefix is HSM7J4HI
## zip_type is STD
## uid is HSM7J4HI
## file is HSM7J4HK.tar
## file prefix is HSM7J4HK
## zip_type is STD
## uid is HSM7J4HK
## file is HSM7J4KC.tar
## file prefix is HSM7J4KC
## zip_type is STD
## uid is HSM7J4KC
## file is HSM7J4KG.tar
## file prefix is HSM7J4KG
## zip_type is STD
## uid is HSM7J4KG
## file is HSM7J4KI.tar
## file prefix is HSM7J4KI
## zip_type is STD
## uid is HSM7J4KI
## file is HSM7J4KK.tar
## file prefix is HSM7J4KK
## zip_type is STD
## uid is HSM7J4KK
## file is HSM7J4KM.tar
## file prefix is HSM7J4KM
## zip_type is STD
## uid is HSM7J4KM
## file is SM-7BEZT.tar
## file prefix is SM-7BEZT
## zip_type is VIR
## uid is SM-7BEZT
## file is SM-7CP2K.tar
## file prefix is SM-7CP2K
## zip_type is VIR
## uid is SM-7CP2K
## file is SM-7GYHK.tar
## file prefix is SM-7GYHK
## zip_type is VIR
## uid is SM-7GYHK
## file is SM-7MCUI.tar
## file prefix is SM-7MCUI
## zip_type is VIR
## uid is SM-7MCUI
## file is SM-9NBDN.tar
## file prefix is SM-9NBDN
## zip_type is VIR
## uid is SM-9NBDN
## file is SM-A1C88.tar
## file prefix is SM-A1C88
## zip_type is VIR
## uid is SM-A1C88
## file is SM-AAHYG.tar
## file prefix is SM-AAHYG
## zip_type is VIR
## uid is SM-AAHYG
## file is HSM67VDT_P.fastq.gz
## file prefix is HSM67VDT_P
## zip_type is P
## uid is HSM67VDT_P
## file is HSM67VDT.tar
## file prefix is HSM67VDT
## zip_type is STD
## uid is HSM67VDT
## file is HSM67VI3.tar
## file prefix is HSM67VI3
## zip_type is STD
## uid is HSM67VI3
## file is HSM67VI5.tar
## file prefix is HSM67VI5
## zip_type is STD
## uid is HSM67VI5
## file is HSM67VI7.tar
## file prefix is HSM67VI7
## zip_type is STD
## uid is HSM67VI7
## file is HSM67VI9.tar
## file prefix is HSM67VI9
## zip_type is STD
## uid is HSM67VI9
## file is HSM67VIB.tar
## file prefix is HSM67VIB
## zip_type is STD
## uid is HSM67VIB
## file is HSM7CZ24.tar
## file prefix is HSM7CZ24
## zip_type is STD
## uid is HSM7CZ24
## file is HSM7CZ26.tar
## file prefix is HSM7CZ26
## zip_type is STD
## uid is HSM7CZ26
## file is HSM7CZ28.tar
## file prefix is HSM7CZ28
## zip_type is STD
## uid is HSM7CZ28
## file is HSM7CZ2A.tar
## file prefix is HSM7CZ2A
## zip_type is STD
## uid is HSM7CZ2A
## file is HSM7CZ2E.tar
## file prefix is HSM7CZ2E
## zip_type is STD
## uid is HSM7CZ2E
## file is HSM7J4HW.tar
## file prefix is HSM7J4HW
## zip_type is STD
## uid is HSM7J4HW
## file is HSM7J4HY.tar
## file prefix is HSM7J4HY
## zip_type is STD
## uid is HSM7J4HY
## file is HSM7J4I5.tar
## file prefix is HSM7J4I5
## zip_type is STD
## uid is HSM7J4I5
## file is HSM7J4I7.tar
## file prefix is HSM7J4I7
## zip_type is STD
## uid is HSM7J4I7
## file is HSM7J4I9.tar
## file prefix is HSM7J4I9
## zip_type is STD
## uid is HSM7J4I9
## file is HSM7J4O1.tar
## file prefix is HSM7J4O1
## zip_type is STD
## uid is HSM7J4O1
## file is HSM7J4NY.tar
## file prefix is HSM7J4NY
## zip_type is STD
## uid is HSM7J4NY
## file is HSM7J4O3.tar
## file prefix is HSM7J4O3
## zip_type is STD
## uid is HSM7J4O3
## file is HSM7J4O5.tar
## file prefix is HSM7J4O5
## zip_type is STD
## uid is HSM7J4O5
## file is HSM7J4O7.tar
## file prefix is HSM7J4O7
## zip_type is STD
## uid is HSM7J4O7
## file is HSM7J4O9.tar
## file prefix is HSM7J4O9
## zip_type is STD
## uid is HSM7J4O9
## file is SM-7BYBY.tar
## file prefix is SM-7BYBY
## zip_type is VIR
## uid is SM-7BYBY
## file is SM-7FK3H.tar
## file prefix is SM-7FK3H
## zip_type is VIR
## uid is SM-7FK3H
## file is SM-7I1GK.tar
## file prefix is SM-7I1GK
## zip_type is VIR
## uid is SM-7I1GK
## file is SM-9VWAZ.tar
## file prefix is SM-9VWAZ
## zip_type is VIR
## uid is SM-9VWAZ
## file is SM-A6ZJ6.tar
## file prefix is SM-A6ZJ6
## zip_type is VIR
## uid is SM-A6ZJ6
## file is SM-A7J1K.tar
## file prefix is SM-A7J1K
## zip_type is VIR
## uid is SM-A7J1K
## file is SM-AFCK4.tar
## file prefix is SM-AFCK4
## zip_type is VIR
## uid is SM-AFCK4
## file is HSM7CZ1T_P.fastq.gz
## file prefix is HSM7CZ1T_P
## zip_type is P
## uid is HSM7CZ1T_P
## file is HSM7CZ1V.tar
## file prefix is HSM7CZ1V
## zip_type is STD
## uid is HSM7CZ1V
## file is HSM7CZ1Z.tar
## file prefix is HSM7CZ1Z
## zip_type is STD
## uid is HSM7CZ1Z
## file is HSM7CZ22.tar
## file prefix is HSM7CZ22
## zip_type is STD
## uid is HSM7CZ22
## file is HSM7J4GP.tar
## file prefix is HSM7J4GP
## zip_type is STD
## uid is HSM7J4GP
## file is HSM7J4GR.tar
## file prefix is HSM7J4GR
## zip_type is STD
## uid is HSM7J4GR
## file is HSM7J4NC.tar
## file prefix is HSM7J4NC
## zip_type is STD
## uid is HSM7J4NC
## file is HSM7J4NE.tar
## file prefix is HSM7J4NE
## zip_type is STD
## uid is HSM7J4NE
## file is HSM7J4NM.tar
## file prefix is HSM7J4NM
## zip_type is STD
## uid is HSM7J4NM
## file is HSMA33O1.tar
## file prefix is HSMA33O1
## zip_type is STD
## uid is HSMA33O1
## file is HSMA33O3.tar
## file prefix is HSMA33O3
## zip_type is STD
## uid is HSMA33O3
## file is SM-7LDJL.tar
## file prefix is SM-7LDJL
## zip_type is VIR
## uid is SM-7LDJL
## file is SM-7R3A2.tar
## file prefix is SM-7R3A2
## zip_type is VIR
## uid is SM-7R3A2
## file is SM-9SIZO.tar
## file prefix is SM-9SIZO
## zip_type is VIR
## uid is SM-9SIZO
## file is SM-A2BG5.tar
## file prefix is SM-A2BG5
## zip_type is VIR
## uid is SM-A2BG5
## file is HSM67VID.tar
## file prefix is HSM67VID
## zip_type is STD
## uid is HSM67VID
## file is HSM67VIF.tar
## file prefix is HSM67VIF
## zip_type is STD
## uid is HSM67VIF
## file is HSM67VIJ_P.fastq.gz
## file prefix is HSM67VIJ_P
## zip_type is P
## uid is HSM67VIJ_P
## file is HSM67VIL.tar
## file prefix is HSM67VIL
## zip_type is STD
## uid is HSM67VIL
## file is HSM7J4GD.tar
## file prefix is HSM7J4GD
## zip_type is STD
## uid is HSM7J4GD
## file is HSM7J4MY.tar
## file prefix is HSM7J4MY
## zip_type is STD
## uid is HSM7J4MY
## file is HSM7J4N4.tar
## file prefix is HSM7J4N4
## zip_type is STD
## uid is HSM7J4N4
## file is HSM7J4N6.tar
## file prefix is HSM7J4N6
## zip_type is STD
## uid is HSM7J4N6
## file is HSM7J4N8.tar
## file prefix is HSM7J4N8
## zip_type is STD
## uid is HSM7J4N8
## file is HSM7J4NA.tar
## file prefix is HSM7J4NA
## zip_type is STD
## uid is HSM7J4NA
## file is SM-9IT17.tar
## file prefix is SM-9IT17
## zip_type is VIR
## uid is SM-9IT17
## file is SM-9QVJQ.tar
## file prefix is SM-9QVJQ
## zip_type is VIR
## uid is SM-9QVJQ
## file is SM-9ZK9E.tar
## file prefix is SM-9ZK9E
## zip_type is VIR
## uid is SM-9ZK9E
## file is SM-A5A1U.tar
## file prefix is SM-A5A1U
## zip_type is VIR
## uid is SM-A5A1U
## file is SM-A9J8X.tar
## file prefix is SM-A9J8X
## zip_type is VIR
## uid is SM-A9J8X
## file is HSM7CYWS_P.fastq.gz
## file prefix is HSM7CYWS_P
## zip_type is P
## uid is HSM7CYWS_P
## file is HSM7CZ2Z.tar
## file prefix is HSM7CZ2Z
## zip_type is STD
## uid is HSM7CZ2Z
## file is HSM7CZ32.tar
## file prefix is HSM7CZ32
## zip_type is STD
## uid is HSM7CZ32
## file is HSM7CZ36.tar
## file prefix is HSM7CZ36
## zip_type is STD
## uid is HSM7CZ36
## file is HSM7CZ38.tar
## file prefix is HSM7CZ38
## zip_type is STD
## uid is HSM7CZ38
## file is HSM7J4L5.tar
## file prefix is HSM7J4L5
## zip_type is STD
## uid is HSM7J4L5
## file is HSM7J4L9.tar
## file prefix is HSM7J4L9
## zip_type is STD
## uid is HSM7J4L9
## file is HSM7J4ON.tar
## file prefix is HSM7J4ON
## zip_type is STD
## uid is HSM7J4ON
## file is HSM7J4OP.tar
## file prefix is HSM7J4OP
## zip_type is STD
## uid is HSM7J4OP
## file is HSM7J4OT.tar
## file prefix is HSM7J4OT
## zip_type is STD
## uid is HSM7J4OT
## file is HSM7J4OV.tar
## file prefix is HSM7J4OV
## zip_type is STD
## uid is HSM7J4OV
## file is HSM7J4OX.tar
## file prefix is HSM7J4OX
## zip_type is STD
## uid is HSM7J4OX
## file is HSMA33JZ.tar
## file prefix is HSMA33JZ
## zip_type is STD
## uid is HSMA33JZ
## file is SM-7ORQ1.tar
## file prefix is SM-7ORQ1
## zip_type is VIR
## uid is SM-7ORQ1
## file is SM-9HYXY.tar
## file prefix is SM-9HYXY
## zip_type is VIR
## uid is SM-9HYXY
## file is SM-9MPH1.tar
## file prefix is SM-9MPH1
## zip_type is VIR
## uid is SM-9MPH1
## file is SM-9OS7O.tar
## file prefix is SM-9OS7O
## zip_type is VIR
## uid is SM-9OS7O
## file is SM-9WYEO.tar
## file prefix is SM-9WYEO
## zip_type is VIR
## uid is SM-9WYEO
## file is SM-A12KG.tar
## file prefix is SM-A12KG
## zip_type is VIR
## uid is SM-A12KG
## file is SM-A8VEW.tar
## file prefix is SM-A8VEW
## zip_type is VIR
## uid is SM-A8VEW
## file is SM-AFJ1O.tar
## file prefix is SM-AFJ1O
## zip_type is VIR
## uid is SM-AFJ1O
## file is HSM7CZ2V_P.fastq.gz
## file prefix is HSM7CZ2V_P
## zip_type is P
## uid is HSM7CZ2V_P
## file is HSM7CZ3A.tar
## file prefix is HSM7CZ3A
## zip_type is STD
## uid is HSM7CZ3A
## file is HSM7CZ3C.tar
## file prefix is HSM7CZ3C
## zip_type is STD
## uid is HSM7CZ3C
## file is HSM7CZ3E.tar
## file prefix is HSM7CZ3E
## zip_type is STD
## uid is HSM7CZ3E
## file is HSM7CZ3G.tar
## file prefix is HSM7CZ3G
## zip_type is STD
## uid is HSM7CZ3G
## file is HSM7J4PA.tar
## file prefix is HSM7J4PA
## zip_type is STD
## uid is HSM7J4PA
## file is HSM7J4MK.tar
## file prefix is HSM7J4MK
## zip_type is STD
## uid is HSM7J4MK
## file is HSM7J4OB.tar
## file prefix is HSM7J4OB
## zip_type is STD
## uid is HSM7J4OB
## file is HSM7J4OE.tar
## file prefix is HSM7J4OE
## zip_type is STD
## uid is HSM7J4OE
## file is HSM7J4OL.tar
## file prefix is HSM7J4OL
## zip_type is STD
## uid is HSM7J4OL
## file is HSMA33JB.tar
## file prefix is HSMA33JB
## zip_type is STD
## uid is HSMA33JB
## file is HSMA33JD.tar
## file prefix is HSMA33JD
## zip_type is STD
## uid is HSMA33JD
## file is SM-7R3A9.tar
## file prefix is SM-7R3A9
## zip_type is VIR
## uid is SM-7R3A9
## file is SM-9IMBK.tar
## file prefix is SM-9IMBK
## zip_type is VIR
## uid is SM-9IMBK
## file is SM-9LB5Z.tar
## file prefix is SM-9LB5Z
## zip_type is VIR
## uid is SM-9LB5Z
## file is SM-9XKWY.tar
## file prefix is SM-9XKWY
## zip_type is VIR
## uid is SM-9XKWY
## file is SM-A9F5A.tar
## file prefix is SM-A9F5A
## zip_type is VIR
## uid is SM-A9F5A
## file is SM-AKTYZ.tar
## file prefix is SM-AKTYZ
## zip_type is VIR
## uid is SM-AKTYZ
## file is SM-AMR2H.tar
## file prefix is SM-AMR2H
## zip_type is VIR
## uid is SM-AMR2H
## file is HSM7J4Q1.tar
## file prefix is HSM7J4Q1
## zip_type is STD
## uid is HSM7J4Q1
## file is HSM7J4Q3_P.fastq.gz
## file prefix is HSM7J4Q3_P
## zip_type is P
## uid is HSM7J4Q3_P
## file is HSM7J4Q7.tar
## file prefix is HSM7J4Q7
## zip_type is STD
## uid is HSM7J4Q7
## file is HSM7J4Q9.tar
## file prefix is HSM7J4Q9
## zip_type is STD
## uid is HSM7J4Q9
## file is HSM7J4MS.tar
## file prefix is HSM7J4MS
## zip_type is STD
## uid is HSM7J4MS
## file is HSM7J4MW.tar
## file prefix is HSM7J4MW
## zip_type is STD
## uid is HSM7J4MW
## file is HSM7J4IC.tar
## file prefix is HSM7J4IC
## zip_type is STD
## uid is HSM7J4IC
## file is HSM7J4OZ.tar
## file prefix is HSM7J4OZ
## zip_type is STD
## uid is HSM7J4OZ
## file is HSM7J4P2.tar
## file prefix is HSM7J4P2
## zip_type is STD
## uid is HSM7J4P2
## file is HSMA33MV.tar
## file prefix is HSMA33MV
## zip_type is STD
## uid is HSMA33MV
## file is HSMA33JN.tar
## file prefix is HSMA33JN
## zip_type is STD
## uid is HSMA33JN
## file is HSMA33JP.tar
## file prefix is HSMA33JP
## zip_type is STD
## uid is HSMA33JP
## file is HSMA33JR.tar
## file prefix is HSMA33JR
## zip_type is STD
## uid is HSMA33JR
## file is SM-9SNJA.tar
## file prefix is SM-9SNJA
## zip_type is VIR
## uid is SM-9SNJA
## file is SM-9ZE8G.tar
## file prefix is SM-9ZE8G
## zip_type is VIR
## uid is SM-9ZE8G
## file is SM-A77UJ.tar
## file prefix is SM-A77UJ
## zip_type is VIR
## uid is SM-A77UJ
## file is SM-APR9B.tar
## file prefix is SM-APR9B
## zip_type is VIR
## uid is SM-APR9B
## file is HSM7CZ2X_P.fastq.gz
## file prefix is HSM7CZ2X_P
## zip_type is P
## uid is HSM7CZ2X_P
## file is HSM7J4HM.tar
## file prefix is HSM7J4HM
## zip_type is STD
## uid is HSM7J4HM
## file is HSM7J4HO.tar
## file prefix is HSM7J4HO
## zip_type is STD
## uid is HSM7J4HO
## file is HSM7J4HQ.tar
## file prefix is HSM7J4HQ
## zip_type is STD
## uid is HSM7J4HQ
## file is HSM7J4HS.tar
## file prefix is HSM7J4HS
## zip_type is STD
## uid is HSM7J4HS
## file is HSM7J4HU.tar
## file prefix is HSM7J4HU
## zip_type is STD
## uid is HSM7J4HU
## file is HSM7J4JZ.tar
## file prefix is HSM7J4JZ
## zip_type is STD
## uid is HSM7J4JZ
## file is HSM7J4K2.tar
## file prefix is HSM7J4K2
## zip_type is STD
## uid is HSM7J4K2
## file is HSM7J4K4.tar
## file prefix is HSM7J4K4
## zip_type is STD
## uid is HSM7J4K4
## file is HSM7J4K6.tar
## file prefix is HSM7J4K6
## zip_type is STD
## uid is HSM7J4K6
## file is HSM7J4K8.tar
## file prefix is HSM7J4K8
## zip_type is STD
## uid is HSM7J4K8
## file is HSM7J4KA.tar
## file prefix is HSM7J4KA
## zip_type is STD
## uid is HSM7J4KA
## file is HSMA33OP.tar
## file prefix is HSMA33OP
## zip_type is STD
## uid is HSMA33OP
## file is HSMA33OR.tar
## file prefix is HSMA33OR
## zip_type is STD
## uid is HSMA33OR
## file is HSMA33OT.tar
## file prefix is HSMA33OT
## zip_type is STD
## uid is HSMA33OT
## file is HSMA33OV.tar
## file prefix is HSMA33OV
## zip_type is STD
## uid is HSMA33OV
## file is HSMA33OX.tar
## file prefix is HSMA33OX
## zip_type is STD
## uid is HSMA33OX
## file is HSMA33OZ.tar
## file prefix is HSMA33OZ
## zip_type is STD
## uid is HSMA33OZ
## file is HSMA33MA.tar
## file prefix is HSMA33MA
## zip_type is STD
## uid is HSMA33MA
## file is HSMA33MC.tar
## file prefix is HSMA33MC
## zip_type is STD
## uid is HSMA33MC
## file is HSMA33ME.tar
## file prefix is HSMA33ME
## zip_type is STD
## uid is HSMA33ME
## file is HSMA33MG.tar
## file prefix is HSMA33MG
## zip_type is STD
## uid is HSMA33MG
## file is HSMA33MI.tar
## file prefix is HSMA33MI
## zip_type is STD
## uid is HSMA33MI
## file is HSMA33MK.tar
## file prefix is HSMA33MK
## zip_type is STD
## uid is HSMA33MK
## file is SM-9NBEI.tar
## file prefix is SM-9NBEI
## zip_type is VIR
## uid is SM-9NBEI
## file is SM-9QMQ8.tar
## file prefix is SM-9QMQ8
## zip_type is VIR
## uid is SM-9QMQ8
## file is SM-9SIJH.tar
## file prefix is SM-9SIJH
## zip_type is VIR
## uid is SM-9SIJH
## file is SM-9VWCF.tar
## file prefix is SM-9VWCF
## zip_type is VIR
## uid is SM-9VWCF
## file is SM-9XIO2.tar
## file prefix is SM-9XIO2
## zip_type is VIR
## uid is SM-9XIO2
## file is SM-9ZA5K.tar
## file prefix is SM-9ZA5K
## zip_type is VIR
## uid is SM-9ZA5K
## file is SM-A3753.tar
## file prefix is SM-A3753
## zip_type is VIR
## uid is SM-A3753
## file is SM-A4GB9.tar
## file prefix is SM-A4GB9
## zip_type is VIR
## uid is SM-A4GB9
## file is SM-A77WD.tar
## file prefix is SM-A77WD
## zip_type is VIR
## uid is SM-A77WD
## file is SM-A9F6Z.tar
## file prefix is SM-A9F6Z
## zip_type is VIR
## uid is SM-A9F6Z
## file is SM-AGULX.tar
## file prefix is SM-AGULX
## zip_type is VIR
## uid is SM-AGULX
## file is SM-AIG5H.tar
## file prefix is SM-AIG5H
## zip_type is VIR
## uid is SM-AIG5H
## file is SM-AKTXH.tar
## file prefix is SM-AKTXH
## zip_type is VIR
## uid is SM-AKTXH
## file is SM-APDXE.tar
## file prefix is SM-APDXE
## zip_type is VIR
## uid is SM-APDXE
## file is SM-ARMBT.tar
## file prefix is SM-ARMBT
## zip_type is VIR
## uid is SM-ARMBT
## file is SM-AUIUT.tar
## file prefix is SM-AUIUT
## zip_type is VIR
## uid is SM-AUIUT
## file is SM-AVR74.tar
## file prefix is SM-AVR74
## zip_type is VIR
## uid is SM-AVR74
## file is SM-AZYBB.tar
## file prefix is SM-AZYBB
## zip_type is VIR
## uid is SM-AZYBB
## file is HSM7J4J7.tar
## file prefix is HSM7J4J7
## zip_type is STD
## uid is HSM7J4J7
## file is HSM7J4J9.tar
## file prefix is HSM7J4J9
## zip_type is STD
## uid is HSM7J4J9
## file is HSM7J4JD.tar
## file prefix is HSM7J4JD
## zip_type is STD
## uid is HSM7J4JD
## file is HSM7J4JF.tar
## file prefix is HSM7J4JF
## zip_type is STD
## uid is HSM7J4JF
## file is HSMA33NA.tar
## file prefix is HSMA33NA
## zip_type is STD
## uid is HSMA33NA
## file is HSMA33NC.tar
## file prefix is HSMA33NC
## zip_type is STD
## uid is HSMA33NC
## file is HSMA33NG_P.fastq.gz
## file prefix is HSMA33NG_P
## zip_type is P
## uid is HSMA33NG_P
## file is HSMA33KE.tar
## file prefix is HSMA33KE
## zip_type is STD
## uid is HSMA33KE
## file is HSMA33KM.tar
## file prefix is HSMA33KM
## zip_type is STD
## uid is HSMA33KM
## file is HSMA33KO.tar
## file prefix is HSMA33KO
## zip_type is STD
## uid is HSMA33KO
## file is HSMA33PX.tar
## file prefix is HSMA33PX
## zip_type is STD
## uid is HSMA33PX
## file is HSMA33PZ.tar
## file prefix is HSMA33PZ
## zip_type is STD
## uid is HSMA33PZ
## file is HSMA33Q6.tar
## file prefix is HSMA33Q6
## zip_type is STD
## uid is HSMA33Q6
## file is SM-9WOBE.tar
## file prefix is SM-9WOBE
## zip_type is VIR
## uid is SM-9WOBE
## file is SM-A3754.tar
## file prefix is SM-A3754
## zip_type is VIR
## uid is SM-A3754
## file is SM-AADL4.tar
## file prefix is SM-AADL4
## zip_type is VIR
## uid is SM-AADL4
## file is SM-B2OQM.tar
## file prefix is SM-B2OQM
## zip_type is VIR
## uid is SM-B2OQM
## file is SM-CE9V6.tar
## file prefix is SM-CE9V6
## zip_type is VIR
## uid is SM-CE9V6
## file is HSM7J4LP.tar
## file prefix is HSM7J4LP
## zip_type is STD
## uid is HSM7J4LP
## file is HSM7J4JH.tar
## file prefix is HSM7J4JH
## zip_type is STD
## uid is HSM7J4JH
## file is HSM7J4JJ.tar
## file prefix is HSM7J4JJ
## zip_type is STD
## uid is HSM7J4JJ
## file is HSM7J4JN.tar
## file prefix is HSM7J4JN
## zip_type is STD
## uid is HSM7J4JN
## file is HSM7J4JP.tar
## file prefix is HSM7J4JP
## zip_type is STD
## uid is HSM7J4JP
## file is HSMA33NO.tar
## file prefix is HSMA33NO
## zip_type is STD
## uid is HSMA33NO
## file is HSMA33NQ.tar
## file prefix is HSMA33NQ
## zip_type is STD
## uid is HSMA33NQ
## file is HSMA33KQ.tar
## file prefix is HSMA33KQ
## zip_type is STD
## uid is HSMA33KQ
## file is HSMA33KS.tar
## file prefix is HSMA33KS
## zip_type is STD
## uid is HSMA33KS
## file is HSMA33KU.tar
## file prefix is HSMA33KU
## zip_type is STD
## uid is HSMA33KU
## file is HSMA33L1.tar
## file prefix is HSMA33L1
## zip_type is STD
## uid is HSMA33L1
## file is HSMA33PL.tar
## file prefix is HSMA33PL
## zip_type is STD
## uid is HSMA33PL
## file is HSMA33PN.tar
## file prefix is HSMA33PN
## zip_type is STD
## uid is HSMA33PN
## file is SM-9VWC7.tar
## file prefix is SM-9VWC7
## zip_type is VIR
## uid is SM-9VWC7
## file is SM-9XIDV.tar
## file prefix is SM-9XIDV
## zip_type is VIR
## uid is SM-9XIDV
## file is SM-9ZA1P.tar
## file prefix is SM-9ZA1P
## zip_type is VIR
## uid is SM-9ZA1P
## file is SM-A4O3V.tar
## file prefix is SM-A4O3V
## zip_type is VIR
## uid is SM-A4O3V
## file is SM-AF6MT.tar
## file prefix is SM-AF6MT
## zip_type is VIR
## uid is SM-AF6MT
## file is SM-AGUMZ.tar
## file prefix is SM-AGUMZ
## zip_type is VIR
## uid is SM-AGUMZ
## file is SM-ARMCW.tar
## file prefix is SM-ARMCW
## zip_type is VIR
## uid is SM-ARMCW
## file is SM-AVR7A.tar
## file prefix is SM-AVR7A
## zip_type is VIR
## uid is SM-AVR7A
## file is SM-AXQR3.tar
## file prefix is SM-AXQR3
## zip_type is VIR
## uid is SM-AXQR3
## file is HSM7J4JT_P.fastq.gz
## file prefix is HSM7J4JT_P
## zip_type is P
## uid is HSM7J4JT_P
## file is HSM7J4NO.tar
## file prefix is HSM7J4NO
## zip_type is STD
## uid is HSM7J4NO
## file is HSM7J4NS.tar
## file prefix is HSM7J4NS
## zip_type is STD
## uid is HSM7J4NS
## file is HSM7J4NU.tar
## file prefix is HSM7J4NU
## zip_type is STD
## uid is HSM7J4NU
## file is HSMA33OD.tar
## file prefix is HSMA33OD
## zip_type is STD
## uid is HSMA33OD
## file is HSMA33OJ.tar
## file prefix is HSMA33OJ
## zip_type is STD
## uid is HSMA33OJ
## file is HSMA33OL.tar
## file prefix is HSMA33OL
## zip_type is STD
## uid is HSMA33OL
## file is HSMA33LX.tar
## file prefix is HSMA33LX
## zip_type is STD
## uid is HSMA33LX
## file is HSMA33LZ.tar
## file prefix is HSMA33LZ
## zip_type is STD
## uid is HSMA33LZ
## file is HSMA33M2.tar
## file prefix is HSMA33M2
## zip_type is STD
## uid is HSMA33M2
## file is HSMA33M8.tar
## file prefix is HSMA33M8
## zip_type is STD
## uid is HSMA33M8
## file is SM-A373M.tar
## file prefix is SM-A373M
## zip_type is VIR
## uid is SM-A373M
## file is SM-A8VEX.tar
## file prefix is SM-A8VEX
## zip_type is VIR
## uid is SM-A8VEX
## file is SM-ADDEV.tar
## file prefix is SM-ADDEV
## zip_type is VIR
## uid is SM-ADDEV
## file is SM-AJ651.tar
## file prefix is SM-AJ651
## zip_type is VIR
## uid is SM-AJ651
## file is SM-ARB6B.tar
## file prefix is SM-ARB6B
## zip_type is VIR
## uid is SM-ARB6B
## file is SM-B5ORM.tar
## file prefix is SM-B5ORM
## zip_type is VIR
## uid is SM-B5ORM
## file is SM-CG7AX.tar
## file prefix is SM-CG7AX
## zip_type is VIR
## uid is SM-CG7AX
## file is HSM7J4JV_P.fastq.gz
## file prefix is HSM7J4JV_P
## zip_type is P
## uid is HSM7J4JV_P
## file is HSMA33MX.tar
## file prefix is HSMA33MX
## zip_type is STD
## uid is HSMA33MX
## file is HSMA33MZ.tar
## file prefix is HSMA33MZ
## zip_type is STD
## uid is HSMA33MZ
## file is HSMA33N4.tar
## file prefix is HSMA33N4
## zip_type is STD
## uid is HSMA33N4
## file is HSMA33IO.tar
## file prefix is HSMA33IO
## zip_type is STD
## uid is HSMA33IO
## file is HSMA33IS.tar
## file prefix is HSMA33IS
## zip_type is STD
## uid is HSMA33IS
## file is HSMA33RD.tar
## file prefix is HSMA33RD
## zip_type is STD
## uid is HSMA33RD
## file is HSMA33RF.tar
## file prefix is HSMA33RF
## zip_type is STD
## uid is HSMA33RF
## file is SM-ADID6.tar
## file prefix is SM-ADID6
## zip_type is VIR
## uid is SM-ADID6
## file is SM-C1MZ9.tar
## file prefix is SM-C1MZ9
## zip_type is VIR
## uid is SM-C1MZ9
## file is SM-CG9XT.tar
## file prefix is SM-CG9XT
## zip_type is VIR
## uid is SM-CG9XT
## file is HSMA33IE.tar
## file prefix is HSMA33IE
## zip_type is STD
## uid is HSMA33IE
## file is HSMA33IG.tar
## file prefix is HSMA33IG
## zip_type is STD
## uid is HSMA33IG
## file is HSMA33IK.tar
## file prefix is HSMA33IK
## zip_type is STD
## uid is HSMA33IK
## file is HSMA33LP.tar
## file prefix is HSMA33LP
## zip_type is STD
## uid is HSMA33LP
## file is HSMA33QY.tar
## file prefix is HSMA33QY
## zip_type is STD
## uid is HSMA33QY
## file is HSMA33R1.tar
## file prefix is HSMA33R1
## zip_type is STD
## uid is HSMA33R1
## file is HSMA33R5.tar
## file prefix is HSMA33R5
## zip_type is STD
## uid is HSMA33R5
## file is HSMA33R7.tar
## file prefix is HSMA33R7
## zip_type is STD
## uid is HSMA33R7
## file is HSMA33R9.tar
## file prefix is HSMA33R9
## zip_type is STD
## uid is HSMA33R9
## file is HSMA33S4.tar
## file prefix is HSMA33S4
## zip_type is STD
## uid is HSMA33S4
## file is SM-AF6NM.tar
## file prefix is SM-AF6NM
## zip_type is VIR
## uid is SM-AF6NM
## file is SM-AL5F8.tar
## file prefix is SM-AL5F8
## zip_type is VIR
## uid is SM-AL5F8
## file is SM-B1JON.tar
## file prefix is SM-B1JON
## zip_type is VIR
## uid is SM-B1JON
## file is SM-BYX23.tar
## file prefix is SM-BYX23
## zip_type is VIR
## uid is SM-BYX23
## file is HSMA33NW.tar
## file prefix is HSMA33NW
## zip_type is STD
## uid is HSMA33NW
## file is HSMA33P2_P.fastq.gz
## file prefix is HSMA33P2_P
## zip_type is P
## uid is HSMA33P2_P
## file is HSMA33P6.tar
## file prefix is HSMA33P6
## zip_type is STD
## uid is HSMA33P6
## file is HSMA33IA.tar
## file prefix is HSMA33IA
## zip_type is STD
## uid is HSMA33IA
## file is HSMA33IC.tar
## file prefix is HSMA33IC
## zip_type is STD
## uid is HSMA33IC
## file is HSMA33LB.tar
## file prefix is HSMA33LB
## zip_type is STD
## uid is HSMA33LB
## file is HSMA33LH.tar
## file prefix is HSMA33LH
## zip_type is STD
## uid is HSMA33LH
## file is HSMA33LH_TR.tar
## file prefix is HSMA33LH_TR
## zip_type is TR
## uid is HSMA33LH_TR
## file is HSMA33LJ.tar
## file prefix is HSMA33LJ
## zip_type is STD
## uid is HSMA33LJ
## file is HSMA33RR.tar
## file prefix is HSMA33RR
## zip_type is STD
## uid is HSMA33RR
## file is HSMA33RT.tar
## file prefix is HSMA33RT
## zip_type is STD
## uid is HSMA33RT
## file is HSMA33RX.tar
## file prefix is HSMA33RX
## zip_type is STD
## uid is HSMA33RX
## file is HSMA33RX_TR.tar
## file prefix is HSMA33RX_TR
## zip_type is TR
## uid is HSMA33RX_TR
## file is SM-AJ64W.tar
## file prefix is SM-AJ64W
## zip_type is VIR
## uid is SM-AJ64W
## file is SM-AJDLQ.tar
## file prefix is SM-AJDLQ
## zip_type is VIR
## uid is SM-AJDLQ
## file is SM-ALGUL.tar
## file prefix is SM-ALGUL
## zip_type is VIR
## uid is SM-ALGUL
## file is SM-ASVMD.tar
## file prefix is SM-ASVMD
## zip_type is VIR
## uid is SM-ASVMD
## file is SM-AZAHH.tar
## file prefix is SM-AZAHH
## zip_type is VIR
## uid is SM-AZAHH
## file is SM-B1F5B.tar
## file prefix is SM-B1F5B
## zip_type is VIR
## uid is SM-B1F5B
## file is SM-CJFJI.tar
## file prefix is SM-CJFJI
## zip_type is VIR
## uid is SM-CJFJI
## file is SM-CSF6Z.tar
## file prefix is SM-CSF6Z
## zip_type is VIR
## uid is SM-CSF6Z
## file is HSMA33NY.tar
## file prefix is HSMA33NY
## zip_type is STD
## uid is HSMA33NY
## file is HSMA33J1_P.fastq.gz
## file prefix is HSMA33J1_P
## zip_type is P
## uid is HSMA33J1_P
## file is HSMA33J3.tar
## file prefix is HSMA33J3
## zip_type is STD
## uid is HSMA33J3
## file is HSMA33J5.tar
## file prefix is HSMA33J5
## zip_type is STD
## uid is HSMA33J5
## file is HSMA33J7.tar
## file prefix is HSMA33J7
## zip_type is STD
## uid is HSMA33J7
## file is HSMA33J9.tar
## file prefix is HSMA33J9
## zip_type is STD
## uid is HSMA33J9
## file is HSMA33MS.tar
## file prefix is HSMA33MS
## zip_type is STD
## uid is HSMA33MS
## file is HSMA33QM.tar
## file prefix is HSMA33QM
## zip_type is STD
## uid is HSMA33QM
## file is HSMA33QO.tar
## file prefix is HSMA33QO
## zip_type is STD
## uid is HSMA33QO
## file is HSMA33SE.tar
## file prefix is HSMA33SE
## zip_type is STD
## uid is HSMA33SE
## file is HSMA33SG.tar
## file prefix is HSMA33SG
## zip_type is STD
## uid is HSMA33SG
## file is HSMA33SI.tar
## file prefix is HSMA33SI
## zip_type is STD
## uid is HSMA33SI
## file is HSMA33SK.tar
## file prefix is HSMA33SK
## zip_type is STD
## uid is HSMA33SK
## file is SM-AL5ES.tar
## file prefix is SM-AL5ES
## zip_type is VIR
## uid is SM-AL5ES
## file is SM-AN5FS.tar
## file prefix is SM-AN5FS
## zip_type is VIR
## uid is SM-AN5FS
## file is SM-AWD4T.tar
## file prefix is SM-AWD4T
## zip_type is VIR
## uid is SM-AWD4T
## file is SM-BYMCE.tar
## file prefix is SM-BYMCE
## zip_type is VIR
## uid is SM-BYMCE
## file is SM-CLM5K.tar
## file prefix is SM-CLM5K
## zip_type is VIR
## uid is SM-CLM5K
## file is SM-CNQYZ.tar
## file prefix is SM-CNQYZ
## zip_type is VIR
## uid is SM-CNQYZ
## file is SM-CQEDH.tar
## file prefix is SM-CQEDH
## zip_type is VIR
## uid is SM-CQEDH
## file is MSM5LLDI.tar
## file prefix is MSM5LLDI
## zip_type is STD
## uid is MSM5LLDI
## file is MSM5LLDK.tar
## file prefix is MSM5LLDK
## zip_type is STD
## uid is MSM5LLDK
## file is MSM5LLDM.tar
## file prefix is MSM5LLDM
## zip_type is STD
## uid is MSM5LLDM
## file is MSM5LLDQ.tar
## file prefix is MSM5LLDQ
## zip_type is STD
## uid is MSM5LLDQ
## file is MSM5LLDS.tar
## file prefix is MSM5LLDS
## zip_type is STD
## uid is MSM5LLDS
## file is MSM5LLDU.tar
## file prefix is MSM5LLDU
## zip_type is STD
## uid is MSM5LLDU
## file is MSM5LLHX_P.fastq.gz
## file prefix is MSM5LLHX_P
## zip_type is P
## uid is MSM5LLHX_P
## file is MSM5LLI2_P.fastq.gz
## file prefix is MSM5LLI2_P
## zip_type is P
## uid is MSM5LLI2_P
## file is MSM5LLI8_P.fastq.gz
## file prefix is MSM5LLI8_P
## zip_type is P
## uid is MSM5LLI8_P
## file is MSM5LLI6_P.fastq.gz
## file prefix is MSM5LLI6_P
## zip_type is P
## uid is MSM5LLI6_P
## file is MSM5LLI4_P.fastq.gz
## file prefix is MSM5LLI4_P
## zip_type is P
## uid is MSM5LLI4_P
## file is MSM5LLHQ_P.fastq.gz
## file prefix is MSM5LLHQ_P
## zip_type is P
## uid is MSM5LLHQ_P
## file is MSM5LLD6_P.fastq.gz
## file prefix is MSM5LLD6_P
## zip_type is P
## uid is MSM5LLD6_P
## file is MSM5LLDC.tar
## file prefix is MSM5LLDC
## zip_type is STD
## uid is MSM5LLDC
## file is MSM5LLDE.tar
## file prefix is MSM5LLDE
## zip_type is STD
## uid is MSM5LLDE
## file is MSM5LLDA.tar
## file prefix is MSM5LLDA
## zip_type is STD
## uid is MSM5LLDA
## file is MSM6J2LB.tar
## file prefix is MSM6J2LB
## zip_type is STD
## uid is MSM6J2LB
## file is SM-B4A4F.tar
## file prefix is SM-B4A4F
## zip_type is VIR
## uid is SM-B4A4F
## file is SM-B4A4J.tar
## file prefix is SM-B4A4J
## zip_type is VIR
## uid is SM-B4A4J
## file is SM-B44U2.tar
## file prefix is SM-B44U2
## zip_type is VIR
## uid is SM-B44U2
## file is SM-B44TV.tar
## file prefix is SM-B44TV
## zip_type is VIR
## uid is SM-B44TV
## file is SM-69UMZ.tar
## file prefix is SM-69UMZ
## zip_type is VIR
## uid is SM-69UMZ
## file is SM-B44U1.tar
## file prefix is SM-B44U1
## zip_type is VIR
## uid is SM-B44U1
## file is SM-6VHL9.tar
## file prefix is SM-6VHL9
## zip_type is VIR
## uid is SM-6VHL9
## file is SM-6XTSW.tar
## file prefix is SM-6XTSW
## zip_type is VIR
## uid is SM-6XTSW
## file is SM-7D34L.tar
## file prefix is SM-7D34L
## zip_type is VIR
## uid is SM-7D34L
## file is MSM5FZ9N_P.fastq.gz
## file prefix is MSM5FZ9N_P
## zip_type is P
## uid is MSM5FZ9N_P
## file is MSM5LLE3_P.fastq.gz
## file prefix is MSM5LLE3_P
## zip_type is P
## uid is MSM5LLE3_P
## file is MSM5LLE9_P.fastq.gz
## file prefix is MSM5LLE9_P
## zip_type is P
## uid is MSM5LLE9_P
## file is MSM5FZ9X_P.fastq.gz
## file prefix is MSM5FZ9X_P
## zip_type is P
## uid is MSM5FZ9X_P
## file is MSM5FZA2_P.fastq.gz
## file prefix is MSM5FZA2_P
## zip_type is P
## uid is MSM5FZA2_P
## file is MSM5ZOJY_P.fastq.gz
## file prefix is MSM5ZOJY_P
## zip_type is P
## uid is MSM5ZOJY_P
## file is MSM633FF_P.fastq.gz
## file prefix is MSM633FF_P
## zip_type is P
## uid is MSM633FF_P
## file is CSM5LLGB_P.fastq.gz
## file prefix is CSM5LLGB_P
## zip_type is P
## uid is CSM5LLGB_P
## file is MSM5LLGH_P.fastq.gz
## file prefix is MSM5LLGH_P
## zip_type is P
## uid is MSM5LLGH_P
## file is MSM5LLGF_P.fastq.gz
## file prefix is MSM5LLGF_P
## zip_type is P
## uid is MSM5LLGF_P
## file is MSM5LLGD_P.fastq.gz
## file prefix is MSM5LLGD_P
## zip_type is P
## uid is MSM5LLGD_P
## file is MSM5LLGJ_P.fastq.gz
## file prefix is MSM5LLGJ_P
## zip_type is P
## uid is MSM5LLGJ_P
## file is MSM5LLGL.tar
## file prefix is MSM5LLGL
## zip_type is STD
## uid is MSM5LLGL
## file is MSM6J2IQ.tar
## file prefix is MSM6J2IQ
## zip_type is STD
## uid is MSM6J2IQ
## file is MSM6J2IY.tar
## file prefix is MSM6J2IY
## zip_type is STD
## uid is MSM6J2IY
## file is MSM6J2J1.tar
## file prefix is MSM6J2J1
## zip_type is STD
## uid is MSM6J2J1
## file is MSM6J2QR.tar
## file prefix is MSM6J2QR
## zip_type is STD
## uid is MSM6J2QR
## file is SM-6WOC3.tar
## file prefix is SM-6WOC3
## zip_type is VIR
## uid is SM-6WOC3
## file is SM-77FXO.tar
## file prefix is SM-77FXO
## zip_type is VIR
## uid is SM-77FXO
## file is MSM5LLHR_P.fastq.gz
## file prefix is MSM5LLHR_P
## zip_type is P
## uid is MSM5LLHR_P
## file is MSM5LLIE_P.fastq.gz
## file prefix is MSM5LLIE_P
## zip_type is P
## uid is MSM5LLIE_P
## file is MSM5LLIG_P.fastq.gz
## file prefix is MSM5LLIG_P
## zip_type is P
## uid is MSM5LLIG_P
## file is MSM5LLIK_P.fastq.gz
## file prefix is MSM5LLIK_P
## zip_type is P
## uid is MSM5LLIK_P
## file is MSM5LLIM_P.fastq.gz
## file prefix is MSM5LLIM_P
## zip_type is P
## uid is MSM5LLIM_P
## file is MSM5LLEP.tar
## file prefix is MSM5LLEP
## zip_type is STD
## uid is MSM5LLEP
## file is MSM5LLER.tar
## file prefix is MSM5LLER
## zip_type is STD
## uid is MSM5LLER
## file is MSM6J2LT.tar
## file prefix is MSM6J2LT
## zip_type is STD
## uid is MSM6J2LT
## file is SM-6V47A.tar
## file prefix is SM-6V47A
## zip_type is VIR
## uid is SM-6V47A
## file is MSM5LLHV_P.fastq.gz
## file prefix is MSM5LLHV_P
## zip_type is P
## uid is MSM5LLHV_P
## file is MSM5LLIC_P.fastq.gz
## file prefix is MSM5LLIC_P
## zip_type is P
## uid is MSM5LLIC_P
## file is MSM5LLHE_P.fastq.gz
## file prefix is MSM5LLHE_P
## zip_type is P
## uid is MSM5LLHE_P
## file is MSM5LLHG_P.fastq.gz
## file prefix is MSM5LLHG_P
## zip_type is P
## uid is MSM5LLHG_P
## file is MSM5LLHI_P.fastq.gz
## file prefix is MSM5LLHI_P
## zip_type is P
## uid is MSM5LLHI_P
## file is MSM5LLHO_P.fastq.gz
## file prefix is MSM5LLHO_P
## zip_type is P
## uid is MSM5LLHO_P
## file is MSM6J2J3.tar
## file prefix is MSM6J2J3
## zip_type is STD
## uid is MSM6J2J3
## file is MSM6J2J5.tar
## file prefix is MSM6J2J5
## zip_type is STD
## uid is MSM6J2J5
## file is MSM6J2JB.tar
## file prefix is MSM6J2JB
## zip_type is STD
## uid is MSM6J2JB
## file is MSM6J2JD.tar
## file prefix is MSM6J2JD
## zip_type is STD
## uid is MSM6J2JD
## file is MSM6J2PK.tar
## file prefix is MSM6J2PK
## zip_type is STD
## uid is MSM6J2PK
## file is MSM6J2PM.tar
## file prefix is MSM6J2PM
## zip_type is STD
## uid is MSM6J2PM
## file is SM-6WOCA.tar
## file prefix is SM-6WOCA
## zip_type is VIR
## uid is SM-6WOCA
## file is SM-6Y2UM.tar
## file prefix is SM-6Y2UM
## zip_type is VIR
## uid is SM-6Y2UM
## file is SM-734RH.tar
## file prefix is SM-734RH
## zip_type is VIR
## uid is SM-734RH
## file is SM-7ED9G.tar
## file prefix is SM-7ED9G
## zip_type is VIR
## uid is SM-7ED9G
## file is MSM5LLIQ_P.fastq.gz
## file prefix is MSM5LLIQ_P
## zip_type is P
## uid is MSM5LLIQ_P
## file is MSM5LLFK_P.fastq.gz
## file prefix is MSM5LLFK_P
## zip_type is P
## uid is MSM5LLFK_P
## file is MSM5LLFM_P.fastq.gz
## file prefix is MSM5LLFM_P
## zip_type is P
## uid is MSM5LLFM_P
## file is MSM5LLFO_P.fastq.gz
## file prefix is MSM5LLFO_P
## zip_type is P
## uid is MSM5LLFO_P
## file is MSM5LLFU_P.fastq.gz
## file prefix is MSM5LLFU_P
## zip_type is P
## uid is MSM5LLFU_P
## file is MSM6J2HB.tar
## file prefix is MSM6J2HB
## zip_type is STD
## uid is MSM6J2HB
## file is MSM6J2HD.tar
## file prefix is MSM6J2HD
## zip_type is STD
## uid is MSM6J2HD
## file is MSM6J2HF.tar
## file prefix is MSM6J2HF
## zip_type is STD
## uid is MSM6J2HF
## file is MSM6J2HH.tar
## file prefix is MSM6J2HH
## zip_type is STD
## uid is MSM6J2HH
## file is MSM6J2HJ.tar
## file prefix is MSM6J2HJ
## zip_type is STD
## uid is MSM6J2HJ
## file is MSM6J2HL.tar
## file prefix is MSM6J2HL
## zip_type is STD
## uid is MSM6J2HL
## file is MSM6J2OH.tar
## file prefix is MSM6J2OH
## zip_type is STD
## uid is MSM6J2OH
## file is MSM6J2OJ.tar
## file prefix is MSM6J2OJ
## zip_type is STD
## uid is MSM6J2OJ
## file is MSM6J2OL.tar
## file prefix is MSM6J2OL
## zip_type is STD
## uid is MSM6J2OL
## file is MSM6J2ON.tar
## file prefix is MSM6J2ON
## zip_type is STD
## uid is MSM6J2ON
## file is MSM6J2OP.tar
## file prefix is MSM6J2OP
## zip_type is STD
## uid is MSM6J2OP
## file is MSM79HEA.tar
## file prefix is MSM79HEA
## zip_type is STD
## uid is MSM79HEA
## file is SM-7H4H2.tar
## file prefix is SM-7H4H2
## zip_type is VIR
## uid is SM-7H4H2
## file is MSM5LLIS_P.fastq.gz
## file prefix is MSM5LLIS_P
## zip_type is P
## uid is MSM5LLIS_P
## file is MSM5LLH2_P.fastq.gz
## file prefix is MSM5LLH2_P
## zip_type is P
## uid is MSM5LLH2_P
## file is MSM5LLH4_P.fastq.gz
## file prefix is MSM5LLH4_P
## zip_type is P
## uid is MSM5LLH4_P
## file is MSM5LLH8_P.fastq.gz
## file prefix is MSM5LLH8_P
## zip_type is P
## uid is MSM5LLH8_P
## file is MSM5LLHA_P.fastq.gz
## file prefix is MSM5LLHA_P
## zip_type is P
## uid is MSM5LLHA_P
## file is MSM5LLHC.tar
## file prefix is MSM5LLHC
## zip_type is STD
## uid is MSM5LLHC
## file is MSM6J2KC.tar
## file prefix is MSM6J2KC
## zip_type is STD
## uid is MSM6J2KC
## file is MSM6J2KE.tar
## file prefix is MSM6J2KE
## zip_type is STD
## uid is MSM6J2KE
## file is MSM6J2KM.tar
## file prefix is MSM6J2KM
## zip_type is STD
## uid is MSM6J2KM
## file is MSM6J2R2.tar
## file prefix is MSM6J2R2
## zip_type is STD
## uid is MSM6J2R2
## file is MSM6J2R8.tar
## file prefix is MSM6J2R8
## zip_type is STD
## uid is MSM6J2R8
## file is MSM6J2RA.tar
## file prefix is MSM6J2RA
## zip_type is STD
## uid is MSM6J2RA
## file is MSM6J2RC.tar
## file prefix is MSM6J2RC
## zip_type is STD
## uid is MSM6J2RC
## file is SM-6X9VY.tar
## file prefix is SM-6X9VY
## zip_type is VIR
## uid is SM-6X9VY
## file is SM-6XSKU.tar
## file prefix is SM-6XSKU
## zip_type is VIR
## uid is SM-6XSKU
## file is SM-6ZEV2.tar
## file prefix is SM-6ZEV2
## zip_type is VIR
## uid is SM-6ZEV2
## file is SM-74Y7I.tar
## file prefix is SM-74Y7I
## zip_type is VIR
## uid is SM-74Y7I
## file is SM-77FYE.tar
## file prefix is SM-77FYE
## zip_type is VIR
## uid is SM-77FYE
## file is SM-7DN34.tar
## file prefix is SM-7DN34
## zip_type is VIR
## uid is SM-7DN34
## file is SM-7KMR5.tar
## file prefix is SM-7KMR5
## zip_type is VIR
## uid is SM-7KMR5
## file is MSM5LLFG_P.fastq.gz
## file prefix is MSM5LLFG_P
## zip_type is P
## uid is MSM5LLFG_P
## file is MSM5LLGN_P.fastq.gz
## file prefix is MSM5LLGN_P
## zip_type is P
## uid is MSM5LLGN_P
## file is MSM5LLGR_P.fastq.gz
## file prefix is MSM5LLGR_P
## zip_type is P
## uid is MSM5LLGR_P
## file is MSM6J2IE.tar
## file prefix is MSM6J2IE
## zip_type is STD
## uid is MSM6J2IE
## file is MSM6J2IG.tar
## file prefix is MSM6J2IG
## zip_type is STD
## uid is MSM6J2IG
## file is MSM6J2II.tar
## file prefix is MSM6J2II
## zip_type is STD
## uid is MSM6J2II
## file is MSM6J2IK.tar
## file prefix is MSM6J2IK
## zip_type is STD
## uid is MSM6J2IK
## file is MSM6J2IM.tar
## file prefix is MSM6J2IM
## zip_type is STD
## uid is MSM6J2IM
## file is MSM6J2IO.tar
## file prefix is MSM6J2IO
## zip_type is STD
## uid is MSM6J2IO
## file is MSM6J2Q3.tar
## file prefix is MSM6J2Q3
## zip_type is STD
## uid is MSM6J2Q3
## file is MSM6J2Q5.tar
## file prefix is MSM6J2Q5
## zip_type is STD
## uid is MSM6J2Q5
## file is MSM6J2Q7.tar
## file prefix is MSM6J2Q7
## zip_type is STD
## uid is MSM6J2Q7
## file is MSM6J2Q9.tar
## file prefix is MSM6J2Q9
## zip_type is STD
## uid is MSM6J2Q9
## file is MSM6J2QB.tar
## file prefix is MSM6J2QB
## zip_type is STD
## uid is MSM6J2QB
## file is MSM6J2QD.tar
## file prefix is MSM6J2QD
## zip_type is STD
## uid is MSM6J2QD
## file is MSM79H63.tar
## file prefix is MSM79H63
## zip_type is STD
## uid is MSM79H63
## file is MSM79H65.tar
## file prefix is MSM79H65
## zip_type is STD
## uid is MSM79H65
## file is MSM79H67.tar
## file prefix is MSM79H67
## zip_type is STD
## uid is MSM79H67
## file is MSM79H69.tar
## file prefix is MSM79H69
## zip_type is STD
## uid is MSM79H69
## file is MSM79H6B.tar
## file prefix is MSM79H6B
## zip_type is STD
## uid is MSM79H6B
## file is SM-6V2D9.tar
## file prefix is SM-6V2D9
## zip_type is VIR
## uid is SM-6V2D9
## file is SM-76EOE.tar
## file prefix is SM-76EOE
## zip_type is VIR
## uid is SM-76EOE
## file is SM-7IL11.tar
## file prefix is SM-7IL11
## zip_type is VIR
## uid is SM-7IL11
## file is MSM5LLF2_P.fastq.gz
## file prefix is MSM5LLF2_P
## zip_type is P
## uid is MSM5LLF2_P
## file is MSM5LLF4.tar
## file prefix is MSM5LLF4
## zip_type is STD
## uid is MSM5LLF4
## file is MSM5LLF6.tar
## file prefix is MSM5LLF6
## zip_type is STD
## uid is MSM5LLF6
## file is MSM5LLF8.tar
## file prefix is MSM5LLF8
## zip_type is STD
## uid is MSM5LLF8
## file is MSM6J2LH.tar
## file prefix is MSM6J2LH
## zip_type is STD
## uid is MSM6J2LH
## file is MSM6J2LJ.tar
## file prefix is MSM6J2LJ
## zip_type is STD
## uid is MSM6J2LJ
## file is MSM6J2LL.tar
## file prefix is MSM6J2LL
## zip_type is STD
## uid is MSM6J2LL
## file is MSM6J2LN.tar
## file prefix is MSM6J2LN
## zip_type is STD
## uid is MSM6J2LN
## file is MSM6J2LR.tar
## file prefix is MSM6J2LR
## zip_type is STD
## uid is MSM6J2LR
## file is MSM6J2MB.tar
## file prefix is MSM6J2MB
## zip_type is STD
## uid is MSM6J2MB
## file is MSM6J2MD.tar
## file prefix is MSM6J2MD
## zip_type is STD
## uid is MSM6J2MD
## file is MSM6J2MF.tar
## file prefix is MSM6J2MF
## zip_type is STD
## uid is MSM6J2MF
## file is MSM6J2MH.tar
## file prefix is MSM6J2MH
## zip_type is STD
## uid is MSM6J2MH
## file is MSM6J2MJ.tar
## file prefix is MSM6J2MJ
## zip_type is STD
## uid is MSM6J2MJ
## file is MSM6J2ML.tar
## file prefix is MSM6J2ML
## zip_type is STD
## uid is MSM6J2ML
## file is MSM79HBN.tar
## file prefix is MSM79HBN
## zip_type is STD
## uid is MSM79HBN
## file is MSM79HBP.tar
## file prefix is MSM79HBP
## zip_type is STD
## uid is MSM79HBP
## file is MSM79HBR.tar
## file prefix is MSM79HBR
## zip_type is STD
## uid is MSM79HBR
## file is MSM79HBT.tar
## file prefix is MSM79HBT
## zip_type is STD
## uid is MSM79HBT
## file is MSM79HBV.tar
## file prefix is MSM79HBV
## zip_type is STD
## uid is MSM79HBV
## file is SM-73BNG.tar
## file prefix is SM-73BNG
## zip_type is VIR
## uid is SM-73BNG
## file is SM-76EO3.tar
## file prefix is SM-76EO3
## zip_type is VIR
## uid is SM-76EO3
## file is SM-7E5EJ.tar
## file prefix is SM-7E5EJ
## zip_type is VIR
## uid is SM-7E5EJ
## file is SM-7IWER.tar
## file prefix is SM-7IWER
## zip_type is VIR
## uid is SM-7IWER
## file is SM-7KPV7.tar
## file prefix is SM-7KPV7
## zip_type is VIR
## uid is SM-7KPV7
## file is CSM6J2H9_P.fastq.gz
## file prefix is CSM6J2H9_P
## zip_type is P
## uid is CSM6J2H9_P
## file is MSM6J2HN.tar
## file prefix is MSM6J2HN
## zip_type is STD
## uid is MSM6J2HN
## file is MSM6J2HP.tar
## file prefix is MSM6J2HP
## zip_type is STD
## uid is MSM6J2HP
## file is MSM6J2HR.tar
## file prefix is MSM6J2HR
## zip_type is STD
## uid is MSM6J2HR
## file is MSM6J2HT.tar
## file prefix is MSM6J2HT
## zip_type is STD
## uid is MSM6J2HT
## file is MSM6J2QL.tar
## file prefix is MSM6J2QL
## zip_type is STD
## uid is MSM6J2QL
## file is MSM6J2QH.tar
## file prefix is MSM6J2QH
## zip_type is STD
## uid is MSM6J2QH
## file is MSM6J2QP.tar
## file prefix is MSM6J2QP
## zip_type is STD
## uid is MSM6J2QP
## file is MSM6J2QJ.tar
## file prefix is MSM6J2QJ
## zip_type is STD
## uid is MSM6J2QJ
## file is MSM6J2QF.tar
## file prefix is MSM6J2QF
## zip_type is STD
## uid is MSM6J2QF
## file is MSM79H5Q.tar
## file prefix is MSM79H5Q
## zip_type is STD
## uid is MSM79H5Q
## file is MSM79H5U.tar
## file prefix is MSM79H5U
## zip_type is STD
## uid is MSM79H5U
## file is MSM79H5Y.tar
## file prefix is MSM79H5Y
## zip_type is STD
## uid is MSM79H5Y
## file is MSM79H5S.tar
## file prefix is MSM79H5S
## zip_type is STD
## uid is MSM79H5S
## file is MSM79HAH.tar
## file prefix is MSM79HAH
## zip_type is STD
## uid is MSM79HAH
## file is SM-6XSKY.tar
## file prefix is SM-6XSKY
## zip_type is VIR
## uid is SM-6XSKY
## file is SM-6ZJI3.tar
## file prefix is SM-6ZJI3
## zip_type is VIR
## uid is SM-6ZJI3
## file is SM-71ZVB.tar
## file prefix is SM-71ZVB
## zip_type is VIR
## uid is SM-71ZVB
## file is SM-7B7LG.tar
## file prefix is SM-7B7LG
## zip_type is VIR
## uid is SM-7B7LG
## file is SM-7CP2G.tar
## file prefix is SM-7CP2G
## zip_type is VIR
## uid is SM-7CP2G
## file is SM-7GYIZ.tar
## file prefix is SM-7GYIZ
## zip_type is VIR
## uid is SM-7GYIZ
## file is SM-7M36N.tar
## file prefix is SM-7M36N
## zip_type is VIR
## uid is SM-7M36N
## file is SM-7RFY8.tar
## file prefix is SM-7RFY8
## zip_type is VIR
## uid is SM-7RFY8
## file is SM-9WYEE.tar
## file prefix is SM-9WYEE
## zip_type is VIR
## uid is SM-9WYEE
## file is MSM6J2JF_P.fastq.gz
## file prefix is MSM6J2JF_P
## zip_type is P
## uid is MSM6J2JF_P
## file is MSM6J2JN.tar
## file prefix is MSM6J2JN
## zip_type is STD
## uid is MSM6J2JN
## file is MSM6J2JP.tar
## file prefix is MSM6J2JP
## zip_type is STD
## uid is MSM6J2JP
## file is MSM6J2JR.tar
## file prefix is MSM6J2JR
## zip_type is STD
## uid is MSM6J2JR
## file is MSM6J2JT.tar
## file prefix is MSM6J2JT
## zip_type is STD
## uid is MSM6J2JT
## file is MSM6J2LV.tar
## file prefix is MSM6J2LV
## zip_type is STD
## uid is MSM6J2LV
## file is MSM6J2SE.tar
## file prefix is MSM6J2SE
## zip_type is STD
## uid is MSM6J2SE
## file is MSM6J2SI.tar
## file prefix is MSM6J2SI
## zip_type is STD
## uid is MSM6J2SI
## file is MSM6J2SK.tar
## file prefix is MSM6J2SK
## zip_type is STD
## uid is MSM6J2SK
## file is MSM79HES.tar
## file prefix is MSM79HES
## zip_type is STD
## uid is MSM79HES
## file is MSM79HEU.tar
## file prefix is MSM79HEU
## zip_type is STD
## uid is MSM79HEU
## file is MSM79HEW.tar
## file prefix is MSM79HEW
## zip_type is STD
## uid is MSM79HEW
## file is MSM79HEY.tar
## file prefix is MSM79HEY
## zip_type is STD
## uid is MSM79HEY
## file is MSM79HB6.tar
## file prefix is MSM79HB6
## zip_type is STD
## uid is MSM79HB6
## file is SM-6Y2UU.tar
## file prefix is SM-6Y2UU
## zip_type is VIR
## uid is SM-6Y2UU
## file is SM-6ZEV4.tar
## file prefix is SM-6ZEV4
## zip_type is VIR
## uid is SM-6ZEV4
## file is SM-72ZHH.tar
## file prefix is SM-72ZHH
## zip_type is VIR
## uid is SM-72ZHH
## file is SM-7CP39.tar
## file prefix is SM-7CP39
## zip_type is VIR
## uid is SM-7CP39
## file is SM-7ED9P.tar
## file prefix is SM-7ED9P
## zip_type is VIR
## uid is SM-7ED9P
## file is SM-7K25J.tar
## file prefix is SM-7K25J
## zip_type is VIR
## uid is SM-7K25J
## file is SM-9KONF.tar
## file prefix is SM-9KONF
## zip_type is VIR
## uid is SM-9KONF
## file is SM-9QMNB.tar
## file prefix is SM-9QMNB
## zip_type is VIR
## uid is SM-9QMNB
## file is SM-A312U.tar
## file prefix is SM-A312U
## zip_type is VIR
## uid is SM-A312U
## file is MSM6J2JH_P.fastq.gz
## file prefix is MSM6J2JH_P
## zip_type is P
## uid is MSM6J2JH_P
## file is MSM6J2JZ.tar
## file prefix is MSM6J2JZ
## zip_type is STD
## uid is MSM6J2JZ
## file is MSM6J2K2.tar
## file prefix is MSM6J2K2
## zip_type is STD
## uid is MSM6J2K2
## file is MSM6J2K4.tar
## file prefix is MSM6J2K4
## zip_type is STD
## uid is MSM6J2K4
## file is MSM6J2K6.tar
## file prefix is MSM6J2K6
## zip_type is STD
## uid is MSM6J2K6
## file is MSM6J2K8.tar
## file prefix is MSM6J2K8
## zip_type is STD
## uid is MSM6J2K8
## file is MSM6J2KA.tar
## file prefix is MSM6J2KA
## zip_type is STD
## uid is MSM6J2KA
## file is MSM6J2RK.tar
## file prefix is MSM6J2RK
## zip_type is STD
## uid is MSM6J2RK
## file is MSM6J2RM.tar
## file prefix is MSM6J2RM
## zip_type is STD
## uid is MSM6J2RM
## file is MSM6J2RO.tar
## file prefix is MSM6J2RO
## zip_type is STD
## uid is MSM6J2RO
## file is MSM6J2RQ.tar
## file prefix is MSM6J2RQ
## zip_type is STD
## uid is MSM6J2RQ
## file is MSM6J2RU.tar
## file prefix is MSM6J2RU
## zip_type is STD
## uid is MSM6J2RU
## file is MSM6J2RS.tar
## file prefix is MSM6J2RS
## zip_type is STD
## uid is MSM6J2RS
## file is MSM79H6D.tar
## file prefix is MSM79H6D
## zip_type is STD
## uid is MSM79H6D
## file is MSM79H6F.tar
## file prefix is MSM79H6F
## zip_type is STD
## uid is MSM79H6F
## file is MSM79H6J.tar
## file prefix is MSM79H6J
## zip_type is STD
## uid is MSM79H6J
## file is MSM79H6H.tar
## file prefix is MSM79H6H
## zip_type is STD
## uid is MSM79H6H
## file is MSM79H6L.tar
## file prefix is MSM79H6L
## zip_type is STD
## uid is MSM79H6L
## file is MSM79H6N.tar
## file prefix is MSM79H6N
## zip_type is STD
## uid is MSM79H6N
## file is MSM79H9Y.tar
## file prefix is MSM79H9Y
## zip_type is STD
## uid is MSM79H9Y
## file is MSM79HA1.tar
## file prefix is MSM79HA1
## zip_type is STD
## uid is MSM79HA1
## file is MSM79HA3.tar
## file prefix is MSM79HA3
## zip_type is STD
## uid is MSM79HA3
## file is MSM79HA7.tar
## file prefix is MSM79HA7
## zip_type is STD
## uid is MSM79HA7
## file is SM-6Y87K.tar
## file prefix is SM-6Y87K
## zip_type is VIR
## uid is SM-6Y87K
## file is SM-6ZEUL.tar
## file prefix is SM-6ZEUL
## zip_type is VIR
## uid is SM-6ZEUL
## file is SM-6ZT2P.tar
## file prefix is SM-6ZT2P
## zip_type is VIR
## uid is SM-6ZT2P
## file is SM-7CRWR.tar
## file prefix is SM-7CRWR
## zip_type is VIR
## uid is SM-7CRWR
## file is SM-7FY2O.tar
## file prefix is SM-7FY2O
## zip_type is VIR
## uid is SM-7FY2O
## file is SM-7MCTO.tar
## file prefix is SM-7MCTO
## zip_type is VIR
## uid is SM-7MCTO
## file is SM-9KOON.tar
## file prefix is SM-9KOON
## zip_type is VIR
## uid is SM-9KOON
## file is SM-9SNJ5.tar
## file prefix is SM-9SNJ5
## zip_type is VIR
## uid is SM-9SNJ5
## file is MSM6J2PO.tar
## file prefix is MSM6J2PO
## zip_type is STD
## uid is MSM6J2PO
## file is MSM6J2PQ_P.fastq.gz
## file prefix is MSM6J2PQ_P
## zip_type is P
## uid is MSM6J2PQ_P
## file is MSM6J2PS.tar
## file prefix is MSM6J2PS
## zip_type is STD
## uid is MSM6J2PS
## file is MSM6J2PU.tar
## file prefix is MSM6J2PU
## zip_type is STD
## uid is MSM6J2PU
## file is MSM6J2PW.tar
## file prefix is MSM6J2PW
## zip_type is STD
## uid is MSM6J2PW
## file is MSM79HC4.tar
## file prefix is MSM79HC4
## zip_type is STD
## uid is MSM79HC4
## file is MSM79HC8.tar
## file prefix is MSM79HC8
## zip_type is STD
## uid is MSM79HC8
## file is MSM79H7O.tar
## file prefix is MSM79H7O
## zip_type is STD
## uid is MSM79H7O
## file is MSM79H7Q.tar
## file prefix is MSM79H7Q
## zip_type is STD
## uid is MSM79H7Q
## file is MSM79H7W.tar
## file prefix is MSM79H7W
## zip_type is STD
## uid is MSM79H7W
## file is MSM79H7Y.tar
## file prefix is MSM79H7Y
## zip_type is STD
## uid is MSM79H7Y
## file is MSM9VZNH.tar
## file prefix is MSM9VZNH
## zip_type is STD
## uid is MSM9VZNH
## file is MSM9VZNH_TR.tar
## file prefix is MSM9VZNH_TR
## zip_type is TR
## uid is MSM9VZNH_TR
## file is MSM9VZNL.tar
## file prefix is MSM9VZNL
## zip_type is STD
## uid is MSM9VZNL
## file is SM-7GYI3.tar
## file prefix is SM-7GYI3
## zip_type is VIR
## uid is SM-7GYI3
## file is SM-7K1UP.tar
## file prefix is SM-7K1UP
## zip_type is VIR
## uid is SM-7K1UP
## file is SM-9JM8W.tar
## file prefix is SM-9JM8W
## zip_type is VIR
## uid is SM-9JM8W
## file is SM-9OS61.tar
## file prefix is SM-9OS61
## zip_type is VIR
## uid is SM-9OS61
## file is SM-9WJXI.tar
## file prefix is SM-9WJXI
## zip_type is VIR
## uid is SM-9WJXI
## file is SM-9Y7CX.tar
## file prefix is SM-9Y7CX
## zip_type is VIR
## uid is SM-9Y7CX
## file is SM-A5A1P.tar
## file prefix is SM-A5A1P
## zip_type is VIR
## uid is SM-A5A1P
## file is SM-A8GSW.tar
## file prefix is SM-A8GSW
## zip_type is VIR
## uid is SM-A8GSW
## file is SM-ADYKB.tar
## file prefix is SM-ADYKB
## zip_type is VIR
## uid is SM-ADYKB
## file is MSM6J2Q1.tar
## file prefix is MSM6J2Q1
## zip_type is STD
## uid is MSM6J2Q1
## file is MSM6J2LW.tar
## file prefix is MSM6J2LW
## zip_type is STD
## uid is MSM6J2LW
## file is MSM6J2LY.tar
## file prefix is MSM6J2LY
## zip_type is STD
## uid is MSM6J2LY
## file is MSM6J2M3_P.fastq.gz
## file prefix is MSM6J2M3_P
## zip_type is P
## uid is MSM6J2M3_P
## file is MSM6J2M3.tar
## file prefix is MSM6J2M3
## zip_type is STD
## uid is MSM6J2M3
## file is MSM79HBL.tar
## file prefix is MSM79HBL
## zip_type is STD
## uid is MSM79HBL
## file is MSM79HFW.tar
## file prefix is MSM79HFW
## zip_type is STD
## uid is MSM79HFW
## file is MSM79HFY.tar
## file prefix is MSM79HFY
## zip_type is STD
## uid is MSM79HFY
## file is MSM79H6Y.tar
## file prefix is MSM79H6Y
## zip_type is STD
## uid is MSM79H6Y
## file is MSM79HBB.tar
## file prefix is MSM79HBB
## zip_type is STD
## uid is MSM79HBB
## file is SM-787I5.tar
## file prefix is SM-787I5
## zip_type is VIR
## uid is SM-787I5
## file is SM-7GCHD.tar
## file prefix is SM-7GCHD
## zip_type is VIR
## uid is SM-7GCHD
## file is SM-9P5G4.tar
## file prefix is SM-9P5G4
## zip_type is VIR
## uid is SM-9P5G4
## file is SM-A375E.tar
## file prefix is SM-A375E
## zip_type is VIR
## uid is SM-A375E
## file is MSM6J2RG_P.fastq.gz
## file prefix is MSM6J2RG_P
## zip_type is P
## uid is MSM6J2RG_P
## file is MSM6J2N6_P.fastq.gz
## file prefix is MSM6J2N6_P
## zip_type is P
## uid is MSM6J2N6_P
## file is MSM79H52_P.fastq.gz
## file prefix is MSM79H52_P
## zip_type is P
## uid is MSM79H52_P
## file is MSM79H54.tar
## file prefix is MSM79H54
## zip_type is STD
## uid is MSM79H54
## file is MSM79H58.tar
## file prefix is MSM79H58
## zip_type is STD
## uid is MSM79H58
## file is MSM79H5A.tar
## file prefix is MSM79H5A
## zip_type is STD
## uid is MSM79H5A
## file is MSM79H9A.tar
## file prefix is MSM79H9A
## zip_type is STD
## uid is MSM79H9A
## file is MSM79H9C.tar
## file prefix is MSM79H9C
## zip_type is STD
## uid is MSM79H9C
## file is MSM79H9G.tar
## file prefix is MSM79H9G
## zip_type is STD
## uid is MSM79H9G
## file is MSM79H9K.tar
## file prefix is MSM79H9K
## zip_type is STD
## uid is MSM79H9K
## file is MSM9VZF7.tar
## file prefix is MSM9VZF7
## zip_type is STD
## uid is MSM9VZF7
## file is MSM9VZFF.tar
## file prefix is MSM9VZFF
## zip_type is STD
## uid is MSM9VZFF
## file is MSM9VZFH.tar
## file prefix is MSM9VZFH
## zip_type is STD
## uid is MSM9VZFH
## file is MSM9VZPT.tar
## file prefix is MSM9VZPT
## zip_type is STD
## uid is MSM9VZPT
## file is SM-7K56L.tar
## file prefix is SM-7K56L
## zip_type is VIR
## uid is SM-7K56L
## file is SM-9LB64.tar
## file prefix is SM-9LB64
## zip_type is VIR
## uid is SM-9LB64
## file is SM-9W3B1.tar
## file prefix is SM-9W3B1
## zip_type is VIR
## uid is SM-9W3B1
## file is SM-9WSQP.tar
## file prefix is SM-9WSQP
## zip_type is VIR
## uid is SM-9WSQP
## file is SM-A9F6V.tar
## file prefix is SM-A9F6V
## zip_type is VIR
## uid is SM-A9F6V
## file is MSM79H5E.tar
## file prefix is MSM79H5E
## zip_type is STD
## uid is MSM79H5E
## file is MSM79H5G.tar
## file prefix is MSM79H5G
## zip_type is STD
## uid is MSM79H5G
## file is MSM79H5K.tar
## file prefix is MSM79H5K
## zip_type is STD
## uid is MSM79H5K
## file is MSM79H5M.tar
## file prefix is MSM79H5M
## zip_type is STD
## uid is MSM79H5M
## file is MSM79H9M.tar
## file prefix is MSM79H9M
## zip_type is STD
## uid is MSM79H9M
## file is MSM79H9Q.tar
## file prefix is MSM79H9Q
## zip_type is STD
## uid is MSM79H9Q
## file is MSM79H9W.tar
## file prefix is MSM79H9W
## zip_type is STD
## uid is MSM79H9W
## file is MSM9VZGO.tar
## file prefix is MSM9VZGO
## zip_type is STD
## uid is MSM9VZGO
## file is MSM9VZGS.tar
## file prefix is MSM9VZGS
## zip_type is STD
## uid is MSM9VZGS
## file is MSM9VZGU.tar
## file prefix is MSM9VZGU
## zip_type is STD
## uid is MSM9VZGU
## file is MSM9VZLJ.tar
## file prefix is MSM9VZLJ
## zip_type is STD
## uid is MSM9VZLJ
## file is MSM9VZHB.tar
## file prefix is MSM9VZHB
## zip_type is STD
## uid is MSM9VZHB
## file is MSM9VZHF.tar
## file prefix is MSM9VZHF
## zip_type is STD
## uid is MSM9VZHF
## file is SM-7I1GQ.tar
## file prefix is SM-7I1GQ
## zip_type is VIR
## uid is SM-7I1GQ
## file is SM-7K25B.tar
## file prefix is SM-7K25B
## zip_type is VIR
## uid is SM-7K25B
## file is SM-7MWXI.tar
## file prefix is SM-7MWXI
## zip_type is VIR
## uid is SM-7MWXI
## file is SM-9MPH5.tar
## file prefix is SM-9MPH5
## zip_type is VIR
## uid is SM-9MPH5
## file is SM-9RRZJ.tar
## file prefix is SM-9RRZJ
## zip_type is VIR
## uid is SM-9RRZJ
## file is SM-9Y7DZ.tar
## file prefix is SM-9Y7DZ
## zip_type is VIR
## uid is SM-9Y7DZ
## file is SM-A4UVS.tar
## file prefix is SM-A4UVS
## zip_type is VIR
## uid is SM-A4UVS
## file is SM-AF6NC.tar
## file prefix is SM-AF6NC
## zip_type is VIR
## uid is SM-AF6NC
## file is SM-AMR38.tar
## file prefix is SM-AMR38
## zip_type is VIR
## uid is SM-AMR38
## file is SM-ARMDF.tar
## file prefix is SM-ARMDF
## zip_type is VIR
## uid is SM-ARMDF
## file is MSM79HF1_P.fastq.gz
## file prefix is MSM79HF1_P
## zip_type is P
## uid is MSM79HF1_P
## file is MSM79HF3.tar
## file prefix is MSM79HF3
## zip_type is STD
## uid is MSM79HF3
## file is MSM79HF5.tar
## file prefix is MSM79HF5
## zip_type is STD
## uid is MSM79HF5
## file is MSM79HF7.tar
## file prefix is MSM79HF7
## zip_type is STD
## uid is MSM79HF7
## file is MSM79HF9.tar
## file prefix is MSM79HF9
## zip_type is STD
## uid is MSM79HF9
## file is MSM79HF9_TR.tar
## file prefix is MSM79HF9_TR
## zip_type is TR
## uid is MSM79HF9_TR
## file is MSM79HFB.tar
## file prefix is MSM79HFB
## zip_type is STD
## uid is MSM79HFB
## file is MSM79HAJ.tar
## file prefix is MSM79HAJ
## zip_type is STD
## uid is MSM79HAJ
## file is MSM79HAL.tar
## file prefix is MSM79HAL
## zip_type is STD
## uid is MSM79HAL
## file is MSM79HAN.tar
## file prefix is MSM79HAN
## zip_type is STD
## uid is MSM79HAN
## file is MSM79HAT.tar
## file prefix is MSM79HAT
## zip_type is STD
## uid is MSM79HAT
## file is MSM79HAR.tar
## file prefix is MSM79HAR
## zip_type is STD
## uid is MSM79HAR
## file is MSM9VZLL.tar
## file prefix is MSM9VZLL
## zip_type is STD
## uid is MSM9VZLL
## file is MSM9VZLN.tar
## file prefix is MSM9VZLN
## zip_type is STD
## uid is MSM9VZLN
## file is MSM9VZLP.tar
## file prefix is MSM9VZLP
## zip_type is STD
## uid is MSM9VZLP
## file is MSM9VZLR.tar
## file prefix is MSM9VZLR
## zip_type is STD
## uid is MSM9VZLR
## file is MSM9VZLV.tar
## file prefix is MSM9VZLV
## zip_type is STD
## uid is MSM9VZLV
## file is MSM9VZLT.tar
## file prefix is MSM9VZLT
## zip_type is STD
## uid is MSM9VZLT
## file is MSMA26ER.tar
## file prefix is MSMA26ER
## zip_type is STD
## uid is MSMA26ER
## file is MSMA26ET.tar
## file prefix is MSMA26ET
## zip_type is STD
## uid is MSMA26ET
## file is MSMA26EZ.tar
## file prefix is MSMA26EZ
## zip_type is STD
## uid is MSMA26EZ
## file is MSMA2688.tar
## file prefix is MSMA2688
## zip_type is STD
## uid is MSMA2688
## file is SM-7RYA9.tar
## file prefix is SM-7RYA9
## zip_type is VIR
## uid is SM-7RYA9
## file is SM-9LISD.tar
## file prefix is SM-9LISD
## zip_type is VIR
## uid is SM-9LISD
## file is SM-9PJ15.tar
## file prefix is SM-9PJ15
## zip_type is VIR
## uid is SM-9PJ15
## file is SM-9SNKK.tar
## file prefix is SM-9SNKK
## zip_type is VIR
## uid is SM-9SNKK
## file is SM-A9F5L.tar
## file prefix is SM-A9F5L
## zip_type is VIR
## uid is SM-A9F5L
## file is SM-A9VSD.tar
## file prefix is SM-A9VSD
## zip_type is VIR
## uid is SM-A9VSD
## file is SM-AGUN4.tar
## file prefix is SM-AGUN4
## zip_type is VIR
## uid is SM-AGUN4
## file is SM-AKTZ4.tar
## file prefix is SM-AKTZ4
## zip_type is VIR
## uid is SM-AKTZ4
## file is SM-AW573.tar
## file prefix is SM-AW573
## zip_type is VIR
## uid is SM-AW573
## file is MSM79H9U_P.fastq.gz
## file prefix is MSM79H9U_P
## zip_type is P
## uid is MSM79H9U_P
## file is MSM79HBX_P.fastq.gz
## file prefix is MSM79HBX_P
## zip_type is P
## uid is MSM79HBX_P
## file is MSM79HDK.tar
## file prefix is MSM79HDK
## zip_type is STD
## uid is MSM79HDK
## file is MSM79HDM.tar
## file prefix is MSM79HDM
## zip_type is STD
## uid is MSM79HDM
## file is MSM79HDO.tar
## file prefix is MSM79HDO
## zip_type is STD
## uid is MSM79HDO
## file is MSM79HDQ.tar
## file prefix is MSM79HDQ
## zip_type is STD
## uid is MSM79HDQ
## file is MSM79HDQ_TR.tar
## file prefix is MSM79HDQ_TR
## zip_type is TR
## uid is MSM79HDQ_TR
## file is MSM79HDS.tar
## file prefix is MSM79HDS
## zip_type is STD
## uid is MSM79HDS
## file is MSM79HDU.tar
## file prefix is MSM79HDU
## zip_type is STD
## uid is MSM79HDU
## file is MSM79H98.tar
## file prefix is MSM79H98
## zip_type is STD
## uid is MSM79H98
## file is MSM9VZEK.tar
## file prefix is MSM9VZEK
## zip_type is STD
## uid is MSM9VZEK
## file is MSM9VZEK_TR.tar
## file prefix is MSM9VZEK_TR
## zip_type is TR
## uid is MSM9VZEK_TR
## file is MSM9VZEM.tar
## file prefix is MSM9VZEM
## zip_type is STD
## uid is MSM9VZEM
## file is MSM9VZEO.tar
## file prefix is MSM9VZEO
## zip_type is STD
## uid is MSM9VZEO
## file is MSM9VZEQ.tar
## file prefix is MSM9VZEQ
## zip_type is STD
## uid is MSM9VZEQ
## file is MSM9VZES.tar
## file prefix is MSM9VZES
## zip_type is STD
## uid is MSM9VZES
## file is MSM7J16J.tar
## file prefix is MSM7J16J
## zip_type is STD
## uid is MSM7J16J
## file is MSM7J16L.tar
## file prefix is MSM7J16L
## zip_type is STD
## uid is MSM7J16L
## file is MSM7J16N.tar
## file prefix is MSM7J16N
## zip_type is STD
## uid is MSM7J16N
## file is MSM7J16P.tar
## file prefix is MSM7J16P
## zip_type is STD
## uid is MSM7J16P
## file is MSM7J16R.tar
## file prefix is MSM7J16R
## zip_type is STD
## uid is MSM7J16R
## file is MSM9VZJB.tar
## file prefix is MSM9VZJB
## zip_type is STD
## uid is MSM9VZJB
## file is MSMA26AL.tar
## file prefix is MSMA26AL
## zip_type is STD
## uid is MSMA26AL
## file is MSMA26AN.tar
## file prefix is MSMA26AN
## zip_type is STD
## uid is MSMA26AN
## file is MSMA26AP.tar
## file prefix is MSMA26AP
## zip_type is STD
## uid is MSMA26AP
## file is MSMA26AR.tar
## file prefix is MSMA26AR
## zip_type is STD
## uid is MSMA26AR
## file is MSMA26AT.tar
## file prefix is MSMA26AT
## zip_type is STD
## uid is MSMA26AT
## file is SM-9RTUY.tar
## file prefix is SM-9RTUY
## zip_type is VIR
## uid is SM-9RTUY
## file is SM-A12KK.tar
## file prefix is SM-A12KK
## zip_type is VIR
## uid is SM-A12KK
## file is SM-AADKV.tar
## file prefix is SM-AADKV
## zip_type is VIR
## uid is SM-AADKV
## file is SM-AEB6F.tar
## file prefix is SM-AEB6F
## zip_type is VIR
## uid is SM-AEB6F
## file is SM-AOMTB.tar
## file prefix is SM-AOMTB
## zip_type is VIR
## uid is SM-AOMTB
## file is MSM79HBZ.tar
## file prefix is MSM79HBZ
## zip_type is STD
## uid is MSM79HBZ
## file is MSM79HD8_P.fastq.gz
## file prefix is MSM79HD8_P
## zip_type is P
## uid is MSM79HD8_P
## file is MSM79HDA.tar
## file prefix is MSM79HDA
## zip_type is STD
## uid is MSM79HDA
## file is MSM79HDC.tar
## file prefix is MSM79HDC
## zip_type is STD
## uid is MSM79HDC
## file is MSM79HDE.tar
## file prefix is MSM79HDE
## zip_type is STD
## uid is MSM79HDE
## file is MSM79HDG.tar
## file prefix is MSM79HDG
## zip_type is STD
## uid is MSM79HDG
## file is MSM79HDG_TR.tar
## file prefix is MSM79HDG_TR
## zip_type is TR
## uid is MSM79HDG_TR
## file is MSM79HDI.tar
## file prefix is MSM79HDI
## zip_type is STD
## uid is MSM79HDI
## file is MSM9VZEU.tar
## file prefix is MSM9VZEU
## zip_type is STD
## uid is MSM9VZEU
## file is MSM9VZEW.tar
## file prefix is MSM9VZEW
## zip_type is STD
## uid is MSM9VZEW
## file is MSM9VZEY.tar
## file prefix is MSM9VZEY
## zip_type is STD
## uid is MSM9VZEY
## file is MSM9VZF1.tar
## file prefix is MSM9VZF1
## zip_type is STD
## uid is MSM9VZF1
## file is MSM9VZF3.tar
## file prefix is MSM9VZF3
## zip_type is STD
## uid is MSM9VZF3
## file is MSM9VZF5.tar
## file prefix is MSM9VZF5
## zip_type is STD
## uid is MSM9VZF5
## file is MSM9VZOY.tar
## file prefix is MSM9VZOY
## zip_type is STD
## uid is MSM9VZOY
## file is MSM9VZOU.tar
## file prefix is MSM9VZOU
## zip_type is STD
## uid is MSM9VZOU
## file is MSM9VZOW.tar
## file prefix is MSM9VZOW
## zip_type is STD
## uid is MSM9VZOW
## file is MSM9VZOS.tar
## file prefix is MSM9VZOS
## zip_type is STD
## uid is MSM9VZOS
## file is MSM9VZP1.tar
## file prefix is MSM9VZP1
## zip_type is STD
## uid is MSM9VZP1
## file is MSM9VZP3.tar
## file prefix is MSM9VZP3
## zip_type is STD
## uid is MSM9VZP3
## file is MSMA26EH.tar
## file prefix is MSMA26EH
## zip_type is STD
## uid is MSMA26EH
## file is MSMA26EJ.tar
## file prefix is MSMA26EJ
## zip_type is STD
## uid is MSMA26EJ
## file is MSMA26EL.tar
## file prefix is MSMA26EL
## zip_type is STD
## uid is MSMA26EL
## file is MSMA26EN.tar
## file prefix is MSMA26EN
## zip_type is STD
## uid is MSMA26EN
## file is MSMA26EP.tar
## file prefix is MSMA26EP
## zip_type is STD
## uid is MSMA26EP
## file is SM-9NBEM.tar
## file prefix is SM-9NBEM
## zip_type is VIR
## uid is SM-9NBEM
## file is SM-9QMOE.tar
## file prefix is SM-9QMOE
## zip_type is VIR
## uid is SM-9QMOE
## file is SM-9SIJC.tar
## file prefix is SM-9SIJC
## zip_type is VIR
## uid is SM-9SIJC
## file is SM-9VW91.tar
## file prefix is SM-9VW91
## zip_type is VIR
## uid is SM-9VW91
## file is SM-9XG8C.tar
## file prefix is SM-9XG8C
## zip_type is VIR
## uid is SM-9XG8C
## file is SM-9ZA1F.tar
## file prefix is SM-9ZA1F
## zip_type is VIR
## uid is SM-9ZA1F
## file is SM-A1ZVI.tar
## file prefix is SM-A1ZVI
## zip_type is VIR
## uid is SM-A1ZVI
## file is SM-A4GB5.tar
## file prefix is SM-A4GB5
## zip_type is VIR
## uid is SM-A4GB5
## file is SM-A77UP.tar
## file prefix is SM-A77UP
## zip_type is VIR
## uid is SM-A77UP
## file is SM-A9F63.tar
## file prefix is SM-A9F63
## zip_type is VIR
## uid is SM-A9F63
## file is SM-AF6MN.tar
## file prefix is SM-AF6MN
## zip_type is VIR
## uid is SM-AF6MN
## file is SM-AGUNF.tar
## file prefix is SM-AGUNF
## zip_type is VIR
## uid is SM-AGUNF
## file is SM-AIG6M.tar
## file prefix is SM-AIG6M
## zip_type is VIR
## uid is SM-AIG6M
## file is SM-AKTYV.tar
## file prefix is SM-AKTYV
## zip_type is VIR
## uid is SM-AKTYV
## file is SM-AN8MO.tar
## file prefix is SM-AN8MO
## zip_type is VIR
## uid is SM-AN8MO
## file is SM-APDY6.tar
## file prefix is SM-APDY6
## zip_type is VIR
## uid is SM-APDY6
## file is SM-ARMBH.tar
## file prefix is SM-ARMBH
## zip_type is VIR
## uid is SM-ARMBH
## file is SM-AVACV.tar
## file prefix is SM-AVACV
## zip_type is VIR
## uid is SM-AVACV
## file is SM-AXQRB.tar
## file prefix is SM-AXQRB
## zip_type is VIR
## uid is SM-AXQRB
## file is MSM79HD6_P.fastq.gz
## file prefix is MSM79HD6_P
## zip_type is P
## uid is MSM79HD6_P
## file is MSM79H7C.tar
## file prefix is MSM79H7C
## zip_type is STD
## uid is MSM79H7C
## file is MSM79H7E.tar
## file prefix is MSM79H7E
## zip_type is STD
## uid is MSM79H7E
## file is MSM79H7M.tar
## file prefix is MSM79H7M
## zip_type is STD
## uid is MSM79H7M
## file is MSM79H7G.tar
## file prefix is MSM79H7G
## zip_type is STD
## uid is MSM79H7G
## file is MSMA26AX.tar
## file prefix is MSMA26AX
## zip_type is STD
## uid is MSMA26AX
## file is MSMA26AZ.tar
## file prefix is MSMA26AZ
## zip_type is STD
## uid is MSMA26AZ
## file is MSMA26AZ_TR.tar
## file prefix is MSMA26AZ_TR
## zip_type is TR
## uid is MSMA26AZ_TR
## file is MSMA26AV.tar
## file prefix is MSMA26AV
## zip_type is STD
## uid is MSMA26AV
## file is MSMB4LZ4.tar
## file prefix is MSMB4LZ4
## zip_type is STD
## uid is MSMB4LZ4
## file is SM-9VWCC.tar
## file prefix is SM-9VWCC
## zip_type is VIR
## uid is SM-9VWCC
## file is SM-9X47T.tar
## file prefix is SM-9X47T
## zip_type is VIR
## uid is SM-9X47T
## file is SM-AAHYC.tar
## file prefix is SM-AAHYC
## zip_type is VIR
## uid is SM-AAHYC
## file is SM-AUNQ4.tar
## file prefix is SM-AUNQ4
## zip_type is VIR
## uid is SM-AUNQ4
## file is SM-AYIRK.tar
## file prefix is SM-AYIRK
## zip_type is VIR
## uid is SM-AYIRK
## file is SM-CMUKM.tar
## file prefix is SM-CMUKM
## zip_type is VIR
## uid is SM-CMUKM
## file is MSM79HCG.tar
## file prefix is MSM79HCG
## zip_type is STD
## uid is MSM79HCG
## file is MSM79HCI.tar
## file prefix is MSM79HCI
## zip_type is STD
## uid is MSM79HCI
## file is MSM79HCK.tar
## file prefix is MSM79HCK
## zip_type is STD
## uid is MSM79HCK
## file is MSM79HCN_P.fastq.gz
## file prefix is MSM79HCN_P
## zip_type is P
## uid is MSM79HCN_P
## file is MSM79HCP.tar
## file prefix is MSM79HCP
## zip_type is STD
## uid is MSM79HCP
## file is MSM79HCR.tar
## file prefix is MSM79HCR
## zip_type is STD
## uid is MSM79HCR
## file is MSM79H81.tar
## file prefix is MSM79H81
## zip_type is STD
## uid is MSM79H81
## file is MSM79H83.tar
## file prefix is MSM79H83
## zip_type is STD
## uid is MSM79H83
## file is MSM79H85.tar
## file prefix is MSM79H85
## zip_type is STD
## uid is MSM79H85
## file is MSM79H87.tar
## file prefix is MSM79H87
## zip_type is STD
## uid is MSM79H87
## file is MSM79H89.tar
## file prefix is MSM79H89
## zip_type is STD
## uid is MSM79H89
## file is MSM79H8B.tar
## file prefix is MSM79H8B
## zip_type is STD
## uid is MSM79H8B
## file is MSM9VZOG.tar
## file prefix is MSM9VZOG
## zip_type is STD
## uid is MSM9VZOG
## file is MSM9VZOI.tar
## file prefix is MSM9VZOI
## zip_type is STD
## uid is MSM9VZOI
## file is MSM9VZOK.tar
## file prefix is MSM9VZOK
## zip_type is STD
## uid is MSM9VZOK
## file is MSM9VZOM.tar
## file prefix is MSM9VZOM
## zip_type is STD
## uid is MSM9VZOM
## file is MSM9VZOO.tar
## file prefix is MSM9VZOO
## zip_type is STD
## uid is MSM9VZOO
## file is MSM9VZOQ.tar
## file prefix is MSM9VZOQ
## zip_type is STD
## uid is MSM9VZOQ
## file is MSM9VZHJ.tar
## file prefix is MSM9VZHJ
## zip_type is STD
## uid is MSM9VZHJ
## file is MSM9VZHL.tar
## file prefix is MSM9VZHL
## zip_type is STD
## uid is MSM9VZHL
## file is MSM9VZHN.tar
## file prefix is MSM9VZHN
## zip_type is STD
## uid is MSM9VZHN
## file is MSM9VZHP.tar
## file prefix is MSM9VZHP
## zip_type is STD
## uid is MSM9VZHP
## file is MSM9VZHR.tar
## file prefix is MSM9VZHR
## zip_type is STD
## uid is MSM9VZHR
## file is MSM9VZHT.tar
## file prefix is MSM9VZHT
## zip_type is STD
## uid is MSM9VZHT
## file is SM-9N78Z.tar
## file prefix is SM-9N78Z
## zip_type is VIR
## uid is SM-9N78Z
## file is SM-9RRZ6.tar
## file prefix is SM-9RRZ6
## zip_type is VIR
## uid is SM-9RRZ6
## file is SM-9U25W.tar
## file prefix is SM-9U25W
## zip_type is VIR
## uid is SM-9U25W
## file is SM-9WO1B.tar
## file prefix is SM-9WO1B
## zip_type is VIR
## uid is SM-9WO1B
## file is SM-9Y491.tar
## file prefix is SM-9Y491
## zip_type is VIR
## uid is SM-9Y491
## file is SM-A373I.tar
## file prefix is SM-A373I
## zip_type is VIR
## uid is SM-A373I
## file is SM-ACTMW.tar
## file prefix is SM-ACTMW
## zip_type is VIR
## uid is SM-ACTMW
## file is SM-AFSHO.tar
## file prefix is SM-AFSHO
## zip_type is VIR
## uid is SM-AFSHO
## file is SM-ALM8R.tar
## file prefix is SM-ALM8R
## zip_type is VIR
## uid is SM-ALM8R
## file is SM-ANJ7Y.tar
## file prefix is SM-ANJ7Y
## zip_type is VIR
## uid is SM-ANJ7Y
## file is SM-ARHO4.tar
## file prefix is SM-ARHO4
## zip_type is VIR
## uid is SM-ARHO4
## file is SM-AWQJ7.tar
## file prefix is SM-AWQJ7
## zip_type is VIR
## uid is SM-AWQJ7
## file is MSM79H8D.tar
## file prefix is MSM79H8D
## zip_type is STD
## uid is MSM79H8D
## file is MSM79H8F.tar
## file prefix is MSM79H8F
## zip_type is STD
## uid is MSM79H8F
## file is MSM79H8H.tar
## file prefix is MSM79H8H
## zip_type is STD
## uid is MSM79H8H
## file is MSM79H8J_P.fastq.gz
## file prefix is MSM79H8J_P
## zip_type is P
## uid is MSM79H8J_P
## file is MSM79H8L.tar
## file prefix is MSM79H8L
## zip_type is STD
## uid is MSM79H8L
## file is MSM79H8N.tar
## file prefix is MSM79H8N
## zip_type is STD
## uid is MSM79H8N
## file is MSM9VZJZ.tar
## file prefix is MSM9VZJZ
## zip_type is STD
## uid is MSM9VZJZ
## file is MSM9VZGW.tar
## file prefix is MSM9VZGW
## zip_type is STD
## uid is MSM9VZGW
## file is MSM9VZGY.tar
## file prefix is MSM9VZGY
## zip_type is STD
## uid is MSM9VZGY
## file is MSM9VZH7.tar
## file prefix is MSM9VZH7
## zip_type is STD
## uid is MSM9VZH7
## file is MSMAPC7J.tar
## file prefix is MSMAPC7J
## zip_type is STD
## uid is MSMAPC7J
## file is SM-9ZJI9.tar
## file prefix is SM-9ZJI9
## zip_type is VIR
## uid is SM-9ZJI9
## file is SM-A9F61.tar
## file prefix is SM-A9F61
## zip_type is VIR
## uid is SM-A9F61
## file is SM-AI38I.tar
## file prefix is SM-AI38I
## zip_type is VIR
## uid is SM-AI38I
## file is SM-B5F8Z.tar
## file prefix is SM-B5F8Z
## zip_type is VIR
## uid is SM-B5F8Z
## file is MSM79H94_P.fastq.gz
## file prefix is MSM79H94_P
## zip_type is P
## uid is MSM79H94_P
## file is MSM9VZMA.tar
## file prefix is MSM9VZMA
## zip_type is STD
## uid is MSM9VZMA
## file is MSM9VZMA_TR.tar
## file prefix is MSM9VZMA_TR
## zip_type is TR
## uid is MSM9VZMA_TR
## file is MSM9VZME.tar
## file prefix is MSM9VZME
## zip_type is STD
## uid is MSM9VZME
## file is MSM9VZMC.tar
## file prefix is MSM9VZMC
## zip_type is STD
## uid is MSM9VZMC
## file is MSM9VZMI.tar
## file prefix is MSM9VZMI
## zip_type is STD
## uid is MSM9VZMI
## file is MSM9VZHX.tar
## file prefix is MSM9VZHX
## zip_type is STD
## uid is MSM9VZHX
## file is MSM9VZHZ.tar
## file prefix is MSM9VZHZ
## zip_type is STD
## uid is MSM9VZHZ
## file is MSM9VZI2.tar
## file prefix is MSM9VZI2
## zip_type is STD
## uid is MSM9VZI2
## file is MSM9VZI6.tar
## file prefix is MSM9VZI6
## zip_type is STD
## uid is MSM9VZI6
## file is MSMAPC7P.tar
## file prefix is MSMAPC7P
## zip_type is STD
## uid is MSMAPC7P
## file is MSMAPC7R.tar
## file prefix is MSMAPC7R
## zip_type is STD
## uid is MSMAPC7R
## file is MSMAPC7T.tar
## file prefix is MSMAPC7T
## zip_type is STD
## uid is MSMAPC7T
## file is MSMB4LZC.tar
## file prefix is MSMB4LZC
## zip_type is STD
## uid is MSMB4LZC
## file is SM-A8GSK.tar
## file prefix is SM-A8GSK
## zip_type is VIR
## uid is SM-A8GSK
## file is SM-AK8KQ.tar
## file prefix is SM-AK8KQ
## zip_type is VIR
## uid is SM-AK8KQ
## file is SM-AQJXU.tar
## file prefix is SM-AQJXU
## zip_type is VIR
## uid is SM-AQJXU
## file is SM-AU2MI.tar
## file prefix is SM-AU2MI
## zip_type is VIR
## uid is SM-AU2MI
## file is SM-AYMHU.tar
## file prefix is SM-AYMHU
## zip_type is VIR
## uid is SM-AYMHU
## file is SM-BYB94.tar
## file prefix is SM-BYB94
## zip_type is VIR
## uid is SM-BYB94
## file is SM-BZ5M7.tar
## file prefix is SM-BZ5M7
## zip_type is VIR
## uid is SM-BZ5M7
## file is SM-CNW2F.tar
## file prefix is SM-CNW2F
## zip_type is VIR
## uid is SM-CNW2F
## file is MSM9VZFJ_P.fastq.gz
## file prefix is MSM9VZFJ_P
## zip_type is P
## uid is MSM9VZFJ_P
## file is MSM9VZFL.tar
## file prefix is MSM9VZFL
## zip_type is STD
## uid is MSM9VZFL
## file is MSM9VZFN.tar
## file prefix is MSM9VZFN
## zip_type is STD
## uid is MSM9VZFN
## file is MSM9VZFR.tar
## file prefix is MSM9VZFR
## zip_type is STD
## uid is MSM9VZFR
## file is MSM9VZFT.tar
## file prefix is MSM9VZFT
## zip_type is STD
## uid is MSM9VZFT
## file is MSM9VZKC.tar
## file prefix is MSM9VZKC
## zip_type is STD
## uid is MSM9VZKC
## file is MSM9VZKE.tar
## file prefix is MSM9VZKE
## zip_type is STD
## uid is MSM9VZKE
## file is MSM9VZKI.tar
## file prefix is MSM9VZKI
## zip_type is STD
## uid is MSM9VZKI
## file is MSMA26BN.tar
## file prefix is MSMA26BN
## zip_type is STD
## uid is MSMA26BN
## file is MSMA26BR.tar
## file prefix is MSMA26BR
## zip_type is STD
## uid is MSMA26BR
## file is MSMA26BT.tar
## file prefix is MSMA26BT
## zip_type is STD
## uid is MSMA26BT
## file is MSMA26BV.tar
## file prefix is MSMA26BV
## zip_type is STD
## uid is MSMA26BV
## file is MSMA26BX.tar
## file prefix is MSMA26BX
## zip_type is STD
## uid is MSMA26BX
## file is MSMB4LZ8.tar
## file prefix is MSMB4LZ8
## zip_type is STD
## uid is MSMB4LZ8
## file is SM-A1ZV6.tar
## file prefix is SM-A1ZV6
## zip_type is VIR
## uid is SM-A1ZV6
## file is SM-A9F52.tar
## file prefix is SM-A9F52
## zip_type is VIR
## uid is SM-A9F52
## file is SM-AIG6I.tar
## file prefix is SM-AIG6I
## zip_type is VIR
## uid is SM-AIG6I
## file is SM-AKTYG.tar
## file prefix is SM-AKTYG
## zip_type is VIR
## uid is SM-AKTYG
## file is SM-APDX2.tar
## file prefix is SM-APDX2
## zip_type is VIR
## uid is SM-APDX2
## file is SM-AXTT2.tar
## file prefix is SM-AXTT2
## zip_type is VIR
## uid is SM-AXTT2
## file is SM-B5F98.tar
## file prefix is SM-B5F98
## zip_type is VIR
## uid is SM-B5F98
## file is SM-CP6EN.tar
## file prefix is SM-CP6EN
## zip_type is VIR
## uid is SM-CP6EN
## file is MSM9VZN4_P.fastq.gz
## file prefix is MSM9VZN4_P
## zip_type is P
## uid is MSM9VZN4_P
## file is MSM9VZLX_P.fastq.gz
## file prefix is MSM9VZLX_P
## zip_type is P
## uid is MSM9VZLX_P
## file is MSM9VZLZ.tar
## file prefix is MSM9VZLZ
## zip_type is STD
## uid is MSM9VZLZ
## file is MSM9VZNR.tar
## file prefix is MSM9VZNR
## zip_type is STD
## uid is MSM9VZNR
## file is MSM9VZNX.tar
## file prefix is MSM9VZNX
## zip_type is STD
## uid is MSM9VZNX
## file is MSM9VZNZ.tar
## file prefix is MSM9VZNZ
## zip_type is STD
## uid is MSM9VZNZ
## file is MSM9VZO2.tar
## file prefix is MSM9VZO2
## zip_type is STD
## uid is MSM9VZO2
## file is MSM9VZIM.tar
## file prefix is MSM9VZIM
## zip_type is STD
## uid is MSM9VZIM
## file is MSM9VZIO.tar
## file prefix is MSM9VZIO
## zip_type is STD
## uid is MSM9VZIO
## file is MSM9VZIQ.tar
## file prefix is MSM9VZIQ
## zip_type is STD
## uid is MSM9VZIQ
## file is MSM9VZIS.tar
## file prefix is MSM9VZIS
## zip_type is STD
## uid is MSM9VZIS
## file is MSM9VZIU.tar
## file prefix is MSM9VZIU
## zip_type is STD
## uid is MSM9VZIU
## file is MSM9VZIW.tar
## file prefix is MSM9VZIW
## zip_type is STD
## uid is MSM9VZIW
## file is MSMAPC55.tar
## file prefix is MSMAPC55
## zip_type is STD
## uid is MSMAPC55
## file is MSMAPC57.tar
## file prefix is MSMAPC57
## zip_type is STD
## uid is MSMAPC57
## file is MSMAPC59.tar
## file prefix is MSMAPC59
## zip_type is STD
## uid is MSMAPC59
## file is MSMAPC5B.tar
## file prefix is MSMAPC5B
## zip_type is STD
## uid is MSMAPC5B
## file is MSMAPC5D.tar
## file prefix is MSMAPC5D
## zip_type is STD
## uid is MSMAPC5D
## file is MSMB4LZP.tar
## file prefix is MSMB4LZP
## zip_type is STD
## uid is MSMB4LZP
## file is MSMB4LZR.tar
## file prefix is MSMB4LZR
## zip_type is STD
## uid is MSMB4LZR
## file is SM-AFSI1.tar
## file prefix is SM-AFSI1
## zip_type is VIR
## uid is SM-AFSI1
## file is SM-AN5G2.tar
## file prefix is SM-AN5G2
## zip_type is VIR
## uid is SM-AN5G2
## file is SM-ARMC1.tar
## file prefix is SM-ARMC1
## zip_type is VIR
## uid is SM-ARMC1
## file is SM-BZNKE.tar
## file prefix is SM-BZNKE
## zip_type is VIR
## uid is SM-BZNKE
## file is SM-CKJ6U.tar
## file prefix is SM-CKJ6U
## zip_type is VIR
## uid is SM-CKJ6U
## file is SM-COSYA.tar
## file prefix is SM-COSYA
## zip_type is VIR
## uid is SM-COSYA
## file is MSM9VZMM.tar
## file prefix is MSM9VZMM
## zip_type is STD
## uid is MSM9VZMM
## file is MSM9VZMO.tar
## file prefix is MSM9VZMO
## zip_type is STD
## uid is MSM9VZMO
## file is MSM9VZMS.tar
## file prefix is MSM9VZMS
## zip_type is STD
## uid is MSM9VZMS
## file is MSM9VZMU.tar
## file prefix is MSM9VZMU
## zip_type is STD
## uid is MSM9VZMU
## file is MSM9VZMW.tar
## file prefix is MSM9VZMW
## zip_type is STD
## uid is MSM9VZMW
## file is MSM9VZL7.tar
## file prefix is MSM9VZL7
## zip_type is STD
## uid is MSM9VZL7
## file is MSM9VZL9.tar
## file prefix is MSM9VZL9
## zip_type is STD
## uid is MSM9VZL9
## file is MSM9VZLB.tar
## file prefix is MSM9VZLB
## zip_type is STD
## uid is MSM9VZLB
## file is MSM9VZLD.tar
## file prefix is MSM9VZLD
## zip_type is STD
## uid is MSM9VZLD
## file is MSM9VZLF.tar
## file prefix is MSM9VZLF
## zip_type is STD
## uid is MSM9VZLF
## file is MSM9VZLH.tar
## file prefix is MSM9VZLH
## zip_type is STD
## uid is MSM9VZLH
## file is MSMA26BB.tar
## file prefix is MSMA26BB
## zip_type is STD
## uid is MSMA26BB
## file is MSMA26BD.tar
## file prefix is MSMA26BD
## zip_type is STD
## uid is MSMA26BD
## file is MSMA26BF.tar
## file prefix is MSMA26BF
## zip_type is STD
## uid is MSMA26BF
## file is MSMA26BH.tar
## file prefix is MSMA26BH
## zip_type is STD
## uid is MSMA26BH
## file is MSMA26BJ.tar
## file prefix is MSMA26BJ
## zip_type is STD
## uid is MSMA26BJ
## file is MSMA26BL.tar
## file prefix is MSMA26BL
## zip_type is STD
## uid is MSMA26BL
## file is MSMAPC6K.tar
## file prefix is MSMAPC6K
## zip_type is STD
## uid is MSMAPC6K
## file is MSMAPC6M.tar
## file prefix is MSMAPC6M
## zip_type is STD
## uid is MSMAPC6M
## file is MSMAPC6O.tar
## file prefix is MSMAPC6O
## zip_type is STD
## uid is MSMAPC6O
## file is MSMB4LZZ.tar
## file prefix is MSMB4LZZ
## zip_type is STD
## uid is MSMB4LZZ
## file is MSMB4LZV.tar
## file prefix is MSMB4LZV
## zip_type is STD
## uid is MSMB4LZV
## file is MSMB4LZX.tar
## file prefix is MSMB4LZX
## zip_type is STD
## uid is MSMB4LZX
## file is SM-A47RR.tar
## file prefix is SM-A47RR
## zip_type is VIR
## uid is SM-A47RR
## file is SM-A77W1.tar
## file prefix is SM-A77W1
## zip_type is VIR
## uid is SM-A77W1
## file is SM-AF6MH.tar
## file prefix is SM-AF6MH
## zip_type is VIR
## uid is SM-AF6MH
## file is SM-AKTXR.tar
## file prefix is SM-AKTXR
## zip_type is VIR
## uid is SM-AKTXR
## file is SM-AMVF5.tar
## file prefix is SM-AMVF5
## zip_type is VIR
## uid is SM-AMVF5
## file is SM-BZWMY.tar
## file prefix is SM-BZWMY
## zip_type is VIR
## uid is SM-BZWMY
## file is SM-CMUKQ.tar
## file prefix is SM-CMUKQ
## zip_type is VIR
## uid is SM-CMUKQ
## file is MSM9VZM4_P.fastq.gz
## file prefix is MSM9VZM4_P
## zip_type is P
## uid is MSM9VZM4_P
## file is MSM9VZPH.tar
## file prefix is MSM9VZPH
## zip_type is STD
## uid is MSM9VZPH
## file is MSM9VZPL.tar
## file prefix is MSM9VZPL
## zip_type is STD
## uid is MSM9VZPL
## file is MSM9VZPN.tar
## file prefix is MSM9VZPN
## zip_type is STD
## uid is MSM9VZPN
## file is MSM9VZIY.tar
## file prefix is MSM9VZIY
## zip_type is STD
## uid is MSM9VZIY
## file is MSM9VZJ3.tar
## file prefix is MSM9VZJ3
## zip_type is STD
## uid is MSM9VZJ3
## file is MSMA26CV.tar
## file prefix is MSMA26CV
## zip_type is STD
## uid is MSMA26CV
## file is MSMA26CX.tar
## file prefix is MSMA26CX
## zip_type is STD
## uid is MSMA26CX
## file is MSMAPC6A.tar
## file prefix is MSMAPC6A
## zip_type is STD
## uid is MSMAPC6A
## file is MSMAPC6C.tar
## file prefix is MSMAPC6C
## zip_type is STD
## uid is MSMAPC6C
## file is MSMAPC6E.tar
## file prefix is MSMAPC6E
## zip_type is STD
## uid is MSMAPC6E
## file is MSMAPC6G.tar
## file prefix is MSMAPC6G
## zip_type is STD
## uid is MSMAPC6G
## file is MSMB4LXY.tar
## file prefix is MSMB4LXY
## zip_type is STD
## uid is MSMB4LXY
## file is MSMB4LYB.tar
## file prefix is MSMB4LYB
## zip_type is STD
## uid is MSMB4LYB
## file is SM-ACKUZ.tar
## file prefix is SM-ACKUZ
## zip_type is VIR
## uid is SM-ACKUZ
## file is SM-AG2E2.tar
## file prefix is SM-AG2E2
## zip_type is VIR
## uid is SM-AG2E2
## file is SM-ANU6G.tar
## file prefix is SM-ANU6G
## zip_type is VIR
## uid is SM-ANU6G
## file is SM-ARMBP.tar
## file prefix is SM-ARMBP
## zip_type is VIR
## uid is SM-ARMBP
## file is SM-B58FT.tar
## file prefix is SM-B58FT
## zip_type is VIR
## uid is SM-B58FT
## file is SM-C1MZ2.tar
## file prefix is SM-C1MZ2
## zip_type is VIR
## uid is SM-C1MZ2
## file is SM-CHS71.tar
## file prefix is SM-CHS71
## zip_type is VIR
## uid is SM-CHS71
## file is SM-CTLFW.tar
## file prefix is SM-CTLFW
## zip_type is VIR
## uid is SM-CTLFW
## file is MSM9VZJJ_P.fastq.gz
## file prefix is MSM9VZJJ_P
## zip_type is P
## uid is MSM9VZJJ_P
## file is MSM9VZJF_P.fastq.gz
## file prefix is MSM9VZJF_P
## zip_type is P
## uid is MSM9VZJF_P
## file is MSM9VZL5.tar
## file prefix is MSM9VZL5
## zip_type is STD
## uid is MSM9VZL5
## file is MSMA26DM.tar
## file prefix is MSMA26DM
## zip_type is STD
## uid is MSMA26DM
## file is MSMA26DG.tar
## file prefix is MSMA26DG
## zip_type is STD
## uid is MSMA26DG
## file is MSMA26DI.tar
## file prefix is MSMA26DI
## zip_type is STD
## uid is MSMA26DI
## file is MSMA26DO.tar
## file prefix is MSMA26DO
## zip_type is STD
## uid is MSMA26DO
## file is MSMA26DK.tar
## file prefix is MSMA26DK
## zip_type is STD
## uid is MSMA26DK
## file is MSMAPC5H.tar
## file prefix is MSMAPC5H
## zip_type is STD
## uid is MSMAPC5H
## file is MSMAPC5L.tar
## file prefix is MSMAPC5L
## zip_type is STD
## uid is MSMAPC5L
## file is MSMB4LXW.tar
## file prefix is MSMB4LXW
## zip_type is STD
## uid is MSMB4LXW
## file is MSMB4LXS.tar
## file prefix is MSMB4LXS
## zip_type is STD
## uid is MSMB4LXS
## file is MSMB4LZK.tar
## file prefix is MSMB4LZK
## zip_type is STD
## uid is MSMB4LZK
## file is SM-AKTY2.tar
## file prefix is SM-AKTY2
## zip_type is VIR
## uid is SM-AKTY2
## file is SM-APDXU.tar
## file prefix is SM-APDXU
## zip_type is VIR
## uid is SM-APDXU
## file is SM-AVDBD.tar
## file prefix is SM-AVDBD
## zip_type is VIR
## uid is SM-AVDBD
## file is SM-AVR6V.tar
## file prefix is SM-AVR6V
## zip_type is VIR
## uid is SM-AVR6V
## file is SM-C1MYX.tar
## file prefix is SM-C1MYX
## zip_type is VIR
## uid is SM-C1MYX
## file is SM-CJNBG.tar
## file prefix is SM-CJNBG
## zip_type is VIR
## uid is SM-CJNBG
## file is SM-CTTKX.tar
## file prefix is SM-CTTKX
## zip_type is VIR
## uid is SM-CTTKX
## file is MSMA267V.tar
## file prefix is MSMA267V
## zip_type is STD
## uid is MSMA267V
## file is MSMA267X.tar
## file prefix is MSMA267X
## zip_type is STD
## uid is MSMA267X
## file is MSMA2684.tar
## file prefix is MSMA2684
## zip_type is STD
## uid is MSMA2684
## file is MSMAPC5Z.tar
## file prefix is MSMAPC5Z
## zip_type is STD
## uid is MSMAPC5Z
## file is MSMAPC66.tar
## file prefix is MSMAPC66
## zip_type is STD
## uid is MSMAPC66
## file is MSMAPC64.tar
## file prefix is MSMAPC64
## zip_type is STD
## uid is MSMAPC64
## file is MSMB4LYH.tar
## file prefix is MSMB4LYH
## zip_type is STD
## uid is MSMB4LYH
## file is SM-APDWR.tar
## file prefix is SM-APDWR
## zip_type is VIR
## uid is SM-APDWR
## file is SM-B58FX.tar
## file prefix is SM-B58FX
## zip_type is VIR
## uid is SM-B58FX
## file is SM-CRP6D.tar
## file prefix is SM-CRP6D
## zip_type is VIR
## uid is SM-CRP6D
## file is PSM6XBQM_P.fastq.gz
## file prefix is PSM6XBQM_P
## zip_type is P
## uid is PSM6XBQM_P
## file is PSM6XBQS.tar
## file prefix is PSM6XBQS
## zip_type is STD
## uid is PSM6XBQS
## file is PSM6XBQU.tar
## file prefix is PSM6XBQU
## zip_type is STD
## uid is PSM6XBQU
## file is PSM6XBQY.tar
## file prefix is PSM6XBQY
## zip_type is STD
## uid is PSM6XBQY
## file is PSM6XBQY_TR.tar
## file prefix is PSM6XBQY_TR
## zip_type is TR
## uid is PSM6XBQY_TR
## file is PSM6XBR1.tar
## file prefix is PSM6XBR1
## zip_type is STD
## uid is PSM6XBR1
## file is PSM6XBSS.tar
## file prefix is PSM6XBSS
## zip_type is STD
## uid is PSM6XBSS
## file is PSM6XBSU.tar
## file prefix is PSM6XBSU
## zip_type is STD
## uid is PSM6XBSU
## file is PSM6XBSU_TR.tar
## file prefix is PSM6XBSU_TR
## zip_type is TR
## uid is PSM6XBSU_TR
## file is PSM6XBT1.tar
## file prefix is PSM6XBT1
## zip_type is STD
## uid is PSM6XBT1
## file is PSM7J1C8.tar
## file prefix is PSM7J1C8
## zip_type is STD
## uid is PSM7J1C8
## file is PSM7J1CC.tar
## file prefix is PSM7J1CC
## zip_type is STD
## uid is PSM7J1CC
## file is PSM7J1CG.tar
## file prefix is PSM7J1CG
## zip_type is STD
## uid is PSM7J1CG
## file is PSM7J1CI.tar
## file prefix is PSM7J1CI
## zip_type is STD
## uid is PSM7J1CI
## file is PSM7J136.tar
## file prefix is PSM7J136
## zip_type is STD
## uid is PSM7J136
## file is PSM7J13E.tar
## file prefix is PSM7J13E
## zip_type is STD
## uid is PSM7J13E
## file is SM-791AV.tar
## file prefix is SM-791AV
## zip_type is VIR
## uid is SM-791AV
## file is SM-7CRWY.tar
## file prefix is SM-7CRWY
## zip_type is VIR
## uid is SM-7CRWY
## file is SM-7KPV9.tar
## file prefix is SM-7KPV9
## zip_type is VIR
## uid is SM-7KPV9
## file is SM-9P5FP.tar
## file prefix is SM-9P5FP
## zip_type is VIR
## uid is SM-9P5FP
## file is SM-9XIP8.tar
## file prefix is SM-9XIP8
## zip_type is VIR
## uid is SM-9XIP8
## file is SM-A2HYV.tar
## file prefix is SM-A2HYV
## zip_type is VIR
## uid is SM-A2HYV
## file is SM-ACTNV.tar
## file prefix is SM-ACTNV
## zip_type is VIR
## uid is SM-ACTNV
## file is PSM6XBRK_P.fastq.gz
## file prefix is PSM6XBRK_P
## zip_type is P
## uid is PSM6XBRK_P
## file is PSM6XBRK.tar
## file prefix is PSM6XBRK
## zip_type is STD
## uid is PSM6XBRK
## file is PSM6XBRK_TR.tar
## file prefix is PSM6XBRK_TR
## zip_type is TR
## uid is PSM6XBRK_TR
## file is PSM6XBS2_P.fastq.gz
## file prefix is PSM6XBS2_P
## zip_type is P
## uid is PSM6XBS2_P
## file is PSM6XBS4_P.fastq.gz
## file prefix is PSM6XBS4_P
## zip_type is P
## uid is PSM6XBS4_P
## file is PSM6XBS4.tar
## file prefix is PSM6XBS4
## zip_type is STD
## uid is PSM6XBS4
## file is PSM6XBS8.tar
## file prefix is PSM6XBS8
## zip_type is STD
## uid is PSM6XBS8
## file is PSM6XBSA.tar
## file prefix is PSM6XBSA
## zip_type is STD
## uid is PSM6XBSA
## file is PSM6XBSC_P.fastq.gz
## file prefix is PSM6XBSC_P
## zip_type is P
## uid is PSM6XBSC_P
## file is PSM6XBTF_P.fastq.gz
## file prefix is PSM6XBTF_P
## zip_type is P
## uid is PSM6XBTF_P
## file is PSM6XBTH_P.fastq.gz
## file prefix is PSM6XBTH_P
## zip_type is P
## uid is PSM6XBTH_P
## file is PSM6XBTL.tar
## file prefix is PSM6XBTL
## zip_type is STD
## uid is PSM6XBTL
## file is PSM6XBTN_P.fastq.gz
## file prefix is PSM6XBTN_P
## zip_type is P
## uid is PSM6XBTN_P
## file is PSM6XBTP.tar
## file prefix is PSM6XBTP
## zip_type is STD
## uid is PSM6XBTP
## file is PSM7J1CK_P.fastq.gz
## file prefix is PSM7J1CK_P
## zip_type is P
## uid is PSM7J1CK_P
## file is PSM7J1CS_P.fastq.gz
## file prefix is PSM7J1CS_P
## zip_type is P
## uid is PSM7J1CS_P
## file is PSM7J1CU.tar
## file prefix is PSM7J1CU
## zip_type is STD
## uid is PSM7J1CU
## file is PSM7J13U_P.fastq.gz
## file prefix is PSM7J13U_P
## zip_type is P
## uid is PSM7J13U_P
## file is PSM7J13Y_P.fastq.gz
## file prefix is PSM7J13Y_P
## zip_type is P
## uid is PSM7J13Y_P
## file is PSM7J13Y.tar
## file prefix is PSM7J13Y
## zip_type is STD
## uid is PSM7J13Y
## file is PSM7J141.tar
## file prefix is PSM7J141
## zip_type is STD
## uid is PSM7J141
## file is PSM7J143_P.fastq.gz
## file prefix is PSM7J143_P
## zip_type is P
## uid is PSM7J143_P
## file is SM-791B4.tar
## file prefix is SM-791B4
## zip_type is VIR
## uid is SM-791B4
## file is SM-7CP2S.tar
## file prefix is SM-7CP2S
## zip_type is VIR
## uid is SM-7CP2S
## file is SM-9MNBF.tar
## file prefix is SM-9MNBF
## zip_type is VIR
## uid is SM-9MNBF
## file is SM-A1ZUX.tar
## file prefix is SM-A1ZUX
## zip_type is VIR
## uid is SM-A1ZUX
## file is SM-ACTNZ.tar
## file prefix is SM-ACTNZ
## zip_type is VIR
## uid is SM-ACTNZ
## file is PSM6XBRM_P.fastq.gz
## file prefix is PSM6XBRM_P
## zip_type is P
## uid is PSM6XBRM_P
## file is PSM6XBSI.tar
## file prefix is PSM6XBSI
## zip_type is STD
## uid is PSM6XBSI
## file is PSM6XBSO.tar
## file prefix is PSM6XBSO
## zip_type is STD
## uid is PSM6XBSO
## file is PSM6XBSM.tar
## file prefix is PSM6XBSM
## zip_type is STD
## uid is PSM6XBSM
## file is PSM6XBSK.tar
## file prefix is PSM6XBSK
## zip_type is STD
## uid is PSM6XBSK
## file is PSM6XBV2.tar
## file prefix is PSM6XBV2
## zip_type is STD
## uid is PSM6XBV2
## file is PSM6XBV4.tar
## file prefix is PSM6XBV4
## zip_type is STD
## uid is PSM6XBV4
## file is PSM6XBUG.tar
## file prefix is PSM6XBUG
## zip_type is STD
## uid is PSM6XBUG
## file is PSM6XBUI.tar
## file prefix is PSM6XBUI
## zip_type is STD
## uid is PSM6XBUI
## file is PSM6XBUM.tar
## file prefix is PSM6XBUM
## zip_type is STD
## uid is PSM6XBUM
## file is PSM6XBUK.tar
## file prefix is PSM6XBUK
## zip_type is STD
## uid is PSM6XBUK
## file is PSM6XBUQ.tar
## file prefix is PSM6XBUQ
## zip_type is STD
## uid is PSM6XBUQ
## file is PSM6XBUO.tar
## file prefix is PSM6XBUO
## zip_type is STD
## uid is PSM6XBUO
## file is PSM7J18E.tar
## file prefix is PSM7J18E
## zip_type is STD
## uid is PSM7J18E
## file is PSM7J18G.tar
## file prefix is PSM7J18G
## zip_type is STD
## uid is PSM7J18G
## file is PSM7J18I.tar
## file prefix is PSM7J18I
## zip_type is STD
## uid is PSM7J18I
## file is PSM7J18K.tar
## file prefix is PSM7J18K
## zip_type is STD
## uid is PSM7J18K
## file is PSM7J18M.tar
## file prefix is PSM7J18M
## zip_type is STD
## uid is PSM7J18M
## file is PSM7J14L.tar
## file prefix is PSM7J14L
## zip_type is STD
## uid is PSM7J14L
## file is PSM7J14N.tar
## file prefix is PSM7J14N
## zip_type is STD
## uid is PSM7J14N
## file is PSM7J14P.tar
## file prefix is PSM7J14P
## zip_type is STD
## uid is PSM7J14P
## file is PSM7J14R.tar
## file prefix is PSM7J14R
## zip_type is STD
## uid is PSM7J14R
## file is PSM7J14T.tar
## file prefix is PSM7J14T
## zip_type is STD
## uid is PSM7J14T
## file is SM-7DN38.tar
## file prefix is SM-7DN38
## zip_type is VIR
## uid is SM-7DN38
## file is SM-7EWSN.tar
## file prefix is SM-7EWSN
## zip_type is VIR
## uid is SM-7EWSN
## file is SM-7GYIW.tar
## file prefix is SM-7GYIW
## zip_type is VIR
## uid is SM-7GYIW
## file is SM-7KMRD.tar
## file prefix is SM-7KMRD
## zip_type is VIR
## uid is SM-7KMRD
## file is SM-7M8RV.tar
## file prefix is SM-7M8RV
## zip_type is VIR
## uid is SM-7M8RV
## file is SM-7R3AF.tar
## file prefix is SM-7R3AF
## zip_type is VIR
## uid is SM-7R3AF
## file is SM-9IMBP.tar
## file prefix is SM-9IMBP
## zip_type is VIR
## uid is SM-9IMBP
## file is SM-9NBDT.tar
## file prefix is SM-9NBDT
## zip_type is VIR
## uid is SM-9NBDT
## file is SM-9QMQ2.tar
## file prefix is SM-9QMQ2
## zip_type is VIR
## uid is SM-9QMQ2
## file is SM-9VWBA.tar
## file prefix is SM-9VWBA
## zip_type is VIR
## uid is SM-9VWBA
## file is SM-9XG81.tar
## file prefix is SM-9XG81
## zip_type is VIR
## uid is SM-9XG81
## file is SM-9ZA1U.tar
## file prefix is SM-9ZA1U
## zip_type is VIR
## uid is SM-9ZA1U
## file is SM-A2J34.tar
## file prefix is SM-A2J34
## zip_type is VIR
## uid is SM-A2J34
## file is SM-A4O3F.tar
## file prefix is SM-A4O3F
## zip_type is VIR
## uid is SM-A4O3F
## file is SM-A9F6H.tar
## file prefix is SM-A9F6H
## zip_type is VIR
## uid is SM-A9F6H
## file is SM-ACTMJ.tar
## file prefix is SM-ACTMJ
## zip_type is VIR
## uid is SM-ACTMJ
## file is SM-AFCJZ.tar
## file prefix is SM-AFCJZ
## zip_type is VIR
## uid is SM-AFCJZ
## file is SM-AGUM9.tar
## file prefix is SM-AGUM9
## zip_type is VIR
## uid is SM-AGUM9
## file is PSM6XBSE_P.fastq.gz
## file prefix is PSM6XBSE_P
## zip_type is P
## uid is PSM6XBSE_P
## file is PSM6XBSE.tar
## file prefix is PSM6XBSE
## zip_type is STD
## uid is PSM6XBSE
## file is PSM6XBVI.tar
## file prefix is PSM6XBVI
## zip_type is STD
## uid is PSM6XBVI
## file is PSM6XBVK.tar
## file prefix is PSM6XBVK
## zip_type is STD
## uid is PSM6XBVK
## file is PSM6XBVM.tar
## file prefix is PSM6XBVM
## zip_type is STD
## uid is PSM6XBVM
## file is PSM6XBVO.tar
## file prefix is PSM6XBVO
## zip_type is STD
## uid is PSM6XBVO
## file is PSM6XBVQ.tar
## file prefix is PSM6XBVQ
## zip_type is STD
## uid is PSM6XBVQ
## file is PSM6XBVS.tar
## file prefix is PSM6XBVS
## zip_type is STD
## uid is PSM6XBVS
## file is PSM7J1B7.tar
## file prefix is PSM7J1B7
## zip_type is STD
## uid is PSM7J1B7
## file is PSM7J1B9.tar
## file prefix is PSM7J1B9
## zip_type is STD
## uid is PSM7J1B9
## file is PSM7J1BB.tar
## file prefix is PSM7J1BB
## zip_type is STD
## uid is PSM7J1BB
## file is PSM7J1BD.tar
## file prefix is PSM7J1BD
## zip_type is STD
## uid is PSM7J1BD
## file is PSM7J1BF.tar
## file prefix is PSM7J1BF
## zip_type is STD
## uid is PSM7J1BF
## file is PSM7J1BH.tar
## file prefix is PSM7J1BH
## zip_type is STD
## uid is PSM7J1BH
## file is PSM7J18Q.tar
## file prefix is PSM7J18Q
## zip_type is STD
## uid is PSM7J18Q
## file is SM-7EWRZ.tar
## file prefix is SM-7EWRZ
## zip_type is VIR
## uid is SM-7EWRZ
## file is SM-7H4H6.tar
## file prefix is SM-7H4H6
## zip_type is VIR
## uid is SM-7H4H6
## file is SM-9NK5X.tar
## file prefix is SM-9NK5X
## zip_type is VIR
## uid is SM-9NK5X
## file is SM-9SIIV.tar
## file prefix is SM-9SIIV
## zip_type is VIR
## uid is SM-9SIIV
## file is SM-A77W9.tar
## file prefix is SM-A77W9
## zip_type is VIR
## uid is SM-A77W9
## file is PSM6XBSG_P.fastq.gz
## file prefix is PSM6XBSG_P
## zip_type is P
## uid is PSM6XBSG_P
## file is PSM6XBT3.tar
## file prefix is PSM6XBT3
## zip_type is STD
## uid is PSM6XBT3
## file is PSM6XBT5.tar
## file prefix is PSM6XBT5
## zip_type is STD
## uid is PSM6XBT5
## file is PSM6XBT7.tar
## file prefix is PSM6XBT7
## zip_type is STD
## uid is PSM6XBT7
## file is PSM6XBT9.tar
## file prefix is PSM6XBT9
## zip_type is STD
## uid is PSM6XBT9
## file is PSM6XBTB.tar
## file prefix is PSM6XBTB
## zip_type is STD
## uid is PSM6XBTB
## file is PSM6XBTD.tar
## file prefix is PSM6XBTD
## zip_type is STD
## uid is PSM6XBTD
## file is PSM7J1AM.tar
## file prefix is PSM7J1AM
## zip_type is STD
## uid is PSM7J1AM
## file is PSM7J1AS.tar
## file prefix is PSM7J1AS
## zip_type is STD
## uid is PSM7J1AS
## file is PSM7J1AU.tar
## file prefix is PSM7J1AU
## zip_type is STD
## uid is PSM7J1AU
## file is PSM7J1AO.tar
## file prefix is PSM7J1AO
## zip_type is STD
## uid is PSM7J1AO
## file is PSM7J1AQ.tar
## file prefix is PSM7J1AQ
## zip_type is STD
## uid is PSM7J1AQ
## file is PSM7J1AW.tar
## file prefix is PSM7J1AW
## zip_type is STD
## uid is PSM7J1AW
## file is PSM7J193.tar
## file prefix is PSM7J193
## zip_type is STD
## uid is PSM7J193
## file is PSM7J127.tar
## file prefix is PSM7J127
## zip_type is STD
## uid is PSM7J127
## file is PSM7J129.tar
## file prefix is PSM7J129
## zip_type is STD
## uid is PSM7J129
## file is PSM7J12B.tar
## file prefix is PSM7J12B
## zip_type is STD
## uid is PSM7J12B
## file is PSM7J12D.tar
## file prefix is PSM7J12D
## zip_type is STD
## uid is PSM7J12D
## file is PSM7J12F.tar
## file prefix is PSM7J12F
## zip_type is STD
## uid is PSM7J12F
## file is PSM7J169.tar
## file prefix is PSM7J169
## zip_type is STD
## uid is PSM7J169
## file is PSM7J16F.tar
## file prefix is PSM7J16F
## zip_type is STD
## uid is PSM7J16F
## file is PSM7J16B.tar
## file prefix is PSM7J16B
## zip_type is STD
## uid is PSM7J16B
## file is PSM7J16H.tar
## file prefix is PSM7J16H
## zip_type is STD
## uid is PSM7J16H
## file is SM-7H4GR.tar
## file prefix is SM-7H4GR
## zip_type is VIR
## uid is SM-7H4GR
## file is SM-7M8RR.tar
## file prefix is SM-7M8RR
## zip_type is VIR
## uid is SM-7M8RR
## file is SM-9YHLA.tar
## file prefix is SM-9YHLA
## zip_type is VIR
## uid is SM-9YHLA
## file is SM-A9U2H.tar
## file prefix is SM-A9U2H
## zip_type is VIR
## uid is SM-A9U2H
## file is PSM6XBTR.tar
## file prefix is PSM6XBTR
## zip_type is STD
## uid is PSM6XBTR
## file is PSM6XBTT.tar
## file prefix is PSM6XBTT
## zip_type is STD
## uid is PSM6XBTT
## file is PSM6XBTX.tar
## file prefix is PSM6XBTX
## zip_type is STD
## uid is PSM6XBTX
## file is PSM6XBTZ_P.fastq.gz
## file prefix is PSM6XBTZ_P
## zip_type is P
## uid is PSM6XBTZ_P
## file is PSM6XBU2.tar
## file prefix is PSM6XBU2
## zip_type is STD
## uid is PSM6XBU2
## file is PSM7J1DF.tar
## file prefix is PSM7J1DF
## zip_type is STD
## uid is PSM7J1DF
## file is PSM7J1DL.tar
## file prefix is PSM7J1DL
## zip_type is STD
## uid is PSM7J1DL
## file is PSM7J14X.tar
## file prefix is PSM7J14X
## zip_type is STD
## uid is PSM7J14X
## file is PSM7J154.tar
## file prefix is PSM7J154
## zip_type is STD
## uid is PSM7J154
## file is PSM7J156.tar
## file prefix is PSM7J156
## zip_type is STD
## uid is PSM7J156
## file is SM-A8MJR.tar
## file prefix is SM-A8MJR
## zip_type is VIR
## uid is SM-A8MJR
## file is SM-AF6NI.tar
## file prefix is SM-AF6NI
## zip_type is VIR
## uid is SM-AF6NI
## file is PSM6XBVY_P.fastq.gz
## file prefix is PSM6XBVY_P
## zip_type is P
## uid is PSM6XBVY_P
## file is PSM7J199.tar
## file prefix is PSM7J199
## zip_type is STD
## uid is PSM7J199
## file is PSM7J19B.tar
## file prefix is PSM7J19B
## zip_type is STD
## uid is PSM7J19B
## file is PSM7J19F.tar
## file prefix is PSM7J19F
## zip_type is STD
## uid is PSM7J19F
## file is PSM7J19H.tar
## file prefix is PSM7J19H
## zip_type is STD
## uid is PSM7J19H
## file is PSM7J19J.tar
## file prefix is PSM7J19J
## zip_type is STD
## uid is PSM7J19J
## file is PSM7J17L.tar
## file prefix is PSM7J17L
## zip_type is STD
## uid is PSM7J17L
## file is PSM7J17T.tar
## file prefix is PSM7J17T
## zip_type is STD
## uid is PSM7J17T
## file is PSM7J158.tar
## file prefix is PSM7J158
## zip_type is STD
## uid is PSM7J158
## file is PSM7J15A.tar
## file prefix is PSM7J15A
## zip_type is STD
## uid is PSM7J15A
## file is PSM7J15G.tar
## file prefix is PSM7J15G
## zip_type is STD
## uid is PSM7J15G
## file is PSM7J15I.tar
## file prefix is PSM7J15I
## zip_type is STD
## uid is PSM7J15I
## file is PSMA265N.tar
## file prefix is PSMA265N
## zip_type is STD
## uid is PSMA265N
## file is PSMA265T.tar
## file prefix is PSMA265T
## zip_type is STD
## uid is PSMA265T
## file is SM-7MX7Z.tar
## file prefix is SM-7MX7Z
## zip_type is VIR
## uid is SM-7MX7Z
## file is SM-7RFYI.tar
## file prefix is SM-7RFYI
## zip_type is VIR
## uid is SM-7RFYI
## file is SM-9KOP2.tar
## file prefix is SM-9KOP2
## zip_type is VIR
## uid is SM-9KOP2
## file is SM-9QMNE.tar
## file prefix is SM-9QMNE
## zip_type is VIR
## uid is SM-9QMNE
## file is SM-9VWAV.tar
## file prefix is SM-9VWAV
## zip_type is VIR
## uid is SM-9VWAV
## file is SM-A9F5U.tar
## file prefix is SM-A9F5U
## zip_type is VIR
## uid is SM-A9F5U
## file is SM-ARMCD.tar
## file prefix is SM-ARMCD
## zip_type is VIR
## uid is SM-ARMCD
## file is PSM6XBW1_P.fastq.gz
## file prefix is PSM6XBW1_P
## zip_type is P
## uid is PSM6XBW1_P
## file is PSM7J19N.tar
## file prefix is PSM7J19N
## zip_type is STD
## uid is PSM7J19N
## file is PSM7J19P.tar
## file prefix is PSM7J19P
## zip_type is STD
## uid is PSM7J19P
## file is PSM7J19R.tar
## file prefix is PSM7J19R
## zip_type is STD
## uid is PSM7J19R
## file is PSM7J19T.tar
## file prefix is PSM7J19T
## zip_type is STD
## uid is PSM7J19T
## file is PSM7J16U.tar
## file prefix is PSM7J16U
## zip_type is STD
## uid is PSM7J16U
## file is PSM7J16W.tar
## file prefix is PSM7J16W
## zip_type is STD
## uid is PSM7J16W
## file is PSM7J16Y.tar
## file prefix is PSM7J16Y
## zip_type is STD
## uid is PSM7J16Y
## file is PSM7J171.tar
## file prefix is PSM7J171
## zip_type is STD
## uid is PSM7J171
## file is PSM7J173.tar
## file prefix is PSM7J173
## zip_type is STD
## uid is PSM7J173
## file is SM-7T2LJ.tar
## file prefix is SM-7T2LJ
## zip_type is VIR
## uid is SM-7T2LJ
## file is SM-9HYXU.tar
## file prefix is SM-9HYXU
## zip_type is VIR
## uid is SM-9HYXU
## file is SM-ACTOG.tar
## file prefix is SM-ACTOG
## zip_type is VIR
## uid is SM-ACTOG
## file is SM-AETRB.tar
## file prefix is SM-AETRB
## zip_type is VIR
## uid is SM-AETRB
## file is PSM6XBW3.tar
## file prefix is PSM6XBW3
## zip_type is STD
## uid is PSM6XBW3
## file is PSM7J19X_P.fastq.gz
## file prefix is PSM7J19X_P
## zip_type is P
## uid is PSM7J19X_P
## file is PSM7J19Z.tar
## file prefix is PSM7J19Z
## zip_type is STD
## uid is PSM7J19Z
## file is PSM7J1A2.tar
## file prefix is PSM7J1A2
## zip_type is STD
## uid is PSM7J1A2
## file is PSM7J1A4.tar
## file prefix is PSM7J1A4
## zip_type is STD
## uid is PSM7J1A4
## file is PSM7J1A6.tar
## file prefix is PSM7J1A6
## zip_type is STD
## uid is PSM7J1A6
## file is PSM7J1A8.tar
## file prefix is PSM7J1A8
## zip_type is STD
## uid is PSM7J1A8
## file is PSM7J17V.tar
## file prefix is PSM7J17V
## zip_type is STD
## uid is PSM7J17V
## file is PSM7J17X.tar
## file prefix is PSM7J17X
## zip_type is STD
## uid is PSM7J17X
## file is PSM7J17Z.tar
## file prefix is PSM7J17Z
## zip_type is STD
## uid is PSM7J17Z
## file is PSM7J182.tar
## file prefix is PSM7J182
## zip_type is STD
## uid is PSM7J182
## file is PSM7J184.tar
## file prefix is PSM7J184
## zip_type is STD
## uid is PSM7J184
## file is PSM7J186.tar
## file prefix is PSM7J186
## zip_type is STD
## uid is PSM7J186
## file is PSM7J15K.tar
## file prefix is PSM7J15K
## zip_type is STD
## uid is PSM7J15K
## file is PSM7J15M.tar
## file prefix is PSM7J15M
## zip_type is STD
## uid is PSM7J15M
## file is PSM7J15O.tar
## file prefix is PSM7J15O
## zip_type is STD
## uid is PSM7J15O
## file is PSM7J15Q.tar
## file prefix is PSM7J15Q
## zip_type is STD
## uid is PSM7J15Q
## file is PSM7J15S.tar
## file prefix is PSM7J15S
## zip_type is STD
## uid is PSM7J15S
## file is PSM7J15U.tar
## file prefix is PSM7J15U
## zip_type is STD
## uid is PSM7J15U
## file is PSMA265D.tar
## file prefix is PSMA265D
## zip_type is STD
## uid is PSMA265D
## file is PSMA265F.tar
## file prefix is PSMA265F
## zip_type is STD
## uid is PSMA265F
## file is PSMA265L.tar
## file prefix is PSMA265L
## zip_type is STD
## uid is PSMA265L
## file is PSMA265J.tar
## file prefix is PSMA265J
## zip_type is STD
## uid is PSMA265J
## file is PSMA265J_TR.tar
## file prefix is PSMA265J_TR
## zip_type is TR
## uid is PSMA265J_TR
## file is PSMA265H.tar
## file prefix is PSMA265H
## zip_type is STD
## uid is PSMA265H
## file is SM-9J5IS.tar
## file prefix is SM-9J5IS
## zip_type is VIR
## uid is SM-9J5IS
## file is SM-9NBDP.tar
## file prefix is SM-9NBDP
## zip_type is VIR
## uid is SM-9NBDP
## file is SM-9SNKC.tar
## file prefix is SM-9SNKC
## zip_type is VIR
## uid is SM-9SNKC
## file is SM-9VW8U.tar
## file prefix is SM-9VW8U
## zip_type is VIR
## uid is SM-9VW8U
## file is SM-9XG7U.tar
## file prefix is SM-9XG7U
## zip_type is VIR
## uid is SM-9XG7U
## file is SM-9ZA5Q.tar
## file prefix is SM-9ZA5Q
## zip_type is VIR
## uid is SM-9ZA5Q
## file is SM-A1ZWL.tar
## file prefix is SM-A1ZWL
## zip_type is VIR
## uid is SM-A1ZWL
## file is SM-A4O3L.tar
## file prefix is SM-A4O3L
## zip_type is VIR
## uid is SM-A4O3L
## file is SM-ACTMR.tar
## file prefix is SM-ACTMR
## zip_type is VIR
## uid is SM-ACTMR
## file is SM-AF6MB.tar
## file prefix is SM-AF6MB
## zip_type is VIR
## uid is SM-AF6MB
## file is SM-AGUM1.tar
## file prefix is SM-AGUM1
## zip_type is VIR
## uid is SM-AGUM1
## file is SM-AMR1R.tar
## file prefix is SM-AMR1R
## zip_type is VIR
## uid is SM-AMR1R
## file is SM-APDXY.tar
## file prefix is SM-APDXY
## zip_type is VIR
## uid is SM-APDXY
## file is SM-ARGG8.tar
## file prefix is SM-ARGG8
## zip_type is VIR
## uid is SM-ARGG8
## file is SM-ASN2T.tar
## file prefix is SM-ASN2T
## zip_type is VIR
## uid is SM-ASN2T
## file is PSM7J1BJ.tar
## file prefix is PSM7J1BJ
## zip_type is STD
## uid is PSM7J1BJ
## file is PSM7J1BL.tar
## file prefix is PSM7J1BL
## zip_type is STD
## uid is PSM7J1BL
## file is PSM7J1BN_P.fastq.gz
## file prefix is PSM7J1BN_P
## zip_type is P
## uid is PSM7J1BN_P
## file is PSM7J1BP.tar
## file prefix is PSM7J1BP
## zip_type is STD
## uid is PSM7J1BP
## file is PSM7J1BR.tar
## file prefix is PSM7J1BR
## zip_type is STD
## uid is PSM7J1BR
## file is PSM7J12J.tar
## file prefix is PSM7J12J
## zip_type is STD
## uid is PSM7J12J
## file is PSM7J12R.tar
## file prefix is PSM7J12R
## zip_type is STD
## uid is PSM7J12R
## file is PSMA263M.tar
## file prefix is PSMA263M
## zip_type is STD
## uid is PSMA263M
## file is PSMA263S.tar
## file prefix is PSMA263S
## zip_type is STD
## uid is PSMA263S
## file is PSMA263U.tar
## file prefix is PSMA263U
## zip_type is STD
## uid is PSMA263U
## file is PSMA263W.tar
## file prefix is PSMA263W
## zip_type is STD
## uid is PSMA263W
## file is SM-9U1UP.tar
## file prefix is SM-9U1UP
## zip_type is VIR
## uid is SM-9U1UP
## file is SM-9ZA1X.tar
## file prefix is SM-9ZA1X
## zip_type is VIR
## uid is SM-9ZA1X
## file is SM-A9J96.tar
## file prefix is SM-A9J96
## zip_type is VIR
## uid is SM-A9J96
## file is SM-APR9N.tar
## file prefix is SM-APR9N
## zip_type is VIR
## uid is SM-APR9N
## file is PSM7J1BV.tar
## file prefix is PSM7J1BV
## zip_type is STD
## uid is PSM7J1BV
## file is PSM7J1BX.tar
## file prefix is PSM7J1BX
## zip_type is STD
## uid is PSM7J1BX
## file is PSM7J1C2_P.fastq.gz
## file prefix is PSM7J1C2_P
## zip_type is P
## uid is PSM7J1C2_P
## file is PSM7J1C4.tar
## file prefix is PSM7J1C4
## zip_type is STD
## uid is PSM7J1C4
## file is PSM7J12V.tar
## file prefix is PSM7J12V
## zip_type is STD
## uid is PSM7J12V
## file is PSM7J12Z.tar
## file prefix is PSM7J12Z
## zip_type is STD
## uid is PSM7J12Z
## file is PSMA264K.tar
## file prefix is PSMA264K
## zip_type is STD
## uid is PSMA264K
## file is PSM7J1B3_P.fastq.gz
## file prefix is PSM7J1B3_P
## zip_type is P
## uid is PSM7J1B3_P
## file is PSM7J177.tar
## file prefix is PSM7J177
## zip_type is STD
## uid is PSM7J177
## file is PSM7J179.tar
## file prefix is PSM7J179
## zip_type is STD
## uid is PSM7J179
## file is PSM7J17B.tar
## file prefix is PSM7J17B
## zip_type is STD
## uid is PSM7J17B
## file is PSM7J17D.tar
## file prefix is PSM7J17D
## zip_type is STD
## uid is PSM7J17D
## file is PSM7J17F.tar
## file prefix is PSM7J17F
## zip_type is STD
## uid is PSM7J17F
## file is PSM7J15W.tar
## file prefix is PSM7J15W
## zip_type is STD
## uid is PSM7J15W
## file is PSM7J161.tar
## file prefix is PSM7J161
## zip_type is STD
## uid is PSM7J161
## file is PSM7J163.tar
## file prefix is PSM7J163
## zip_type is STD
## uid is PSM7J163
## file is PSMA2668.tar
## file prefix is PSMA2668
## zip_type is STD
## uid is PSMA2668
## file is PSMA266C.tar
## file prefix is PSMA266C
## zip_type is STD
## uid is PSMA266C
## file is SM-9OS66.tar
## file prefix is SM-9OS66
## zip_type is VIR
## uid is SM-9OS66
## file is SM-9RRYZ.tar
## file prefix is SM-9RRYZ
## zip_type is VIR
## uid is SM-9RRYZ
## file is SM-9Y4A1.tar
## file prefix is SM-9Y4A1
## zip_type is VIR
## uid is SM-9Y4A1
## file is SM-A9F6B.tar
## file prefix is SM-A9F6B
## zip_type is VIR
## uid is SM-A9F6B
## file is SM-AFCKA.tar
## file prefix is SM-AFCKA
## zip_type is VIR
## uid is SM-AFCKA
## file is SM-AGUMX.tar
## file prefix is SM-AGUMX
## zip_type is VIR
## uid is SM-AGUMX
## file is SM-APR9F.tar
## file prefix is SM-APR9F
## zip_type is VIR
## uid is SM-APR9F
## file is PSM7J13I.tar
## file prefix is PSM7J13I
## zip_type is STD
## uid is PSM7J13I
## file is PSM7J13K.tar
## file prefix is PSM7J13K
## zip_type is STD
## uid is PSM7J13K
## file is PSM7J13Q.tar
## file prefix is PSM7J13Q
## zip_type is STD
## uid is PSM7J13Q
## file is PSM7J13M.tar
## file prefix is PSM7J13M
## zip_type is STD
## uid is PSM7J13M
## file is PSMA2651_P.fastq.gz
## file prefix is PSMA2651_P
## zip_type is P
## uid is PSMA2651_P
## file is PSMA2653.tar
## file prefix is PSMA2653
## zip_type is STD
## uid is PSMA2653
## file is PSMA2659.tar
## file prefix is PSMA2659
## zip_type is STD
## uid is PSMA2659
## file is PSMA265B.tar
## file prefix is PSMA265B
## zip_type is STD
## uid is PSMA265B
## file is PSMA267D.tar
## file prefix is PSMA267D
## zip_type is STD
## uid is PSMA267D
## file is PSMA267F.tar
## file prefix is PSMA267F
## zip_type is STD
## uid is PSMA267F
## file is PSMA267H.tar
## file prefix is PSMA267H
## zip_type is STD
## uid is PSMA267H
## file is SM-AADLE.tar
## file prefix is SM-AADLE
## zip_type is VIR
## uid is SM-AADLE
## file is SM-ARMD7.tar
## file prefix is SM-ARMD7
## zip_type is VIR
## uid is SM-ARMD7
## file is PSMA264O.tar
## file prefix is PSMA264O
## zip_type is STD
## uid is PSMA264O
## file is PSMA264Q.tar
## file prefix is PSMA264Q
## zip_type is STD
## uid is PSMA264Q
## file is PSMA264S.tar
## file prefix is PSMA264S
## zip_type is STD
## uid is PSMA264S
## file is PSMA264U.tar
## file prefix is PSMA264U
## zip_type is STD
## uid is PSMA264U
## file is PSMA264W.tar
## file prefix is PSMA264W
## zip_type is STD
## uid is PSMA264W
## file is PSMA267J.tar
## file prefix is PSMA267J
## zip_type is STD
## uid is PSMA267J
## file is PSMA267P.tar
## file prefix is PSMA267P
## zip_type is STD
## uid is PSMA267P
## file is PSMA267R.tar
## file prefix is PSMA267R
## zip_type is STD
## uid is PSMA267R
## file is PSMB4MBK.tar
## file prefix is PSMB4MBK
## zip_type is STD
## uid is PSMB4MBK
## file is PSMB4MBI.tar
## file prefix is PSMB4MBI
## zip_type is STD
## uid is PSMB4MBI
## file is PSMB4MC7.tar
## file prefix is PSMB4MC7
## zip_type is STD
## uid is PSMB4MC7
## file is SM-AKTYQ.tar
## file prefix is SM-AKTYQ
## zip_type is VIR
## uid is SM-AKTYQ
## file is SM-AMR2R.tar
## file prefix is SM-AMR2R
## zip_type is VIR
## uid is SM-AMR2R
## file is SM-ATLL3.tar
## file prefix is SM-ATLL3
## zip_type is VIR
## uid is SM-ATLL3
## file is SM-B1JO6.tar
## file prefix is SM-B1JO6
## zip_type is VIR
## uid is SM-B1JO6
## file is SM-B3Z96.tar
## file prefix is SM-B3Z96
## zip_type is VIR
## uid is SM-B3Z96
## file is SM-CHS75.tar
## file prefix is SM-CHS75
## zip_type is VIR
## uid is SM-CHS75
## file is SM-CL4HM.tar
## file prefix is SM-CL4HM
## zip_type is VIR
## uid is SM-CL4HM
## file is PSM7J4EF.tar
## file prefix is PSM7J4EF
## zip_type is STD
## uid is PSM7J4EF
## file is PSMA266I.tar
## file prefix is PSMA266I
## zip_type is STD
## uid is PSMA266I
## file is PSMA266M.tar
## file prefix is PSMA266M
## zip_type is STD
## uid is PSMA266M
## file is PSMA266O.tar
## file prefix is PSMA266O
## zip_type is STD
## uid is PSMA266O
## file is PSMA266Q.tar
## file prefix is PSMA266Q
## zip_type is STD
## uid is PSMA266Q
## file is PSMA269G.tar
## file prefix is PSMA269G
## zip_type is STD
## uid is PSMA269G
## file is PSMA269O.tar
## file prefix is PSMA269O
## zip_type is STD
## uid is PSMA269O
## file is PSMB4MBS.tar
## file prefix is PSMB4MBS
## zip_type is STD
## uid is PSMB4MBS
## file is SM-AMR3O.tar
## file prefix is SM-AMR3O
## zip_type is VIR
## uid is SM-AMR3O
## file is SM-ARMCP.tar
## file prefix is SM-ARMCP
## zip_type is VIR
## uid is SM-ARMCP
## file is SM-AVADS.tar
## file prefix is SM-AVADS
## zip_type is VIR
## uid is SM-AVADS
## file is SM-AYWP4.tar
## file prefix is SM-AYWP4
## zip_type is VIR
## uid is SM-AYWP4
## file is SM-BZNL3.tar
## file prefix is SM-BZNL3
## zip_type is VIR
## uid is SM-BZNL3
## file is SM-CL4HI.tar
## file prefix is SM-CL4HI
## zip_type is VIR
## uid is SM-CL4HI
## file is PSMA265X.tar
## file prefix is PSMA265X
## zip_type is STD
## uid is PSMA265X
## file is PSMA266U.tar
## file prefix is PSMA266U
## zip_type is STD
## uid is PSMA266U
## file is PSMA266Y.tar
## file prefix is PSMA266Y
## zip_type is STD
## uid is PSMA266Y
## file is PSMA2671.tar
## file prefix is PSMA2671
## zip_type is STD
## uid is PSMA2671
## file is PSMA2675.tar
## file prefix is PSMA2675
## zip_type is STD
## uid is PSMA2675
## file is PSMA269S.tar
## file prefix is PSMA269S
## zip_type is STD
## uid is PSMA269S
## file is PSMA269W.tar
## file prefix is PSMA269W
## zip_type is STD
## uid is PSMA269W
## file is PSMA26A1.tar
## file prefix is PSMA26A1
## zip_type is STD
## uid is PSMA26A1
## file is PSMA26A3.tar
## file prefix is PSMA26A3
## zip_type is STD
## uid is PSMA26A3
## file is PSMB4MC1.tar
## file prefix is PSMB4MC1
## zip_type is STD
## uid is PSMB4MC1
## file is PSMB4MC3.tar
## file prefix is PSMB4MC3
## zip_type is STD
## uid is PSMB4MC3
## file is PSMB4MC5.tar
## file prefix is PSMB4MC5
## zip_type is STD
## uid is PSMB4MC5
## file is SM-AK8KU.tar
## file prefix is SM-AK8KU
## zip_type is VIR
## uid is SM-AK8KU
## file is SM-AMR1K.tar
## file prefix is SM-AMR1K
## zip_type is VIR
## uid is SM-AMR1K
## file is SM-AWQK6.tar
## file prefix is SM-AWQK6
## zip_type is VIR
## uid is SM-AWQK6
## file is SM-AYWOV.tar
## file prefix is SM-AYWOV
## zip_type is VIR
## uid is SM-AYWOV
## file is SM-B4LHR.tar
## file prefix is SM-B4LHR
## zip_type is VIR
## uid is SM-B4LHR
## file is SM-CEFO8.tar
## file prefix is SM-CEFO8
## zip_type is VIR
## uid is SM-CEFO8
## file is SM-CRT7Q.tar
## file prefix is SM-CRT7Q
## zip_type is VIR
## uid is SM-CRT7Q
```




```bash
while IFS="," read -r file prefix zip_type uid
do
echo "on file $file with type $zip_type .... "
  case $zip_type in
    "VIR")
      echo "virus"
      ;;
    "P"|"TR"|"STD")
      echo "mtg"
      ;;
  esac
done < <(tail -n +2 csvs_and_other_metadata/file_zip_uid.csv)
```

```
## on file CSM5FZ3N_P.fastq.gz with type P .... 
## mtg
## on file CSM5FZ3R_P.fastq.gz with type P .... 
## mtg
## on file CSM5YRY7_P.fastq.gz with type P .... 
## mtg
## on file CSM5FZ3V_P.fastq.gz with type P .... 
## mtg
## on file CSM5FZ4C_P.fastq.gz with type P .... 
## mtg
## on file CSM5MCVD_P.fastq.gz with type P .... 
## mtg
## on file CSM5MCVF_P.fastq.gz with type P .... 
## mtg
## on file CSM5MCVV_P.fastq.gz with type P .... 
## mtg
## on file CSM5MCWI_P.fastq.gz with type P .... 
## mtg
## on file CSM5MCXD.tar with type STD .... 
## mtg
## on file CSM5MCYS.tar with type STD .... 
## mtg
## on file CSM67U9J.tar with type STD .... 
## mtg
## on file CSM67UA2.tar with type STD .... 
## mtg
## on file CSM67UGC.tar with type STD .... 
## mtg
## on file CSM79HG5.tar with type STD .... 
## mtg
## on file CSM79HGP.tar with type STD .... 
## mtg
## on file SM-6UG79.tar with type VIR .... 
## virus
## on file SM-71WY3.tar with type VIR .... 
## virus
## on file SM-7CRX4.tar with type VIR .... 
## virus
## on file CSM5FZ3T_P.fastq.gz with type P .... 
## mtg
## on file CSM5FZ3X_P.fastq.gz with type P .... 
## mtg
## on file CSM5FZ3Z_P.fastq.gz with type P .... 
## mtg
## on file CSM5FZ42_P.fastq.gz with type P .... 
## mtg
## on file CSM5FZ44_P.fastq.gz with type P .... 
## mtg
## on file CSM5FZ46_P.fastq.gz with type P .... 
## mtg
## on file CSM5MCVJ_P.fastq.gz with type P .... 
## mtg
## on file CSM5MCVL.tar with type STD .... 
## mtg
## on file CSM5MCVN.tar with type STD .... 
## mtg
## on file CSM67UBF.tar with type STD .... 
## mtg
## on file CSM67UBH.tar with type STD .... 
## mtg
## on file CSM67UBN.tar with type STD .... 
## mtg
## on file CSM67UBR.tar with type STD .... 
## mtg
## on file CSM79HJW.tar with type STD .... 
## mtg
## on file CSM79HJY.tar with type STD .... 
## mtg
## on file SM-6X9WV.tar with type VIR .... 
## virus
## on file SM-76CAU.tar with type VIR .... 
## virus
## on file SM-7CP3H.tar with type VIR .... 
## virus
## on file SM-7EWTT.tar with type VIR .... 
## virus
## on file CSM5FZ4E_P.fastq.gz with type P .... 
## mtg
## on file CSM5FZ4G_P.fastq.gz with type P .... 
## mtg
## on file CSM5FZ4K_P.fastq.gz with type P .... 
## mtg
## on file CSM5FZ4M.tar with type STD .... 
## mtg
## on file CSM5MCWM_P.fastq.gz with type P .... 
## mtg
## on file CSM5MCWQ.tar with type STD .... 
## mtg
## on file CSM67UBX.tar with type STD .... 
## mtg
## on file CSM67UBZ.tar with type STD .... 
## mtg
## on file CSM67UC6.tar with type STD .... 
## mtg
## on file CSM79HLM.tar with type STD .... 
## mtg
## on file SM-6OSOR.tar with type VIR .... 
## virus
## on file SM-6Y2V3.tar with type VIR .... 
## virus
## on file SM-77VQ4.tar with type VIR .... 
## virus
## on file SM-7EWUL.tar with type VIR .... 
## virus
## on file CSM5FZ4A_P.fastq.gz with type P .... 
## mtg
## on file CSM5MCU8_P.fastq.gz with type P .... 
## mtg
## on file CSM5MCUA_P.fastq.gz with type P .... 
## mtg
## on file CSM5MCUC_P.fastq.gz with type P .... 
## mtg
## on file CSM5MCUE_P.fastq.gz with type P .... 
## mtg
## on file CSM5MCXH.tar with type STD .... 
## mtg
## on file CSM5MCXJ.tar with type STD .... 
## mtg
## on file CSM5MCXL.tar with type STD .... 
## mtg
## on file CSM5MCXN.tar with type STD .... 
## mtg
## on file CSM5MCXP.tar with type STD .... 
## mtg
## on file CSM5MCXR.tar with type STD .... 
## mtg
## on file CSM67UDF.tar with type STD .... 
## mtg
## on file CSM67UDJ.tar with type STD .... 
## mtg
## on file CSM67UDN.tar with type STD .... 
## mtg
## on file CSM67UDR_TR.tar with type TR .... 
## mtg
## on file CSM67UDR.tar with type STD .... 
## mtg
## on file CSM67UDY.tar with type STD .... 
## mtg
## on file CSM79HLA_TR.tar with type TR .... 
## mtg
## on file CSM79HLA.tar with type STD .... 
## mtg
## on file CSM79HLG.tar with type STD .... 
## mtg
## on file CSM79HLE.tar with type STD .... 
## mtg
## on file CSM79HLC.tar with type STD .... 
## mtg
## on file CSM79HLI.tar with type STD .... 
## mtg
## on file CSM79HLK.tar with type STD .... 
## mtg
## on file SM-6WJN6.tar with type VIR .... 
## virus
## on file SM-6X9X4.tar with type VIR .... 
## virus
## on file SM-6YAZO.tar with type VIR .... 
## virus
## on file SM-6ZKSM.tar with type VIR .... 
## virus
## on file SM-71WXM.tar with type VIR .... 
## virus
## on file SM-73JY4.tar with type VIR .... 
## virus
## on file SM-76EOJ.tar with type VIR .... 
## virus
## on file SM-791BX.tar with type VIR .... 
## virus
## on file SM-7BF2J.tar with type VIR .... 
## virus
## on file SM-7CRJ8.tar with type VIR .... 
## virus
## on file SM-7EWUH.tar with type VIR .... 
## virus
## on file SM-7GYJO.tar with type VIR .... 
## virus
## on file SM-7IL1E.tar with type VIR .... 
## virus
## on file SM-7KPUR.tar with type VIR .... 
## virus
## on file SM-7MOQJ.tar with type VIR .... 
## virus
## on file CSM5MCUQ_P.fastq.gz with type P .... 
## mtg
## on file CSM5MCUS_P.fastq.gz with type P .... 
## mtg
## on file CSM5MCUW_P.fastq.gz with type P .... 
## mtg
## on file CSM5MCUY_P.fastq.gz with type P .... 
## mtg
## on file CSM5MCY4.tar with type STD .... 
## mtg
## on file CSM5MCY8.tar with type STD .... 
## mtg
## on file CSM67UE3.tar with type STD .... 
## mtg
## on file CSM67UE7.tar with type STD .... 
## mtg
## on file CSM67UEA.tar with type STD .... 
## mtg
## on file CSM67UEM.tar with type STD .... 
## mtg
## on file CSM67UEI.tar with type STD .... 
## mtg
## on file CSM79HO1.tar with type STD .... 
## mtg
## on file SM-6WOCE.tar with type VIR .... 
## virus
## on file SM-6ZEVI.tar with type VIR .... 
## virus
## on file SM-7AA2E.tar with type VIR .... 
## virus
## on file SM-7BP5L.tar with type VIR .... 
## virus
## on file SM-7I1G8.tar with type VIR .... 
## virus
## on file SM-7K25Z.tar with type VIR .... 
## virus
```



## SRA exploration
### Chunk SRA
Read in SRA results:


```r
sra <- clean_names(readr::read_csv("csvs_and_other_metadata/sra_result2979.csv"))
```

```
## Rows: 2979 Columns: 17
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
sra <- sra %>% 
  filter(organism_name == "human gut metagenome" | organism_name == "viral metagenome") %>% 
  filter(library_strategy == "WGS") %>% 
  relocate(sample_accession, library_name, organism_name, experiment_title) %>% 
  arrange(desc(organism_name))

ihmp_sras <- sra %>% select(sample_accession)
write_csv(ihmp_sras, "csvs_and_other_metadata/sras_mgx_vrx_2041.csv")

head(sra)
```

```
## # A tibble: 6 × 17
##   sample_accession library_name organism_name  experiment_title experiment_acce…
##   <chr>            <chr>        <chr>          <chr>            <chr>           
## 1 SRS2907573       C3037C5_MVX  viral metagen… Stool sample C3… SRX3641375      
## 2 SRS2907572       P6010C22_MVX viral metagen… Stool sample P6… SRX3641374      
## 3 SRS2907571       M2084C5_MVX  viral metagen… Stool sample M2… SRX3641373      
## 4 SRS2907570       M2069C14_MVX viral metagen… Stool sample M2… SRX3641372      
## 5 SRS2907568       C3028C17_MVX viral metagen… Stool sample C3… SRX3641371      
## 6 SRS2907569       P6018C17_MVX viral metagen… Stool sample P6… SRX3641370      
## # … with 12 more variables: instrument <chr>, submitter <chr>,
## #   study_accession <chr>, study_title <chr>, sample_title <lgl>,
## #   total_size_mb <dbl>, total_ru_ns <dbl>, total_spots <dbl>,
## #   total_bases <dbl>, library_strategy <chr>, library_source <chr>,
## #   library_selection <chr>
```

```r
sra_cols<- names(sra)
```







does `all_samples$external_id` merge on a trimmed down `sra$library_name`!?!?!


```r
#sum_all_samples <- all_samples %>% select(project, participant_id, week_num, data_type, site_sub_coll, external_id)



sra_mod <- sra %>% 
  separate(library_name, into = c("library_name", NA), sep = "_MGX", remove = T) %>% 
  separate(library_name, into = c("library_name", NA), sep = "_MVX", remove = T)
sra_mod

sra_mgx <- sra %>% 
  filter(organism_name == "human gut metagenome")
  
    


### use a split sra$experiment_title and ihmp$external_id to 
sra_test <- sra_mgx %>% 
  separate(experiment_title, into = c("key", NA), sep = " stool")

#ihmp_mod <- all_samples %>% 
#  select(project, participant_id, week_num, data_type, site_sub_coll, external_id) %>% 
#  separate(project, into = c("project", NA), sep = "_MGX", remove = T)%>% 
#  separate(project, into = c("project", NA), sep = "_MVX", remove = T)

ihmp_mod <- abx_pairs_df %>% 
  select(project, participant_id, week_num, data_type, site_sub_coll, external_id) %>% 
  separate(project, into = c("project", NA), sep = "_MGX", remove = T)%>% 
  separate(project, into = c("project", NA), sep = "_MVX", remove = T)




ihmp_mod_2 <- ihmp_mod %>% 
  separate(external_id, into = c(external_id_trim,""), sep = "",remove = T)



  
left_join(ihmp_mod, sra_mod, by = c("project" = "library_name")) %>%
  drop_na(sample_accession) %>% 
  relocate(sample_accession, project, organism_name, experiment_title)
```


### Chuck more SRA

does `all_samples$project` merge on `sra$library_name`!?!?! only on viromes


```r
sum_all_samples <- all_samples %>% select(project, participant_id, week_num, data_type, site_sub_coll, external_id)

left_join(sum_all_samples, sra, by = c("project" = "library_name")) %>%
  drop_na(sample_accession) %>% 
  relocate(sample_accession, project, organism_name, experiment_title)
#%>% 
#  pull()
```



### end

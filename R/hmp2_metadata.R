library(readr)
library(dplyr)
library(stringr)


metadata <- readr::read_csv("R/VirHMP/hmp2_metadata.csv")
viral_metadata <- filter(metadata, !is.na(reads_viral)) 

barplot(table(metadata$data_type))

colnames(metadata)
table(metadata$data_type)

viral_metadata <- filter(metadata, data_type=="viromics") 
str_detect(string = colnames(viral_metadata), pattern = "ntibiot")

colnames(viral_metadata)[str_detect(string = colnames(viral_metadata), pattern = "ntibiot")]

tmp <- viral_metadata%>%
  select(`External ID`, `Participant ID`, week_num, colnames(viral_metadata)[str_detect(string = colnames(viral_metadata), pattern = "ntibiot")])

tmp <- viral_metadata%>%
  select(`External ID`, `Participant ID`, week_num, Antibiotics, Lomotil, `Dipentum (olsalazine)`, `Cipro (Ciprofloxin)`)


tmp <- metadata%>%
  select(`External ID`, `Participant ID`, week_num, Antibiotics, Lomotil, `Dipentum (olsalazine)`, `Cipro (Ciprofloxin)`)


"Lomotil"                                                                         
[318] "Dipentum (olsalazine)"                                                           
[319] "Reason for stopping Dipentum:"                                                   
[320] "Rowasa enemas (mesalamine enemas)"                                               
[321] "Reason for stopping rowasa enemas:"                                              
[322] "Canasa suppositories (mesalamine suppositories)"                                 
[323] "Reason for stopping canasa suppositories:"                                       
[324] "Flagyl (Metronidazole)"                                                          
[325] "Reason for stopping Flagyl:"                                                     
[326] "Cipro (Ciprofloxin)"                                                             
[327] "Reason for stopping Cipro:"                                                      
[328] "Xifaxin (rifaxamin)"                                                             
[329] "Reason for stopping Xifaxin:"                                                    
[330] "Levaquin"                                                                        
[331] "Reason for stopping Levaquin:"  
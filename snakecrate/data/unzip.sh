
while IFS="," read -r file prefix zip_type uid
do
echo "on file $file... "
  case $zip_type in
    "VIR")
      tar xf raw_data/unaltered_raw_files/$file
      rm     raw_data/unaltered_raw_files/$file
      
      bzip2 -d ${prefix}_1_sequence.txt.bz2
      gzip -c  ${prefix}_1_sequence.txt.bz2 > raw_data/gzips/${prefix}_1.fq.gz
      rm ${prefix}_1_sequence.txt
      
      bzip2 -d ${prefix}_2_sequence.txt.bz2
      gzip -c  ${prefix}_2_sequence.txt.bz2 > raw_data/gzips/${prefix}_2.fq.gz
      rm ${prefix}_2_sequence.txt
      
      fastp\
          --in1  raw_data/gzips/${prefix}_1.fq.gz \
          --in2  raw_data/gzips/${prefix}_2.fq.gz \
          --out1 fastp/${prefix}_1.trim.fq.gz\
          --out2 fastp/${prefix}_2.trim.fq.gz\
          --detect_adapter_for_pe \
          --qualified_quality_phred 4 \
          --length_required 31 --correction \
          --json fastp/reports/${prefix}.trim.json \
          --html fastp/reports/${prefix}.trim.html      

      interleave-reads.py fastp/${prefix}_1.trim.fq.gz fastp/${prefix}_2.trim.fq.gz | \
        trim-low-abund.py --gzip -C 3 -Z 18 -M 20e9 -V - -o abundtrim/${prefix}.abundtrim.fq.gz
      echo "VIR done"

      ;;
    "P"|"TR"|"STD")
      echo "not yet"
      ;;
  esac
  
done < <(tail -n +2 ../../R/csvs_and_other_metadata/file_zip_uid.csv)
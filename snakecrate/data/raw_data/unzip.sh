while IFS="," read -r file zip_type uid
do
echo "on file $file with type $zip_type .... "
  case $zip_type in
    "VIR")
      tar xf $file
      rm $file
      
      bzip2 -d ${prefix}_1_sequence.txt.bz2
      gzip -c  ${prefix}_1_sequence.txt.bz2 > ${prefix}_1.fq.gz
      rm ${prefix}_1_sequence.txt
      
      bzip2 -d ${prefix}_2_sequence.txt.bz2
      gzip -c  ${prefix}_2_sequence.txt.bz2 > ${prefix}_2.fq.gz
      rm ${prefix}_2_sequence.txt
      
      fastp\
          --in1 ../raw_data/${prefix}_1.fq.gz \
          --in2 ../raw_data/${prefix}_2.fq.gz \
          --out1 ../raw_data/${prefix}_1.trim.fastq.gz\
          --out2 ../raw_data/${prefix}_2.trim.fastq.gz\
          --detect_adapter_for_pe \
          --qualified_quality_phred 4 \
          --length_required 31 --correction \
          --json SRR1976948.trim.json \
          --html SRR1976948.trim.html      

      interleave-reads.py ../trim/SRR1976948_1.trim.fastq.gz ../trim/SRR1976948_2.trim.fastq.gz | \
        trim-low-abund.py --gzip -C 3 -Z 18 -M 20e9 -V - -o SRR1976948.abundtrim.fq.gz
           
      
      fastp --in1 ${prefix}_1.fq.gz --in2 ${prefix}_2.fq.gz
      ;;
    "P"|"TR"|"STD")
      echo "mtg"
      ;;
  esac
  
  fastp --in1 ../raw_data/SRR1976948_1.fastq.gz \
  --in2 ../raw_data/SRR1976948_2.fastq.gz \
  --out1 SRR1976948_1.trim.fastq.gz \
  --out2 SRR1976948_2.trim.fastq.gz \
  --detect_adapter_for_pe \
  --qualified_quality_phred 4 \
  --length_required 31 --correction \
  --json SRR1976948.trim.json \
  --html SRR1976948.trim.html
  
  
  
done < <(tail -n +2 ../R/csvs_and_other_metadata/file_zip_uid.csv)
while IFS= read -r file
do	
    #echo "suffix: $file, old: SM-$file.abundtrim.fq.gz, new: SM_$file.abundtrim.fq.gz"
    #
    mv abundtrim/SM-${file}.abundtrim.fq.gz abundtrim/SM_${file}.abundtrim.fq.gz 
    echo: "moved SM-$file.abundtrim.fq.gz to SM_$file.abundtrim.fq.gz"
done < s5_vir_suffixes.txt
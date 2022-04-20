while IFS= read -r file
do
	cp ../../../snakecrate/data/raw_data/tr/abundtrim/$file.abundtrim.fq.gz abundtrim/
    echo "$file"
done < "s5_tr_samplenames.txt"

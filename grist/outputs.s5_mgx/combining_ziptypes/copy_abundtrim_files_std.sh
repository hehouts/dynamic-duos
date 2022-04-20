while IFS= read -r file
do
	cp ../../../snakecrate/data/raw_data/std/abundtrim/$file.abundtrim.fq.gz abundtrim/
    echo "$file"
done < "s5_std_samplenames.txt"

while IFS= read -r file
do
	#cp ../../snakecrate/data/raw_data/vir/abundtrim/$file.abundtrim.fq.gz abundtrim/
    echo "file"
done < "s5_vir_samplenames.txt"
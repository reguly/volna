#! bin/bash

for i in bathy*.txt; do 
	numStr=$(echo ${i} | tr -d '[A-Za-z.]'); 
	num=$((numStr*2))
	new=$(printf "bathy%04d.txt" ${num}); 
	echo "Renaming ${i} to ${new}"
	mv ${i} ${new}
done

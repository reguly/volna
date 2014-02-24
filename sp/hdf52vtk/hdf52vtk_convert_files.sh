#! /bin/bash

#This tool takes all the files in the directory and converts

echo "Converting HDF5 datafiles to VTK files..."

if [ -z $1 ]; then
	echo "Please specify the geometry h5 file!"
	exit 1
fi

for i in *.h5; do
	# If the geometry filename is different from the snapshot filename, convert the file
	if [ $i != $1 ]; then
		if [ -n `echo ${i} | grep '\.h5$'` ]; then
			newname=`echo ${i} | sed -e 's/\.h5$/\.vtk/g'`
       		 	echo "Converting ${i} to ${newname} ..."
			if [ -z $2 ]; then
				./hdf52vtk $1 $i 0 $newname 
			else
				./hdf52vtk $1 $i $2 $newname 
			fi
		fi
	fi
done


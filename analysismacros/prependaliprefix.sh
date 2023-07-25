#!/bin/bash

declare prefix="alien:\/\/"
echo "Prepending prefix ${prefix}"
sleep 1;

#define the file with the list of files to download
databasefile=$1
if [ -z "${databasefile}" ]; 
then
  databasefile="myfile";
  echo "Input file not provided! using the standard myfile"
fi

#check if the input file exists
if [ -s ${databasefile} ];
then
  echo "Input file ${databasefile} is found!";
else
  echo "Input file ${databasefile} does not exist!!! aborting";
  return; 
fi

outputfile="p_${databasefile}"
cat ${databasefile} | sed "s/^/${prefix}/" > "${outputfile}"
echo 'Done prepending!'
sleep 2;
cat "${outputfile}"

echo "#################################"
echo "To download from Grid please run:"
echo "#################################"

echo "alisync @${outputfile}"

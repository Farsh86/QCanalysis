#!/bin/bash

# $1 is run number

root -l -b -q "readdigitaltree_summary.C(\"filemerge_0.dat\",\"$1\", \"DIG SCAN - run$1\",\"L0_05\",0)" &
root -l -b -q "readdigitaltree_summary.C(\"filemerge_1.dat\",\"$1\", \"DIG SCAN - run$1\",\"L0_05\",1)" &
root -l -b -q "readdigitaltree_summary.C(\"filemerge_2.dat\",\"$1\", \"DIG SCAN - run$1\",\"L0_05\",2)" &
root -l -b -q "readdigitaltree_summary.C(\"filemerge_3.dat\",\"$1\", \"DIG SCAN - run$1\",\"L0_05\",3)" &
root -l -b -q "readdigitaltree_summary.C(\"filemerge_4.dat\",\"$1\", \"DIG SCAN - run$1\",\"L0_05\",4)" &
root -l -b -q "readdigitaltree_summary.C(\"filemerge_5.dat\",\"$1\", \"DIG SCAN - run$1\",\"L0_05\",5)" &
root -l -b -q "readdigitaltree_summary.C(\"filemerge_6.dat\",\"$1\", \"DIG SCAN - run$1\",\"L0_05\",6)" &
root -l -b -q "readdigitaltree_summary.C(\"filemerge_7.dat\",\"$1\", \"DIG SCAN - run$1\",\"L0_05\",7)" &
root -l -b -q "readdigitaltree_summary.C(\"filemerge_8.dat\",\"$1\", \"DIG SCAN - run$1\",\"L0_05\",8)" &
root -l -b -q "readdigitaltree_summary.C(\"filemerge_9.dat\",\"$1\", \"DIG SCAN - run$1\",\"L0_05\",9)" &
root -l -b -q "readdigitaltree_summary.C(\"filemerge_10.dat\",\"$1\", \"DIG SCAN - run$1\",\"L0_05\",10)" &
root -l -b -q "readdigitaltree_summary.C(\"filemerge_11.dat\",\"$1\", \"DIG SCAN - run$1\",\"L0_05\",11)" &
root -l -b -q "readdigitaltree_summary.C(\"filemerge_12.dat\",\"$1\", \"DIG SCAN - run$1\",\"L0_05\",12)" &
root -l -b -q "readdigitaltree_summary.C(\"filemerge_13.dat\",\"$1\", \"DIG SCAN - run$1\",\"L0_05\",13)" &
root -l -b -q "readdigitaltree_summary.C(\"filemerge_14.dat\",\"$1\", \"DIG SCAN - run$1\",\"L0_05\",14)" &
root -l -b -q "readdigitaltree_summary.C(\"filemerge_15.dat\",\"$1\", \"DIG SCAN - run$1\",\"L0_05\",15)" &
root -l -b -q "readdigitaltree_summary.C(\"filemerge_16.dat\",\"$1\", \"DIG SCAN - run$1\",\"L0_05\",16)" &
root -l -b -q "readdigitaltree_summary.C(\"filemerge_17.dat\",\"$1\", \"DIG SCAN - run$1\",\"L0_05\",17)" &
root -l -b -q "readdigitaltree_summary.C(\"filemerge_18.dat\",\"$1\", \"DIG SCAN - run$1\",\"L0_05\",18)" &
root -l -b -q "readdigitaltree_summary.C(\"filemerge_19.dat\",\"$1\", \"DIG SCAN - run$1\",\"L0_05\",19)" &
root -l -b -q "readdigitaltree_summary.C(\"filemerge_20.dat\",\"$1\", \"DIG SCAN - run$1\",\"L0_05\",20)" &
root -l -b -q "readdigitaltree_summary.C(\"filemerge_21.dat\",\"$1\", \"DIG SCAN - run$1\",\"L0_05\",21)" &
root -l -b -q "readdigitaltree_summary.C(\"filemerge_22.dat\",\"$1\", \"DIG SCAN - run$1\",\"L0_05\",22)" &
root -l -b -q "readdigitaltree_summary.C(\"filemerge_23.dat\",\"$1\", \"DIG SCAN - run$1\",\"L0_05\",23)" &
root -l -b -q "readdigitaltree_summary.C(\"filemerge_24.dat\",\"$1\", \"DIG SCAN - run$1\",\"L0_05\",24)" &
root -l -b -q "readdigitaltree_summary.C(\"filemerge_25.dat\",\"$1\", \"DIG SCAN - run$1\",\"L0_05\",25)" &
root -l -b -q "readdigitaltree_summary.C(\"filemerge_26.dat\",\"$1\", \"DIG SCAN - run$1\",\"L0_05\",26)"

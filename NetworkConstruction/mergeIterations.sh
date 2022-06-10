#!/usr/bin/bash


# [Number of Iterations]

for I in `seq $1`
	do
		cat *$I"_AllEpiData.tsv" >"Iteration_"$I"_merged.tsv"
	done



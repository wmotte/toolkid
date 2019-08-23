#!/bin/bash
for i in threshold*; do cut -d "," -f 1 $i > clusters_${i}.txt;done 
for i in threshold*; do cut -d "," -f 2 $i > modularity_${i}.txt;done 
for i in threshold*; do cut -d "," -f 3 $i > entropy_${i}.txt;done 
paste clusters_* modularity_threshold_* > all_mod.txt
paste clusters_* entropy_threshold_* > all_entropy.txt

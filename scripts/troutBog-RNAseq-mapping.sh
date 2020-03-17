#! /usr/bin/bash

while read -r a b c; do
	echo $c $a $b;
	kallisto quant -i troutBog-index -o "$c" "$a" "$b";
done < troutBog-geodes-metaTs-queue.txt

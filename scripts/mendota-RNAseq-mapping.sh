#! /usr/bin/bash

while read -r a b c; do
	echo $c $a $b;
	kallisto quant -i mendota-index -o "$c" "$a" "$b";
done < mendota-geodes-metaTs-queue.txt

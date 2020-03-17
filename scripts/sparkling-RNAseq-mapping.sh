#! /usr/bin/bash

while read -r a b c; do
	echo $c $a $b;
	kallisto quant -i sparkling-index -o "$c" "$a" "$b";
done < sparkling-geodes-metaTs-queue.txt

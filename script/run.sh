#!/bin/bash

for mode in slow fast
do

	for name in DMA pci_bridge32 des_perf vga_lcd b19 leon3mp netcard
	#for name in DMA
	do
		bench=$name\_$mode
	
		echo "Running $bench..."
		./usizer ~/benchmarks/sizing/ispd2012/ $bench $1 | tee $bench.out
	done

done


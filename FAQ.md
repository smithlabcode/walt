## Frequently Asked Questions ##

1. Why does WALT get much lower mapping rate (mapping effiency) compared to other bisulfite mappers?
	
	By default, WALT ignores reads shorter than 36. If you get lower mapping rate, please go to mapstats file and check the 			percentage of reads that are shorter than 36. If a large portion of reads are shorter than 36, I suggest you use WALT seed 	pattern 7.

	(1) Open walt/src/walt/Makefile, modify “-D SEEDPATTERN3" to " -D SEEDPATTERN7"

	(2) Rebuild the program using `make all`, `make install`
	
	(3) Rebuild your index using the new `makedb` command
	
	(4) Map reads using the new `walt` command
	
	
	
	


smacofSSUO.c ../common/smacofPrint.c ../common/smacofUtils.c
	clang -o smacofSSUO smacofSSUO.c ../common/smacofPrint.c ../common/smacofUtils.c

rlib: 
	R CMD SHLIB smacofSSUO.c ../common/smacofPrint.c ../common/smacofUtils.c

clean:
	rm -rf *.o ../common/*.o core a.out

pristine:
	rm -rf smacofSSUO smacofSSUO.so

format:
	clang-format -i *.c


smacofSSUP: smacofSSUP.c ../common/smacofPrint.c ../common/smacofUtils.c
	clang -o smacofSSUP smacofSSUP.c ../common/smacofPrint.c ../common/smacofUtils.c

rlib: 
	R CMD SHLIB smacofSSUP.c ../common/smacofPrint.c ../common/smacofUtils.c

clean:
	rm -rf *.o ../common/*.o core a.out

pristine:
	rm -rf smacofSSUP smacofSSUP.so

format:
	clang-format -i *.c


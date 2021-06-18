smacofSSWR: smacofSSWI.c ../common/smacofPrint.c ../common/smacofUtils.c
	clang -o smacofSSWI smacofSSWI.c ../common/smacofPrint.c ../common/smacofUtils.c

rlib: 
	R CMD SHLIB smacofSSWI.c ../common/smacofPrint.c ../common/smacofUtils.c

clean:
	rm -rf *.o ../common/*.o core a.out

pristine:
	rm -rf smacofSSWI smacofSSWI.so

format:
	clang-format -i *.c


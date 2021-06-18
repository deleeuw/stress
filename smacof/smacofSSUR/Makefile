smacofSSUR: smacofSSUR.c ../common/smacofPrint.c ../common/smacofUtils.c
	clang -o smacofSSUR smacofSSUR.c ../common/smacofPrint.c ../common/smacofUtils.c

rlib: 
	R CMD SHLIB smacofSSUR.c ../common/smacofPrint.c ../common/smacofUtils.c

clean:
	rm -rf *.o ../common/*.o core a.out

pristine:
	rm -rf smacofSSUR smacofSSUR.so

format:
	clang-format -i *.c


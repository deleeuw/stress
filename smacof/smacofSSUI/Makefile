smacofSSUI: smacofSSUI.c ../common/smacofPrint.c ../common/smacofUtils.c
	clang -o smacofSSUI smacofSSUI.c ../common/smacofPrint.c ../common/smacofUtils.c

rlib: 
	R CMD SHLIB smacofSSUI.c ../common/smacofPrint.c ../common/smacofUtils.c

clean:
	rm -rf *.o ../common/*.o core a.out

pristine:
	rm -rf smacofSSUI smacofSSUI.so

format:
	clang-format -i *.c



Using GNU compilers: gcc and gfortran

$ gcc -c c-call-f.c
$ gfortran -ffree-form -c f_sub.f
$ gcc -o c-call-f-gcc  c-call-f.o f_sub.o  -lgfortran
$ ./c-call-f-gcc/Volumes/myOrange/Users/deleeuw/myStuff/myCourses/R_fli/05_dotFortran/sweep/c_call_f/script.sh
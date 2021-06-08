Suppose foo.c contains

void foo (void) {
	int i = 2;
	int j = i * i;
}

Compile with

gcc -fPIC -dynamiclib -o foo.so foo.c -g

Suppose main.c contains

int main(void) {
	foo ();
}

Compile with

gcc -o main main.c ./foo.so -g

Then say 

/usr/bin/gdb main

Now if you say 

(gdb) b foo
(gdb) r 

you can step through the function foo from the shared
library

	
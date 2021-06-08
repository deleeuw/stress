
gcc -fpic -shared -o foo.so foo.c -g

int main(void) {
	foo ();
}
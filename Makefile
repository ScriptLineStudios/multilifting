main: example.c	
	gcc -g -o example example.c cubiomes/libcubiomes.a -Wall -Wextra -lm -Ofast -Wno-missing-field-initializers
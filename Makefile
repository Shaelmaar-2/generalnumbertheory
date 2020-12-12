prime: primes.c
	gcc -o prime primes.c -lgmp
debug: primes.c
	gcc -g -o prime primes.c -lgmp

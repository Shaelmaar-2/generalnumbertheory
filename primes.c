#include <gmp.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <stdbool.h>
#include <string.h>

int myPow(int x,int n)
{
    int i; /* Variable used in loop counter */
    int number = 1;

    for (i = 0; i < n; ++i)
        number *= x;

    return(number);
}


void fastModExp(mpz_t a, mpz_t k, mpz_t m, mpz_t out){
    mpz_set_ui(out, 1);
    mpz_t temp, ta, tk;
    mpz_init_set_ui(temp, 0);
	mpz_init_set(ta, a);
	mpz_init_set(tk, k);
    
    while (mpz_cmp_si(tk, 1) >= 0) {	 // while k >= 1
        if (mpz_congruent_ui_p(tk, 1, 2) != 0) {  // if k % 2 == 1
            mpz_mul(temp,ta,out);
            mpz_mod(out,temp,m); // out = (out * a) % m
        }
        mpz_mul(temp,ta,ta);
        mpz_mod(ta,temp,m);	// a = (a**2) % m
        mpz_tdiv_q_ui(tk, tk, 2);  // k //=2
    }
    mpz_clear(temp);
	mpz_clear(ta);
	mpz_clear(tk);
}

bool rabinMiller(mpz_t a, int iters) {
	gmp_randstate_t state;
	gmp_randinit_mt(state);
	unsigned long int seed = (unsigned long int) time(NULL);
	gmp_randseed_ui(state, seed);
	mpz_t q, tq, k, x, r, mod, neg;

	mpz_init(mod);
	mpz_init(x);
	mpz_init(q);
	mpz_init(r);
	mpz_init(neg);
	mpz_init(tq);
	mpz_init_set_ui(k, 0);

    mpz_sub_ui(r, a, 3);
	mpz_sub_ui(neg, a, 1);
    mpz_sub_ui(q, a, 1); // q = a - 1

	while(mpz_congruent_ui_p(q, 0, 2) != 0){ // while q % 2 == 0
		mpz_add_ui(k, k, 1);  // k += 1
		mpz_divexact_ui(q, q, 2);  // q /= 2
	}

	unsigned long int ik = mpz_get_ui(k);
	bool brk_flg;

	for(int i = 0; i < iters; i++) {
		brk_flg = false;
		mpz_urandomm(x, state, r);  // set x to a random value in [0,n-4]
		mpz_add_ui(x, x, 2);  // add 2 to be in the range [2,n-2]
		fastModExp(x, q, a, mod);
		if(mpz_cmp_si(mod, 1) != 0){
			for(int j = 0; j < ik; j++){
				mpz_mul_ui(tq, q, myPow(2,j));
				fastModExp(x,tq,a,mod);
				if(mpz_cmp(mod, neg) == 0){
					brk_flg = true;
					break;
				}
			}
			if(brk_flg == false){
			    mpz_clear(q);
   				mpz_clear(tq);
    			mpz_clear(neg);
    			mpz_clear(x);
    			mpz_clear(k);
    			mpz_clear(mod);
    			mpz_clear(r);
				gmp_randclear(state);
				return false;
			}
		}
	}
	mpz_clear(q);
	mpz_clear(tq);
	mpz_clear(neg);
	mpz_clear(x);
	mpz_clear(k);
	mpz_clear(mod);
	mpz_clear(r);
	gmp_randclear(state);
	return true;
}

void newPrime(mpz_t out, int iters, unsigned int bitCnt){
    gmp_randstate_t state;
    gmp_randinit_mt(state);
    unsigned long int seed = (unsigned long int) time(NULL);
    gmp_randseed_ui(state, seed);
	mpz_urandomb(out, state, bitCnt);
	while (!rabinMiller(out, iters)){
		mpz_urandomb(out, state, bitCnt);
	}
	gmp_randclear(state);
}


int main(int argc, char** argv) {
	bool newFlag = false;
	bool helpFlag = false;
    for(int i = 0; i < argc; i++){
		char *arg = argv[i];
		if('-' == arg[0]){
			if(strncmp(arg,"-n",2) == 0) {
				newFlag = true;
			}
			if(strncmp(arg,"-h",2) == 0){
				printf("Usage: \n");
				printf("for new prime: -n bitlen iters (iters of rabin miller)\n");
			}
		}
	}
	if(newFlag){
		mpz_t prm;
		mpz_init(prm);

		int iters, bits;
		bits = atoi(argv[2]);
		iters = atoi(argv[3]);
		newPrime(prm, iters, bits);

		mpz_out_str(stdout, 10, prm);
		printf("\n");

		mpz_clear(prm);
	}

}

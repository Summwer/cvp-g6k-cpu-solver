#include <stdio.h>
#include <mpfr.h>


int main(){
    mpfr_t x;
    mpfr_get_default_prec();
    mpfr_init2(x,256);
    mpfr_set_d(x,3.14,MPFR_RNDN);
    mpfr_printf("x = %.10Rf\n",x);
}

/* arithmetic_operations.c */

/*
 * Original Fortran90 version written by Henry Neeman, University of Oklahoma
 *
 * This program was converted from a Fortran program used
 * in the OU supercomputing course, "Supercomputing in Plain
 * English" offered between February and April 2009.  This
 * C version and the program was created by Edward "Ted" Doyle
 * Ph.D., between Feb 17 and 22, 2009.
 *
 * It was then modified by Henry Neeman, March 3 2009 and Feb 21 2011.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>

#include "second.h"

static void set_permutation(int permute[], int length);
static void check_permutation(int permute[], int length);
static void initialize_source_array(
                float rsrc1[], float rsrc2[], float rsrc3[], float rsrc4[],
                int   isrc1[], int   isrc2[], int   isrc3[], int   isrc4[],
                int length);
static void time_of_add_real_array(
                float dst[], float src1[], float src2[],
                int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_add_int_array(
                int dst[], int src1[], int src2[],
                int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_sum_real_array(
                float *dst, float src[],
                int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_sum_int_array(
                int dst[], int src[],
                int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_sub_real_array(
                float rdst[], float rsrc1[], float rsrc2[],
                int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_sub_int_array(
                int dst[], int src1[], int src2[],
                int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_mul_real_array(
                float dst[], float src1[], float src2[],
                int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_mul_int_array(
                int dst[], int src1[], int src2[],
                int length,
                double * time_of_operation, double *mflops_of_operation);
static void time_of_mad_real_array(
                float dst[], float src1[], float src2[],
                int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_mad_int_array(
                int dst[], int src1[], int src2[],
                int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_mam_real_array(
                float dst[],
                float src1[], float src2[], float src3[], float src4[],
                int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_mam_int_array(
                int dst[],
                int src1[], int src2[], int src3[], int src4[],
                int length,
                double  *time_of_operation, double *mflops_of_operation);
static void time_of_div_real_array(
                float dst[], float src1[], float src2[],
                int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_div_int_array(
                int dst[], int src1[], int src2[],
                int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_pow_real_array(
                float  dst[], float src1[], float src2[],
                int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_mod_int_array(
                int dst[], int src1[], int src2[],
                int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_sqrt_real_array(
                float dst[], float src[],
                int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_cos_real_array(
                float dst[], float src[],
                int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_exp_real_array(
                float dst[], float src[],
                int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_log_real_array(
                float dst[], float src[],
                int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_euc_real_array(
                float dst[],
                float src1[], float src2[], float src3[], float src4[],
                int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_dot_real_array(
                float *dot, float src1[], float src2[],
                int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_lot8_real_array(
                float *lot,
                float src1[], float src2[], float src3[], float src4[],
                int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_lot10_real_array(
                float *lot,
                float src1[], float src2[], float src3[], float src4[],
                int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_lot12_real_array(
                float *lot,
                float src1[], float src2[], float src3[], float src4[],
                int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_lot16_real_array(
                float *lot,
                float src1[], float src2[], float src3[], float src4[],
                int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_lot20_real_array(
                float *lot,
                float src1[], float src2[], float src3[], float src4[],
                int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_lot24_real_array(
                float *lot,
                float src1[], float src2[], float src3[], float src4[],
                int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_i2r_array(
                float rdst[], int isrc[],
                int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_r2i_array(
                int idst[], float rsrc[],
                int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_add_real_array_perm(
                float dst[], float src1[], float src2[],
                int permute[], int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_add_int_array_perm(
                int dst[], int src1[], int src2[],
                int permute[], int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_sum_real_array_perm(
                float *dst, float src[],
                int permute[], int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_sum_int_array_perm(
                int dst[], int src[],
                int permute[], int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_sub_real_array_perm(
                float dst[], float src1[], float src2[],
                int permute[], int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_sub_int_array_perm(
                int dst[], int src1[], int src2[],
                int permute[], int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_mul_real_array_perm(
                float dst[], float src1[], float src2[],
                int permute[], int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_mul_int_array_perm(
                int dst[], int src1[], int src2[],
                int permute[], int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_mad_real_array_perm(
                float dst[], float src1[], float src2[],
                int permute[], int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_mad_int_array_perm(
                int dst[], int src1[], int src2[],
                int permute[], int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_mam_real_array_perm(
                float dst[],
                float src1[], float src2[], float src3[], float src4[],
                int permute[], int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_mam_int_array_perm(
                int dst[], int src1[], int src2[],
                int src3[], int src4[],
                int permute[], int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_div_real_array_perm(
                float dst[], float src1[], float src2[],
                int permute[], int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_div_int_array_perm(
                int dst[], int src1[], int src2[],
                int permute[], int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_pow_real_array_perm(
                float dst[], float src1[], float src2[],
                int permute[], int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_mod_int_array_perm(
                int dst[], int src1[], int src2[],
                int permute[], int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_sqrt_real_array_perm(
                float dst[], float src[],
                int permute[], int length,
                double * time_of_operation, double *mflops_of_operation);
static void time_of_cos_real_array_perm(
                float dst[], float src[],
                int permute[], int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_exp_real_array_perm(
                float dst[], float src[],
                int permute[], int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_log_real_array_perm(
                float dst[], float src[],
                int permute[], int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_euc_real_array_perm(
                float dst[],
                float src1[], float src2[], float src3[], float src4[],
                int permute[], int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_dot_real_array_perm(
                float *dot, float src1[], float src2[],
                int permute[], int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_lot8_real_array_perm(
                float *lot,
                float src1[], float src2[], float src3[], float src4[],
                int permute[], int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_lot10_real_array_perm(
                float *lot,
                float src1[], float src2[], float src3[], float src4[],
                int permute[], int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_lot12_real_array_perm(
                float *lot,
                float src1[], float src2[], float src3[], float src4[],
                int permute[], int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_lot16_real_array_perm(
                float *lot,
                float src1[], float src2[], float src3[], float src4[],
                int permute[], int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_lot20_real_array_perm(
                float *lot,
                float src1[], float src2[], float src3[], float src4[],
                int permute[], int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_lot24_real_array_perm(
                float *lot,
                float src1[], float src2[], float src3[], float src4[],
                int permute[], int length,
                double *time_of_operation, double *mflops_of_operation);
static void time_of_i2r_array_perm(
                float rdst[], int isrc[],
                int permute[],int length,
                double *time_of_operation,  double *mflops_of_operation);
static void time_of_r2i_array_perm(
                int idst[], float rsrc[],
                int permute[], int length,
                double *time_of_operation, double *mflops_of_operation);

int main (int argc, char *argv[])
{ /* main */
    const int maximum_filename_length = 512;
    const int program_failure_code    =  -1;
    const int program_success_code    =   0;

    char output_filename_with_quotes[maximum_filename_length + 2];
    char output_filename[maximum_filename_length];

    FILE*  output_file = (FILE*) NULL;
    float* rdst        = (float*)NULL;
    float* rsrc1       = (float*)NULL;
    float* rsrc2       = (float*)NULL;
    float* rsrc3       = (float*)NULL;
    float* rsrc4       = (float*)NULL;
    int*   idst        = (int*)  NULL;
    int*   isrc1       = (int*)  NULL;
    int*   isrc2       = (int*)  NULL;
    int*   isrc3       = (int*)  NULL;
    int*   isrc4       = (int*)  NULL;
    int*   permute     = (int*)  NULL;

    double add_real_time,     add_int_time,
           sum_real_time,     sum_int_time,
           sub_real_time,     sub_int_time,
           mul_real_time,     mul_int_time,
           mad_real_time,     mad_int_time,
           mam_real_time,     mam_int_time,
           div_real_time,     div_int_time,
           pow_real_time,     mod_int_time,
           sqrt_real_time,    cos_real_time,
           exp_real_time,     log_real_time,
           euc_real_time,     dot_real_time,
           lot8_real_time,    lot10_real_time,
           lot12_real_time,   lot16_real_time,
           lot20_real_time,   lot24_real_time,
           i2r_time,          r2i_time,
           add_real_mflops,   add_int_mflops,
           sum_real_mflops,   sum_int_mflops,
           sub_real_mflops,   sub_int_mflops,
           mul_real_mflops,   mul_int_mflops,
           mad_real_mflops,   mad_int_mflops,
           mam_real_mflops,   mam_int_mflops,
           div_real_mflops,   div_int_mflops,
           pow_real_mflops,   mod_int_mflops,
           sqrt_real_mflops,  cos_real_mflops,
           exp_real_mflops,   log_real_mflops,
           euc_real_mflops,   dot_real_mflops,
           lot8_real_mflops,  lot10_real_mflops,
           lot12_real_mflops, lot16_real_mflops,
           lot20_real_mflops, lot24_real_mflops,
           i2r_mflops,        r2i_mflops,
           add_real_permuted_time,     add_int_permuted_time,
           sum_real_permuted_time,     sum_int_permuted_time,
           sub_real_permuted_time,     sub_int_permuted_time,
           mul_real_permuted_time,     mul_int_permuted_time,
           mad_real_permuted_time,     mad_int_permuted_time,
           mam_real_permuted_time,     mam_int_permuted_time,
           div_real_permuted_time,     div_int_permuted_time,
           pow_real_permuted_time,     mod_int_permuted_time,
           sqrt_real_permuted_time,    cos_real_permuted_time,
           exp_real_permuted_time,     log_real_permuted_time,
           euc_real_permuted_time,     dot_real_permuted_time,
           lot8_real_permuted_time,    lot10_real_permuted_time,
           lot12_real_permuted_time,   lot16_real_permuted_time,
           lot20_real_permuted_time,   lot24_real_permuted_time,
           i2r_permuted_time,          r2i_permuted_time,
           add_real_permuted_mflops,   add_int_permuted_mflops,
           sum_real_permuted_mflops,   sum_int_permuted_mflops,
           sub_real_permuted_mflops,   sub_int_permuted_mflops,
           mul_real_permuted_mflops,   mul_int_permuted_mflops,
           mad_real_permuted_mflops,   mad_int_permuted_mflops,
           mam_real_permuted_mflops,   mam_int_permuted_mflops,
           div_real_permuted_mflops,   div_int_permuted_mflops,
           pow_real_permuted_mflops,   mod_int_permuted_mflops,
           sqrt_real_permuted_mflops,  cos_real_permuted_mflops,
           exp_real_permuted_mflops,   log_real_permuted_mflops,
           euc_real_permuted_mflops,   dot_real_permuted_mflops,
           lot8_real_permuted_mflops,  lot10_real_permuted_mflops,
           lot12_real_permuted_mflops, lot16_real_permuted_mflops,
           lot20_real_permuted_mflops, lot24_real_permuted_mflops,
           i2r_permuted_mflops,        r2i_permuted_mflops;
 
    float dot, lot, rsum;
    int   isum;
    int   length;

    scanf("%d", &length);
    scanf("%s", output_filename_with_quotes);
    strncpy(output_filename, &output_filename_with_quotes[1],
            strlen(output_filename_with_quotes) - 2);

    rdst    = (float*)malloc(length * sizeof(float));
    rsrc1   = (float*)malloc(length * sizeof(float));
    rsrc2   = (float*)malloc(length * sizeof(float));
    rsrc3   = (float*)malloc(length * sizeof(float));
    rsrc4   = (float*)malloc(length * sizeof(float));

    idst    = (int*)  malloc(length * sizeof(int));
    isrc1   = (int*)  malloc(length * sizeof(int));
    isrc2   = (int*)  malloc(length * sizeof(int));
    isrc3   = (int*)  malloc(length * sizeof(int));
    isrc4   = (int*)  malloc(length * sizeof(int));

    permute = (int*)  malloc(length * sizeof(int));

    set_permutation  (permute, length);
    check_permutation(permute, length);

    initialize_source_array(rsrc1, rsrc2, rsrc3, rsrc4,
                            isrc1, isrc2, isrc3, isrc4, length);

    output_file = fopen(output_filename, "w");
    if (output_file == (FILE*)NULL) {
        fprintf(stderr, "ERROR: cannot open output file %s.\n",
            output_filename);
        exit(program_failure_code);
    } /* if (output_file == (FILE*)NULL) */

    setbuf(output_file, (char*)NULL);
    fprintf(output_file,
        "    SIZE             TYPE     op      TIME    MFLOPS    RATIO\n");

    check_permutation(permute, length);
    time_of_add_real_array(rdst, rsrc1, rsrc2, length,
                           &add_real_time, &add_real_mflops);
    time_of_add_int_array (idst, isrc1, isrc2, length,
                           &add_int_time,  &add_int_mflops);

    fprintf(output_file, " %10d          Real     add:  %8.2f  %8.2f \n",
        length, add_real_time, add_real_mflops);
    fprintf(output_file, " %10d          Integer  add:  %8.2f  %8.2f \n",
        length, add_int_time,  add_int_mflops);

    time_of_sum_real_array(&rsum, rsrc3, length,
                           &sum_real_time, &sum_real_mflops);
    time_of_sum_int_array (&isum, isrc3, length,
                           &sum_int_time,  &sum_int_mflops);

    fprintf(output_file, " %10d          Real     sum:  %8.2f  %8.2f \n",
        length, sum_real_time, sum_real_mflops);
    fprintf(output_file, " %10d          Integer  sum:  %8.2f  %8.2f \n",
        length, sum_int_time,  sum_int_mflops);

    time_of_sub_real_array(rdst, rsrc2, rsrc4, length,
                           &sub_real_time, &sub_real_mflops);
    time_of_sub_int_array (idst, isrc2, isrc4, length,
                           &sub_int_time,  &sub_int_mflops);

    fprintf(output_file, " %10d          Real     sub:  %8.2f  %8.2f \n",
        length, sum_real_time, sum_real_mflops);
    fprintf(output_file, " %10d          Integer  sub:  %8.2f  %8.2f \n",
        length, sum_int_time,  sum_int_mflops);

    time_of_mul_real_array(rdst, rsrc1, rsrc3, length,
                           &mul_real_time, &mul_real_mflops);
    time_of_mul_int_array (idst, isrc1, isrc3, length,
                           &mul_int_time,  &mul_int_mflops);

    fprintf(output_file, " %10d          Real     mul:  %8.2f  %8.2f \n",
        length, mul_real_time, mul_real_mflops);
    fprintf(output_file, " %10d          Integer smul:  %8.2f  %8.2f \n",
        length, mul_int_time,  mul_int_mflops);

    time_of_mam_real_array(rdst, rsrc1, rsrc2, rsrc3, rsrc4, length,
                            &mam_real_time, &mam_real_mflops);
    time_of_mam_int_array (idst, isrc1, isrc2, isrc3, isrc4, length,
                            &mam_int_time,  &mam_int_mflops);

    fprintf(output_file, " %10d          Real     mam:  %8.2f  %8.2f \n",
        length, mam_real_time, mam_real_mflops);
    fprintf(output_file, " %10d          Integer  mam:  %8.2f  %8.2f \n",
        length, mam_int_time,  mam_int_mflops);

    time_of_mad_real_array(rdst, rsrc2, rsrc4, length,
                            &mad_real_time, &mad_real_mflops);
    time_of_mad_int_array (idst, isrc2, isrc4, length,
                            &mad_int_time,  &mad_int_mflops);

    fprintf(output_file, " %10d          Real     mad:  %8.2f  %8.2f \n",
        length, mad_real_time, mad_real_mflops);
    fprintf(output_file, " %10d          Integer  mad:  %8.2f  %8.2f \n",
        length, mad_int_time,  mad_int_mflops);

    time_of_div_real_array(rdst, rsrc3, rsrc1, length,
                            &div_real_time, &div_real_mflops);
    time_of_div_int_array (idst, isrc3, isrc1, length,
                            &div_int_time,  &div_int_mflops);

    fprintf(output_file, " %10d          Real     div:  %8.2f  %8.2f \n",
        length, div_real_time, div_real_mflops);
    fprintf(output_file, " %10d          Integer  div:  %8.2f  %8.2f \n",
        length, div_int_time,  div_int_mflops);

    time_of_pow_real_array(rdst, rsrc2, rsrc4, length,
                            &pow_real_time, &pow_real_mflops);
    time_of_mod_int_array (idst, isrc3, isrc1, length,
                            &mod_int_time,  &mod_int_mflops);

    fprintf(output_file, " %10d          Real     pow:  %8.2f  %8.2f \n",
        length, pow_real_time, pow_real_mflops);
    fprintf(output_file, " %10d          Integer  mod:  %8.2f  %8.2f \n",
        length, mod_int_time,  mod_int_mflops);

    time_of_sqrt_real_array(rdst, rsrc3, length,
                             &sqrt_real_time, &sqrt_real_mflops);
    time_of_cos_real_array (rdst, rsrc4, length,
                             &cos_real_time,  &cos_real_mflops);

    fprintf(output_file, " %10d          Real    sqrt:  %8.2f  %8.2f \n",
        length, sqrt_real_time, sqrt_real_mflops);
    fprintf(output_file, " %10d          Real     cos:  %8.2f  %8.2f \n",
        length, cos_real_time,  cos_real_mflops);

    time_of_exp_real_array(rdst, rsrc3, length,
                            &exp_real_time, &exp_real_mflops);
    time_of_log_real_array(rdst, rsrc2, length,
                            &log_real_time, &log_real_mflops);

    fprintf(output_file, " %10d          Real     exp:  %8.2f  %8.2f \n",
        length, sqrt_real_time, sqrt_real_mflops);
    fprintf(output_file, " %10d          Real     log:  %8.2f  %8.2f \n",
        length, cos_real_time, cos_real_mflops);

    time_of_dot_real_array(&dot, rsrc1, rsrc3, length,
                            &dot_real_time, &dot_real_mflops);

    fprintf(output_file, " %10d          Real     dot:  %8.2f  %8.2f \n",
        length, dot_real_time, dot_real_mflops);

    time_of_euc_real_array(rdst, rsrc1, rsrc2, rsrc3, rsrc4, length,
                            &euc_real_time, &euc_real_mflops);

    fprintf(output_file, " %10d          Real     euc:  %8.2f  %8.2f \n",
        length, euc_real_time, euc_real_mflops);

    time_of_lot8_real_array(&lot, rsrc1, rsrc2, rsrc3, rsrc4, length,
                             &lot8_real_time, &lot8_real_mflops);

    fprintf(output_file, " %10d          Real    lot8:  %8.2f  %8.2f \n",
        length, lot8_real_time, lot8_real_mflops);

    time_of_lot10_real_array(&lot, rsrc1, rsrc2, rsrc3, rsrc4, length,
                              &lot10_real_time, &lot10_real_mflops);

    fprintf(output_file, " %10d          Real   lot10:  %8.2f  %8.2f \n",
        length, lot10_real_time, lot10_real_mflops);

    time_of_lot12_real_array(&lot, rsrc1, rsrc2, rsrc3, rsrc4, length,
                              &lot12_real_time, &lot12_real_mflops);

    fprintf(output_file, " %10d          Real   lot12:  %8.2f  %8.2f \n",
        length, lot12_real_time, lot12_real_mflops);

    time_of_lot16_real_array(&lot, rsrc1, rsrc2, rsrc3, rsrc4, length,
                              &lot16_real_time, &lot16_real_mflops);

    fprintf(output_file, " %10d          Real   lot16:  %8.2f  %8.2f \n",
        length, lot16_real_time, lot16_real_mflops);

    time_of_lot20_real_array(&lot, rsrc1, rsrc2, rsrc3, rsrc4, length,
                              &lot20_real_time, &lot20_real_mflops);

    fprintf(output_file, " %10d          Real   lot20:  %8.2f  %8.2f \n",
        length, lot20_real_time, lot20_real_mflops);

    time_of_lot24_real_array(&lot, rsrc1, rsrc2, rsrc3, rsrc4, length,
                              &lot24_real_time, &lot24_real_mflops);

    fprintf(output_file, " %10d          Real   lot24:  %8.2f  %8.2f \n",
        length, lot24_real_time, lot24_real_mflops);

    time_of_i2r_array(rdst, isrc4, length,
                       &i2r_time, &i2r_mflops);
    time_of_r2i_array(idst, rsrc4, length,
                       &r2i_time, &r2i_mflops);

    fprintf(output_file, " %10d                   i2r:  %8.2f  %8.2f \n",
        length, i2r_time, i2r_mflops);
    fprintf(output_file, " %10d                   r2i:  %8.2f  %8.2f \n",
        length, r2i_time, r2i_mflops);

    check_permutation(permute, length);

    time_of_add_real_array_perm(rdst, rsrc1, rsrc2, permute, length,
                                 &add_real_permuted_time,
                                 &add_real_permuted_mflops);
    time_of_add_int_array_perm (idst, isrc1, isrc2, permute, length,
                                 &add_int_permuted_time,
                                 &add_int_permuted_mflops);

    fprintf(output_file, " %10d permuted Real     add:  %8.2f  %8.2f %8.2f\n",
        length, add_real_permuted_time, add_real_permuted_mflops,
        add_real_mflops / add_real_permuted_mflops);
    fprintf(output_file, " %10d permuted Integer  add:  %8.2f  %8.2f %8.2f\n",
        length, add_int_permuted_time, add_int_permuted_mflops,
        add_int_mflops  / add_int_permuted_mflops);

    time_of_sum_real_array_perm(&rsum, rsrc3, permute, length,
                                 &sum_real_permuted_time,
                                 &sum_real_permuted_mflops);
    time_of_sum_int_array_perm (&isum, isrc3, permute, length,
                                 &sum_int_permuted_time,
                                 &sum_int_permuted_mflops);

    fprintf(output_file, " %10d permuted Real     sum:  %8.2f  %8.2f %8.2f\n",
        length, sum_real_permuted_time, sum_real_permuted_mflops,
        sum_real_mflops / sum_real_permuted_mflops);
    fprintf(output_file, " %10d permuted Integer  sum:  %8.2f  %8.2f %8.2f\n",
        length, sum_int_permuted_time, sum_int_permuted_mflops,
        sum_int_mflops / sum_int_permuted_mflops);

    time_of_sub_real_array_perm(rdst, rsrc2, rsrc4, permute, length,
                                 &sub_real_permuted_time,
                                 &sub_real_permuted_mflops);
    time_of_sub_int_array_perm (idst, isrc2, isrc4, permute, length,
                                 &sub_int_permuted_time,
                                 &sub_int_permuted_mflops);

    fprintf(output_file, " %10d permuted Real     sub:  %8.2f  %8.2f %8.2f\n",
        length, sub_real_permuted_time, sub_real_permuted_mflops,
        sub_real_mflops / sub_real_permuted_mflops);
    fprintf(output_file, " %10d permuted Integer  sub:  %8.2f  %8.2f %8.2f\n",
        length, sub_int_permuted_time, sub_int_permuted_mflops,
        sub_int_mflops  / sub_int_permuted_mflops);

    time_of_mul_real_array_perm(rdst, rsrc1, rsrc3, permute, length,
                                 &mul_real_permuted_time,
                                 &mul_real_permuted_mflops);
    time_of_mul_int_array_perm (idst, isrc1, isrc3, permute, length,
                                 &mul_int_permuted_time,
                                 &mul_int_permuted_mflops);

    fprintf(output_file, " %10d permuted Real     mul:  %8.2f  %8.2f %8.2f\n",
        length, mul_real_permuted_time, mul_real_permuted_mflops,
        mul_real_mflops / mul_real_permuted_mflops);
    fprintf(output_file, " %10d permuted Integer  mul:  %8.2f  %8.2f %8.2f\n",
        length, mul_int_permuted_time, mul_int_permuted_mflops,
        mul_int_mflops  / mul_int_permuted_mflops);

    time_of_mad_real_array_perm(rdst, rsrc2, rsrc4, permute, length,
                                 &mad_real_permuted_time,
                                 &mad_real_permuted_mflops);
    time_of_mad_int_array_perm (idst, isrc2, isrc4, permute, length,
                                 &mad_int_permuted_time,
                                 &mad_int_permuted_mflops);

    fprintf(output_file, " %10d permuted Real     mad:  %8.2f  %8.2f %8.2f\n",
        length, mad_real_permuted_time, mad_real_permuted_mflops,
        mad_real_mflops / mad_real_permuted_mflops);
    fprintf(output_file, " %10d permuted Integer  mad:  %8.2f  %8.2f %8.2f\n",
        length, mad_int_permuted_time, mad_int_permuted_mflops,
        mad_int_mflops  / mad_int_permuted_mflops);

    time_of_mam_real_array_perm(rdst, rsrc1, rsrc2, rsrc3, rsrc4,
                                 permute, length,
                                 &mam_real_permuted_time,
                                 &mam_real_permuted_mflops);
    time_of_mam_int_array_perm (idst, isrc1, isrc2, isrc3, isrc4,
                                 permute, length,
                                 &mam_int_permuted_time,
                                 &mam_int_permuted_mflops);

    fprintf(output_file, " %10d permuted Real     mam:  %8.2f  %8.2f %8.2f\n",
        length, mad_real_permuted_time, mad_real_permuted_mflops,
        mam_real_mflops / mam_real_permuted_mflops);
    fprintf(output_file, " %10d permuted Integer  mam:  %8.2f  %8.2f %8.2f\n",
        length, mam_int_permuted_time, mam_int_permuted_mflops,
        mam_int_mflops  / mam_int_permuted_mflops   );

    time_of_div_real_array_perm(rdst, rsrc4, rsrc2, permute, length,
                                 &div_real_permuted_time,
                                 &div_real_permuted_mflops);
    time_of_div_int_array_perm (idst, isrc4, isrc2, permute, length,
                                 &div_int_permuted_time,
                                 &div_int_permuted_mflops);

    fprintf(output_file, " %10d permuted Real     div:  %8.2f  %8.2f %8.2f\n",
        length, div_real_permuted_time, div_real_permuted_mflops,
        div_real_mflops / div_real_permuted_mflops);
    fprintf(output_file, " %10d permuted Integer  div:  %8.2f  %8.2f %8.2f\n",
        length, div_int_permuted_time, div_int_permuted_mflops,
        div_int_mflops  / div_int_permuted_mflops);

    time_of_pow_real_array_perm(rdst, rsrc1, rsrc3, permute, length,
                                 &pow_real_permuted_time,
                                 &pow_real_permuted_mflops);
    time_of_mod_int_array_perm (idst, isrc1, isrc2, permute, length,
                                 &mod_int_permuted_time,
                                 &mod_int_permuted_mflops);

    fprintf(output_file, " %10d permuted Real     pow:  %8.2f  %8.2f %8.2f\n",
        length, pow_real_permuted_time, pow_real_permuted_mflops,
        pow_real_mflops / pow_real_permuted_mflops);
    fprintf(output_file, " %10d permuted Integer  mod:  %8.2f  %8.2f %8.2f\n",
        length, mod_int_permuted_time, mod_int_permuted_mflops,
        mod_int_mflops  / mod_int_permuted_mflops);

    time_of_sqrt_real_array_perm(rdst, rsrc4, permute, length,
                                  &sqrt_real_permuted_time,
                                  &sqrt_real_permuted_mflops);
    time_of_cos_real_array_perm (rdst, rsrc3, permute, length,
                                  &cos_real_permuted_time,
                                  &cos_real_permuted_mflops);

    fprintf(output_file,
               " %10d permuted Real    sqrt:  %8.2f  %8.2f %8.2f\n",
        length, sqrt_real_permuted_time, sqrt_real_permuted_mflops,
        sqrt_real_mflops / sqrt_real_permuted_mflops);
    fprintf(output_file,
               " %10d permuted Real     cos:  %8.2f  %8.2f %8.2f\n",
        length, cos_real_permuted_time, cos_real_permuted_mflops,
        cos_real_mflops  / cos_real_permuted_mflops);

    time_of_exp_real_array_perm(rdst, rsrc4, permute, length,
                                 &exp_real_permuted_time,
                                 &exp_real_permuted_mflops);
    time_of_log_real_array_perm(rdst, rsrc2, permute, length,
                                 &log_real_permuted_time,
                                 &log_real_permuted_mflops);

    fprintf(output_file,
        " %10d permuted Real     exp:  %8.2f  %8.2f %8.2f\n",
        length, exp_real_permuted_time, exp_real_permuted_mflops,
        exp_real_mflops / exp_real_permuted_mflops);
    fprintf(output_file,
        " %10d permuted Real     log:  %8.2f  %8.2f %8.2f\n",
        length, log_real_permuted_time, log_real_permuted_mflops,
        log_real_mflops / log_real_permuted_mflops);

    time_of_dot_real_array_perm(&dot, rsrc1, rsrc3, permute, length,
                                 &dot_real_permuted_time,
                                 &dot_real_permuted_mflops);

    fprintf(output_file,
        " %10d permuted Real     dot:  %8.2f  %8.2f %8.2f\n",
        length, dot_real_permuted_time, dot_real_permuted_mflops,
        dot_real_mflops / dot_real_permuted_mflops);

    time_of_euc_real_array_perm(rdst, rsrc1, rsrc2, rsrc3, rsrc4,
                                 permute, length,
                                 &euc_real_permuted_time,
                                 &euc_real_permuted_mflops);

    fprintf(output_file,
        " %10d permuted Real     euc:  %8.2f  %8.2f %8.2f\n",
        length, euc_real_permuted_time, euc_real_permuted_mflops,
        euc_real_mflops / euc_real_permuted_mflops);

    time_of_lot8_real_array_perm(&lot, rsrc1, rsrc2, rsrc3, rsrc4,
                                  permute, length,
                                  &lot8_real_permuted_time,
                                  &lot8_real_permuted_mflops);

    fprintf(output_file,
        " %10d permuted Real    lot8:  %8.2f  %8.2f %8.2f\n",
        length, lot8_real_permuted_time, lot8_real_permuted_mflops,
        lot8_real_mflops / lot8_real_permuted_mflops);

    time_of_lot10_real_array_perm(&lot, rsrc1, rsrc2, rsrc3, rsrc4,
                                   permute, length,
                                   &lot10_real_permuted_time,
                                   &lot10_real_permuted_mflops);

    fprintf(output_file,
        " %10d permuted Real   lot10:  %8.2f  %8.2f %8.2f\n",
        length, lot10_real_permuted_time, lot10_real_permuted_mflops,
        lot10_real_mflops / lot10_real_permuted_mflops);

    time_of_lot12_real_array_perm(&lot, rsrc1, rsrc2, rsrc3, rsrc4,
                                   permute, length,
                                   &lot12_real_permuted_time,
                                   &lot12_real_permuted_mflops);

    fprintf(output_file,
        " %10d permuted Real   lot12:  %8.2f  %8.2f %8.2f\n",
        length, lot12_real_permuted_time, lot12_real_permuted_mflops,
        lot12_real_mflops / lot12_real_permuted_mflops);

    time_of_lot16_real_array_perm(&lot, rsrc1, rsrc2, rsrc3, rsrc4,
                                   permute, length,
                                   &lot16_real_permuted_time,
                                   &lot16_real_permuted_mflops);

    fprintf(output_file,
        " %10d permuted Real   lot16:  %8.2f  %8.2f %8.2f\n",
        length, lot16_real_permuted_time, lot16_real_permuted_mflops,
        lot16_real_mflops / lot16_real_permuted_mflops);

    time_of_lot20_real_array_perm(&lot, rsrc1, rsrc2, rsrc3, rsrc4,
                                   permute, length,
                                   &lot20_real_permuted_time,
                                   &lot20_real_permuted_mflops);

    fprintf(output_file,
        " %10d permuted Real   lot20:  %8.2f  %8.2f %8.2f\n",
        length, lot20_real_permuted_time, lot20_real_permuted_mflops,
        lot20_real_mflops / lot20_real_permuted_mflops);

    time_of_lot24_real_array_perm(&lot, rsrc1, rsrc2, rsrc3, rsrc4,
                                   permute, length,
                                   &lot24_real_permuted_time,
                                   &lot24_real_permuted_mflops);

    fprintf(output_file,
        " %10d permuted Real   lot24:  %8.2f  %8.2f %8.2f\n",
        length, lot24_real_permuted_time, lot24_real_permuted_mflops,
        lot24_real_mflops / lot24_real_permuted_mflops);

    time_of_i2r_array_perm(rdst, isrc3, permute, length,
                            &i2r_permuted_time,
                            &i2r_permuted_mflops);
    time_of_r2i_array_perm(idst, rsrc1, permute, length,
                            &r2i_permuted_time,
                            &r2i_permuted_mflops);

    fprintf(output_file,
        " %10d permuted          i2r:  %8.2f  %8.2f %8.2f\n",
        length, i2r_permuted_time, i2r_permuted_mflops,
        i2r_mflops / i2r_permuted_mflops);
    fprintf(output_file,
        " %10d permuted          r2i:  %8.2f  %8.2f %8.2f\n",
        length, r2i_permuted_time, r2i_permuted_mflops,
        r2i_mflops / r2i_permuted_mflops);

    fclose(output_file);

    free(rdst);
    free(rsrc1);
    free(rsrc2);
    free(rsrc3);
    free(rsrc4);

    free(idst);
    free(isrc1);
    free(isrc2);
    free(isrc3);
    free(isrc4);
    free(permute);

    return program_success_code;

} /* main */

void set_permutation (int permute[], int length)
{ /* set_permutation */
    int random_value;
    int index, random_index, temporary_index;

    for (index = 0; index < length; index++) {

       permute[index] = index;

    } /* for index */

    srand48(mysecond());

    for (index = 0; index < length; index++) {

        random_index = floor(drand48() * length);
        assert(random_index < length);

        if (random_index != index) {
            temporary_index       = permute[index];
            permute[index]        = permute[random_index];
            permute[random_index] = temporary_index;
        } /* if (random_index != index) */

   } /* for index */

} /* set_permutation */

void check_permutation (int permute[], int length)
{ /* check_permutation */
    int index;

    for (index = 0; index < length; index++) {
        assert(permute[index] < length);
   } /* for index */

} /* check_permutation */

void initialize_source_array (
         float rsrc1[], float rsrc2[], float rsrc3[], float rsrc4[],
         int   isrc1[], int   isrc2[], int   isrc3[], int   isrc4[],
         int   length)
{ /* initialize_source_array */
    int index;

    for (index = 0; index < length; index++) {
        rsrc1[index] =  index * 0.375;
        rsrc2[index] =  index * 1.375;
        rsrc3[index] =  ((float)index - 1.0) / length;
        rsrc4[index] =  rsrc3[index] * rsrc3[index];
        isrc1[index] =  index * 2 + 1;
        isrc2[index] =  index * 3 + 2;
        isrc3[index] = -index * 7;
        isrc4[index] =  index * 11 + 7;
   } /* for index */

} /* initialize_source_array */

void time_of_add_real_array (
         float dst[], float src1[], float src2[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_add_real_array */
    double  start_time, end_time;
    int index;

    start_time = mysecond();

    for (index = 0; index < length; index++) {
        dst[index] = src1[index] + src2[index];
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        1.0 * (((float)length) / *time_of_operation) * 1.0E-06;

} /* time_of_add_real_array */

void time_of_add_int_array (
         int dst[], int src1[], int src2[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_add_int_array */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    for (index = 0; index < length; index++) {
        dst[index] = src1[index] + src2[index];
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        1.0 * (((float)length) / *time_of_operation) * 1.0E-06;

} /* time_of_add_int_array */

void time_of_sum_real_array (
                float dst[], float src[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_sum_real_array */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    *dst = 0.0;
    for (index = 0; index < length; index++) {
        *dst = *dst + src[index];
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        1.0 * (((float)length) / *time_of_operation) * 1.0E-06;

} /* time_of_sum_real_array */

void time_of_sum_int_array (
         int *dst, int src[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_sum_int_array */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    *dst = 0;

    for (index = 0; index < length; index++) {
        *dst = *dst + src[index];
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        1.0 * (((float)length) / *time_of_operation) * 1.0E-06;

} /* time_of_sum_int_array */

void time_of_sub_real_array (
                float dst[], float src1[], float src2[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_sub_real_array */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    for (index = 0; index < length; index++) {
        dst[index] = src1[index] - src2[index];
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        1.0 * (((float)length) / *time_of_operation) * 1.0E-06;

} /* time_of_sub_real_array */

void time_of_sub_int_array (
                int dst[], int src1[], int src2[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_sub_int_array */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    for (index = 0; index < length; index++) {
        dst[index] = src1[index] - src2[index];
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        1.0 * (((float)length) / *time_of_operation) * 1.0E-06;

} /* time_of_sub_int_array */

void time_of_mul_real_array (
                float dst[], float src1[], float src2[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_mul_real_array */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    for (index = 0; index < length; index++) {
        dst[index] = src1[index] * src2[index];
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        1.0 * (((float)length) / *time_of_operation) * 1.0E-06;
} /* time_of_mul_real_array */

void time_of_mul_int_array (
                int dst[], int src1[], int src2[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_mul_int_array */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    for (index = 0; index < length; index++) {
        dst[index] = src1[index] * src2[index];
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        1.0 * (((float)length) / *time_of_operation) * 1.0E-06;

} /* time_of_mul_int_array */

void time_of_mad_real_array (
                float dst[], float src1[], float src2[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_mad_real_array */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    for (index = 0; index < length; index++) {
        dst[index] = src1[index] + 5.0 * src2[index];
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        2.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_mad_real_array */

void time_of_mad_int_array (
                int dst[], int src1[], int src2[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_mad_int_array */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    for (index = 0; index < length; index++) {
        dst[index] = src1[index] + 5 * src2[index];
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        2.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_mad_int_array */

void time_of_mam_real_array (
                float dst[],
         float src1[], float src2[], float src3[], float src4[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_mam_real_array */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    for (index = 0; index < length; index++) {
        dst[index] =
            src1[index] * src2[index] + src3[index] * src4[index];
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        3.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_mam_real_array */

void time_of_mam_int_array (
                int dst[],
         int src1[], int src2[], int src3[], int src4[], int length,
                double  *time_of_operation, double *mflops_of_operation)
{ /* time_of_mam_int_array */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    for (index = 0; index < length; index++) {
        dst[index] =
            src1[index] * src2[index] + src3[index] * src4[index];
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        3.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_mam_int_array */

void time_of_div_real_array (
                float dst[], float src1[], float src2[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_div_real_array */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    for (index = 0; index < length; index++) {
        dst[index] = src1[index] / src2[index];
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        1.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_div_real_array */

void time_of_div_int_array (
                int dst[], int src1[], int src2[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_div_int_array */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    for (index = 0; index < length; index++) {
        dst[index] = src1[index] / src2[index];
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        1.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_div_int_array */

void time_of_pow_real_array (
                float dst[], float src1[], float src2[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_pow_real_array */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    if (sizeof(float) == 4) {
        for (index = 0; index < length; index++) {
            dst[index] = powf(src1[index], src2[index]);
        } /* for index */
    } /* if (sizeof(float) == 4) */
    else {
        for (index = 0; index < length; index++) {
            dst[index] = pow(src1[index], src2[index]);
        } /* for index */
    } /* if (sizeof(float) == 4)...else */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        1.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_pow_real_array */

void time_of_mod_int_array (
                int dst[], int src1[], int src2[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_mod_int_array */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    for (index = 0; index < length; index++) {
       dst[index] = src1[index] % src2[index];
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        1.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_mod_int_array */

void time_of_sqrt_real_array (
                float dst[], float src[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_sqrt_real_array */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    if (sizeof(float) == 4) {
        for (index = 0; index < length; index++) {
            dst[index] = sqrtf(src[index]);
        } /* for index */
    } /* if (sizeof(float) == 4) */
    else {
        for (index = 0; index < length; index++) {
            dst[index] = sqrt(src[index]);
        } /* for index */
    } /* if (sizeof(float) == 4)...else */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        1.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_sqrt_real_array */

void time_of_cos_real_array (
                float dst[], float src[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_cos_real_array */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    if (sizeof(float) == 4) {
        for (index = 0; index < length; index++) {
            dst[index] = cosf(src[index]);
        } /* for index */
    } /* if (sizeof(float) == 4) */
    else {
        for (index = 0; index < length; index++) {
            dst[index] = cos(src[index]);
        } /* for index */
    } /* if (sizeof(float) == 4)...else */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        1.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_cos_real_array */

void time_of_exp_real_array (
                float dst[], float src[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_exp_real_array */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    if (sizeof(float) == 4) {
       for (index = 0; index < length; index++) {
           dst[index] = expf(src[index]);
       } /* for index */
    } /* if (sizeof(float) == 4) */
    else {
       for (index = 0; index < length; index++) {
           dst[index] = exp(src[index]);
       } /* for index */
    } /* if (sizeof(float) == 4)...else */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        1.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_exp_real_array */

void time_of_log_real_array (
                float dst[], float src[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_log_real_array */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    if (sizeof(float) == 4) {
        for (index = 0; index < length; index++) {
            dst[index] = logf(src[index]);
        } /* for index */
    } /* if (sizeof(float) == 4) */
    else {
        for (index = 0; index < length; index++) {
            dst[index] = log(src[index]);
        } /* for index */
    } /* if (sizeof(float) == 4)...else */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        ((float)length / *time_of_operation) * 1.0E-06;

} /* time_of_log_real_array */

void time_of_euc_real_array (
         float dst[], float src1[], float src2[], float src3[], float src4[],
         int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_euc_real_array */
    double start_time, end_time;
    float diff12, diff34;
    int index;

    start_time = mysecond();

    if (sizeof(float) == 4) {
        for (index = 0; index < length; index++) {
            diff12 = src1[index] - src2[index];
            diff34 = src3[index] - src4[index];
            dst[index] = sqrtf(diff12 * diff12 + diff34 * diff34);
        } /* for index */
    } /* if (sizeof(float) == 4) */
    else {
        for (index = 0; index < length; index++) {
            diff12 = src1[index] - src2[index];
            diff34 = src3[index] - src4[index];
            dst[index] = sqrt(diff12 * diff12 + diff34 * diff34);
        } /* for index */
    } /* if (sizeof(float) == 4)...else */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        6.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_euc_real_array */

void time_of_dot_real_array (
         float *dot, float src1[], float src2[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_dot_real_array */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    *dot = 0.0;

    for (index = 0; index < length; index++) {
        *dot = *dot + src1[index] * src2[index];
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        2.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_dot_real_array */

void time_of_lot8_real_array (
         float *lot, float src1[], float src2[], float src3[], float src4[],
         int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_lot8_real_array */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    *lot = 0.0;
    for (index = 0; index < length; index++) {
        *lot +=
            (src1[index] * src2[index]) +
            (src3[index] * src4[index]) +
            (src1[index] + src2[index]) *
            (src3[index] + src4[index]);
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        8.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_lot8_real_array */

void time_of_lot10_real_array(
         float *lot, float src1[], float src2[], float src3[], float src4[],
         int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_lot10_real_array */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    *lot = 0.0;
    for (index = 0; index < length; index++) {
        *lot +=
            (src1[index] * src2[index]) +
            (src3[index] * src4[index]) +
            (src1[index] + src2[index]) *
            (src3[index] + src4[index]) *
            (src1[index] - src2[index]);
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        10.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_lot10_real_array */

void time_of_lot12_real_array (
         float *lot, float src1[], float src2[], float src3[], float src4[],
         int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_lot12_real_array */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    *lot = 0.0;
    for (index = 0; index < length; index++) {
        *lot +=
            (src1[index] * src2[index]) +
            (src3[index] * src4[index]) +
            (src1[index] + src2[index]) *
            (src3[index] + src4[index]) *
            (src1[index] - src2[index]) *
            (src3[index] - src4[index]);
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        12.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_lot12_real_array */

void time_of_lot16_real_array (
         float *lot, float src1[], float src2[], float src3[], float src4[],
         int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_lot16_real_array */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    *lot = 0.0;
    for (index = 0; index < length; index++) {
        *lot +=
            (src1[index] * src2[index]) +
            (src3[index] * src4[index]) +
            (src1[index] + src2[index]) *
            (src3[index] + src4[index]) *
            (src1[index] - src2[index]) *
            (src3[index] - src4[index]) *
            ((src1[index] - src3[index]) +
             (src2[index] - src4[index]));
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        16.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_lot16_real_array */

void time_of_lot20_real_array (
         float *lot, float src1[], float src2[], float src3[], float src4[],
         int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_lot20_real_array */

    double start_time, end_time;
    int index;

    start_time = mysecond();

    *lot = 0.0;
    for (index = 0; index < length; index++) {
        *lot +=
           (src1[index] * src2[index]) +
           (src3[index] * src4[index]) +
           (src1[index] + src2[index]) *
           (src3[index] + src4[index]) *
           (src1[index] - src2[index]) *
           (src3[index] - src4[index]) *
           ((src1[index] - src3[index]) +
            (src2[index] - src4[index])) *
           ((src1[index] + src3[index]) -
            (src2[index] + src4[index]));
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        20.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_lot20_real_array */

void time_of_lot24_real_array (
         float *lot, float src1[], float src2[], float src3[], float src4[],
         int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_lot24_real_array */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    *lot = 0.0;
    for (index = 0; index < length; index++) {
        *lot +=
            (src1[index] * src2[index]) +
            (src3[index] * src4[index]) +
            (src1[index] + src2[index]) *
            (src3[index] + src4[index]) *
            (src1[index] - src2[index]) *
            (src3[index] - src4[index]) *
            ((src1[index] - src3[index]) +
             (src2[index] - src4[index])) *
            ((src1[index] + src3[index]) -
             (src2[index] + src4[index])) +
            (src1[index] * src3[index]) +
            (src2[index] * src4[index]);
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        24.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_lot24_real_array */

void time_of_i2r_array (
         float rdst[], int isrc[],
         int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_i2r_array */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    for (index = 0; index < length; index++) {
         rdst[index] = (float)isrc[index];
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        1.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_i2r_array */

void time_of_r2i_array (
         int idst[], float rsrc[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_r2i_array */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    for (index = 0; index < length; index++) {
      idst[index] = rsrc[index];
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        1.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_r2i_array */

void time_of_add_real_array_perm (
         float dst[], float src1[], float src2[],
         int permute[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_add_real_array_perm */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    for (index = 0; index < length; index++) {
        dst[permute[index]] = src1[permute[index]] + src2[permute[index]];
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        1.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_add_real_array_perm */

void time_of_add_int_array_perm (
         int dst[], int src1[], int src2[],
         int permute[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_add_int_array_perm */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    for (index = 0; index < length; index++) {
      dst[permute[index]] = src1[permute[index]] + src2[permute[index]];
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        1.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_add_int_array_perm */

void time_of_sum_real_array_perm (
         float *dst, float src[],
         int permute[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_sum_real_array_perm */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    *dst = 0.0;
    for (index = 0; index < length; index++) {
        *dst = *dst + src[permute[index]];
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        1.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_sum_real_array_perm */

void time_of_sum_int_array_perm (
         int dst[], int src[],
         int permute[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_sum_int_array_perm */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    *dst = 0;
    for (index = 0; index < length; index++) {
        *dst = *dst + src[permute[index]];
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        1.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /*  time_of_sum_int_array_perm */

void time_of_sub_real_array_perm (
         float dst[], float src1[], float src2[],
         int permute[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_sub_real_array_perm */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    for (index = 0; index < length; index++) {
        dst[permute[index]] = src1[permute[index]] - src2[permute[index]];
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        1.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_sub_real_array_perm */

void time_of_sub_int_array_perm (
         int dst[], int src1[], int src2[],
         int permute[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_sub_int_array_perm */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    for (index = 0; index < length; index++) {
        dst[permute[index]] = src1[permute[index]] - src2[permute[index]];
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        1.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_sub_int_array_perm */

void time_of_mul_real_array_perm (
         float dst[], float src1[], float src2[],
         int permute[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_mul_real_array_perm */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    for (index = 0; index < length; index++) {
        dst[permute[index]] = src1[permute[index]] * src2[permute[index]];
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        1.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_mul_real_array_perm */

void time_of_mul_int_array_perm (
         int dst[], int src1[], int src2[],
         int permute[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_mul_int_array_perm */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    for (index = 0; index < length; index++) {
        dst[permute[index]] = src1[permute[index]] * src2[permute[index]];
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        1.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_mul_int_array_perm */

void time_of_mad_real_array_perm (
         float dst[], float src1[], float src2[],
         int permute[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_mad_real_array_perm */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    for (index = 0; index < length; index++) {
        dst[permute[index]] =
            src1[permute[index]] + 5.0 * src2[permute[index]];
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        2.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_mad_real_array_perm */

void time_of_mad_int_array_perm (
         int dst[], int src1[], int src2[],
         int permute[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_mad_int_array_perm */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    for (index = 0; index < length; index++) {
        dst[permute[index]] =
            src1[permute[index]] + 5 * src2[permute[index]];
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        2.0 * (float)(length / *time_of_operation) * 1.0E-06;
} /* time_of_mad_int_array_perm */

void time_of_mam_real_array_perm (
         float dst[], float src1[], float src2[], float src3[], float src4[],
         int permute[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_mam_real_array_perm */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    for (index = 0; index < length; index++) {
        dst[permute[index]] =
            src1[permute[index]] * src2[permute[index]] +
            src3[permute[index]] * src4[permute[index]];
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        3.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_mam_real_array_perm */

void time_of_mam_int_array_perm (
         int dst[], int src1[], int src2[], int src3[], int src4[],
         int permute[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_mam_int_array_perm */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    for (index = 0; index < length; index++) {
        dst[permute[index]] =
            src1[permute[index]] * src2[permute[index]] +
            src3[permute[index]] * src4[permute[index]];
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        3.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /*  time_of_mam_int_array_perm */

void time_of_div_real_array_perm (
         float dst[], float src1[], float src2[],
         int permute[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_mam_int_array_perm */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    for (index = 0; index < length; index++) {
        dst[permute[index]] = src1[permute[index]] / src2[permute[index]];
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        1.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_div_real_array_perm */

void time_of_div_int_array_perm (
         int dst[], int src1[], int src2[],
         int permute[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_div_int_array_perm */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    for (index = 0; index < length; index++) {
        dst[permute[index]] = src1[permute[index]] / src2[permute[index]];
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        1.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /*  time_of_div_int_array_perm */

void time_of_pow_real_array_perm (
         float dst[], float src1[], float src2[],
         int permute[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_pow_real_array_perm */
    double start_time, end_time;
    int index;

    start_time = mysecond();


//          DO index = 1, length
//            dst(permute(index)) =  &
//  &            src1(permute(index)) ** src2(permute(index))
//      DO !! index

    if (sizeof(float) == 4) {
        for (index = 0; index < length; index++) {
            dst[permute[index]] =
                powf(src1[permute[index]], src2[permute[index]]);
        } /* for index */
    } /* if (sizeof(float) == 4) */
    else {
        for (index = 0; index < length; index++) {
            dst[permute[index]] =
                pow(src1[permute[index]], src2[permute[index]]);
        } /* for index */
    } /* if (sizeof(float) == 4)...else */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        1.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_pow_real_array_perm */

void time_of_mod_int_array_perm (
         int dst[], int src1[], int src2[],
         int permute[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_mod_int_array_perm */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    for (index = 0; index < length; index++) {
        dst[permute[index]] = src1[permute[index]] % src2[permute[index]];
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        1.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_mod_int_array_perm */

void time_of_sqrt_real_array_perm (
         float dst[], float src[],
         int permute[], int length,
         double * time_of_operation, double *mflops_of_operation)
{ /* time_of_sqrt_real_array_perm */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    if (sizeof(float) == 4) {
        for (index = 0; index < length; index++) {
            dst[permute[index]] = sqrtf(src[permute[index]]);
        } /* for index */
    } /* if (sizeof(float) == 4) */
    else {
        for (index = 0; index < length; index++) {
            dst[permute[index]] = sqrt(src[permute[index]]);
        } /* for index */
    } /* if (sizeof(float) == 4)...else */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        1.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_sqrt_real_array_perm */

void time_of_cos_real_array_perm (
         float dst[], float src[],
         int permute[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_cos_real_array_perm */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    if (sizeof(float) == 4) {
        for (index = 0; index < length; index++) {
            dst[permute[index]] = cosf(src[permute[index]]);
        } /* for index */
    } /* if (sizeof(float) == 4) */
    else {
        for (index = 0; index < length; index++) {
            dst[permute[index]] = cos(src[permute[index]]);
        } /* for index */
    } /* if (sizeof(float) == 4)...else */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        1.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_cos_real_array_perm */

void time_of_exp_real_array_perm (
         float dst[], float src[],
         int permute[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_exp_real_array_perm */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    if (sizeof(float) == 4) {
        for (index = 0; index < length; index++) {
            dst[permute[index]] = expf(src[permute[index]]);
        } /* for index */
    } /* if (sizeof(float) == 4) */
    else {
        for (index = 0; index < length; index++) {
            dst[permute[index]] = exp(src[permute[index]]);
        } /* for index */
    } /* if (sizeof(float) == 4)...else */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        1.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_exp_real_array_perm */

void time_of_log_real_array_perm (
         float dst[], float src[],
         int permute[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_log_real_array_perm */
    double start_time, end_time;
    int index;

    start_time = mysecond();


   //       DO index = 1, length
   //           dst(permute(index)) = LOG(src(permute(index)))
   //      DO !! index

    if (sizeof(float) == 4) {
       for (index = 0; index < length; index++) {
           dst[permute[index]] = logf(src[permute[index]]);
       } /* for index */
    } /* if (sizeof(float) == 4) */
    else {
       for (index = 0; index < length; index++) {
           dst[permute[index]] = log(src[permute[index]]);
       } /* for index */
    } /* if (sizeof(float) == 4)...else */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        1.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_log_real_array_perm */

void time_of_euc_real_array_perm (
         float dst[], float src1[], float src2[], float src3[], float src4[],
         int permute[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_euc_real_array_perm */
    double diff12, diff34;
    double start_time, end_time;
    int index;

    start_time = mysecond();

    if (sizeof(float) == 4) {
        for (index = 0; index < length; index++) {
            diff12 = src1[permute[index]] - src2[permute[index]];
            diff34 = src3[permute[index]] - src4[permute[index]];
            dst[permute[index]] = sqrtf(diff12 * diff12 + diff34 * diff34);
        } /* for index */
    } /* if (sizeof(float) == 4) */
    else {
        for (index = 0; index < length; index++) {
            diff12 = src1[permute[index]] - src2[permute[index]];
            diff34 = src3[permute[index]] - src4[permute[index]];
            dst[permute[index]] = sqrt(diff12 * diff12 + diff34 * diff34);
        } /* for index */
    } /* if (sizeof(float) == 4)...else */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        6.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_euc_real_array_perm */

void time_of_dot_real_array_perm (
         float *dot, float src1[], float src2[],
         int permute[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_dot_real_array_perm */
    double diff12, diff34;
    double start_time, end_time;
    int index;

    start_time = mysecond();

    *dot = 0.0;
    for (index =0; index < length; index++) {
        *dot = *dot + src1[permute[index]] * src2[permute[index]];
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        2.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_dot_real_array_perm */

void time_of_lot8_real_array_perm (
         float *lot, float src1[], float src2[], float src3[], float src4[],
         int permute[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_lot8_real_array_perm */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    *lot = 0.0;
    for (index = 0 ; index < length; index++) {
        *lot +=
            (src1[permute[index]] * src2[permute[index]]) +
            (src3[permute[index]] * src4[permute[index]]) +
            (src1[permute[index]] + src2[permute[index]]) *
            (src3[permute[index]] + src4[permute[index]]);
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        8.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /*  time_of_lot8_real_array_perm */

void time_of_lot10_real_array_perm (
         float *lot, float src1[], float src2[], float src3[], float src4[],
         int permute[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_lot10_real_array_perm */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    *lot = 0.0;
    for (index = 0; index < length; index++) {
        *lot +=
            (src1[permute[index]] * src2[permute[index]]) +
            (src3[permute[index]] * src4[permute[index]]) +
            (src1[permute[index]] + src2[permute[index]]) *
            (src3[permute[index]] + src4[permute[index]]) *
            (src1[permute[index]] - src2[permute[index]]);
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        10.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_lot10_real_array_perm */

// stopped here on Feb 18, 2009

void time_of_lot12_real_array_perm (
         float *lot, float src1[], float src2[], float src3[], float src4[],
         int permute[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_lot12_real_array_perm */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    *lot = 0.0;
    for (index = 0; index < length; index++) {
        *lot +=
            (src1[permute[index]] * src2[permute[index]]) +
            (src3[permute[index]] * src4[permute[index]]) +
            (src1[permute[index]] + src2[permute[index]]) *
            (src3[permute[index]] + src4[permute[index]]) *
            (src1[permute[index]] - src2[permute[index]]) *
            (src3[permute[index]] - src4[permute[index]]);
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        12.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_lot12_real_array_perm */

void time_of_lot16_real_array_perm (
         float *lot, float src1[], float src2[], float src3[], float src4[],
         int permute[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_lot16_real_array_perm */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    *lot = 0.0;
    for (index = 0; index < length; index++) {
        *lot +=
            (src1[permute[index]] * src2[permute[index]]) +
            (src3[permute[index]] * src4[permute[index]]) +
            (src1[permute[index]] + src2[permute[index]]) *
            (src3[permute[index]] + src4[permute[index]]) *
            (src1[permute[index]] - src2[permute[index]]) *
            (src3[permute[index]] - src4[permute[index]]) *
            ((src1[permute[index]] - src3[permute[index]]) +
             (src2[permute[index]] - src4[permute[index]]));
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        16.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_lot16_real_array_perm */

void time_of_lot20_real_array_perm (
         float *lot, float src1[], float src2[], float src3[], float src4[],
         int permute[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_lot20_real_array_perm */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    *lot = 0.0;
    for (index = 0; index < length; index++) {
        *lot +=
            (src1[permute[index]] * src2[permute[index]]) +
            (src3[permute[index]] * src4[permute[index]]) +
            (src1[permute[index]] + src2[permute[index]]) *
            (src3[permute[index]] + src4[permute[index]]) *
            (src1[permute[index]] - src2[permute[index]]) *
            (src3[permute[index]] - src4[permute[index]]) *
            ((src1[permute[index]] - src3[permute[index]]) +
             (src2[permute[index]] - src4[permute[index]])) *
            ((src1[permute[index]] + src3[permute[index]]) -
             (src2[permute[index]] + src4[permute[index]]));
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        20.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_lot20_real_array_perm */

void time_of_lot24_real_array_perm (
         float *lot, float src1[], float src2[], float src3[], float src4[],
         int permute[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_lot24_real_array_perm */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    *lot = 0.0;
    for (index = 0; index < length; index++) {
        *lot +=
            (src1[permute[index]] * src2[permute[index]]) +
            (src3[permute[index]] * src4[permute[index]]) +
            (src1[permute[index]] + src2[permute[index]]) *
            (src3[permute[index]] + src4[permute[index]]) *
            (src1[permute[index]] - src2[permute[index]]) *
            (src3[permute[index]] - src4[permute[index]]) *
            ((src1[permute[index]] - src3[permute[index]]) +
             (src2[permute[index]] - src4[permute[index]])) *
            ((src1[permute[index]] + src3[permute[index]]) -
             (src2[permute[index]] + src4[permute[index]])) +
            (src1[permute[index]] * src3[permute[index]]) +
            (src2[permute[index]] * src4[permute[index]]);
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        24.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_lot24_real_array_perm */

void time_of_i2r_array_perm (
         float rdst[], int isrc[],
         int permute[], int length,
         double *time_of_operation,  double *mflops_of_operation)
{ /* time_of_i2r_array_perm */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    for (index = 0; index < length; index++) {
      rdst[permute[index]] = isrc[permute[index]];
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        1.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_i2r_array_perm */

void time_of_r2i_array_perm (
         int idst[], float rsrc[],
         int permute[], int length,
         double *time_of_operation, double *mflops_of_operation)
{ /* time_of_r2i_array_perm */
    double start_time, end_time;
    int index;

    start_time = mysecond();

    for (index = 0; index < length; index++) {
       idst[permute[index]] = rsrc[permute[index]];
    } /* for index */

    end_time = mysecond();

    *time_of_operation = end_time - start_time;
    *mflops_of_operation =
        1.0 * ((float)length / *time_of_operation) * 1.0E-06;
} /* time_of_r2i_array_perm */


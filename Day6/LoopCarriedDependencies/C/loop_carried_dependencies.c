/*
 * loop_carried_dependencies.c
 */

#include <stdio.h>
#include <stdlib.h>
#include "timings.h"

#define DFLT_LENGTH  (1024*1024*16)

static int length = DFLT_LENGTH;
static void get_arguments(int argc, char* argv[]);
static void initialize_sources(double* src1, double* src2, int length);
static void loop_dependency(
                double* dst, double* src1, double* src2, int length);
static void loop_antidependency(
                double* dst, double* src1, double* src2, int length);
extern double second(void);

int main (int argc, char* argv[])
{ /* main */
    double* dst  = (double*)NULL;
    double* src1 = (double*)NULL;
    double* src2 = (double*)NULL;
    double time_at_start, time_at_end;
    double time_dst_src1_src2;
    double time_dst_src1_src1;
    double time_dst_dst_src2;
    double time_dst_src1_dst;
    double time_dst_dst_dst;

    get_arguments(argc, argv);
    dst  = (double*)malloc(sizeof(double) * length);
    src1 = (double*)malloc(sizeof(double) * length);
    src2 = (double*)malloc(sizeof(double) * length);
    time_at_start = second();
    initialize_sources(src1, src2, length);
    time_at_start = second();
    loop_dependency(dst, src1, src2, length);
    time_at_end   = second();
    time_dst_src1_src2 = time_at_end - time_at_start;
    time_at_start = second();
    loop_dependency(dst, src1, src1, length);
    time_at_end   = second();
    time_dst_src1_src1 = time_at_end - time_at_start;
    time_at_start = second();
    loop_dependency(dst, dst,  src2, length);
    time_at_end   = second();
    time_dst_dst_src2  = time_at_end - time_at_start;
    time_at_start = second();
    loop_dependency(dst, src1, dst,  length);
    time_at_end   = second();
    time_dst_src1_dst  = time_at_end - time_at_start;
    time_at_start = second();
    loop_dependency(dst, dst,  dst,  length);
    time_at_end   = second();
    time_dst_dst_dst   = time_at_end - time_at_start;
    printf("CASE                    SEC     MFLOPs   SPEEDUP\n");
    printf("d[i]=s1[i-1]+s2[i]   %7.2f   %7.2f   %7.2f\n",
        time_dst_src1_src2, 1.0E-06 * (double)length / time_dst_src1_src2,
        1.0);
    printf("d[i]=s1[i-1]+s1[i]   %7.2f   %7.2f   %7.2f\n",
        time_dst_src1_src1, 1.0E-06 * (double)length / time_dst_src1_src1,
        (double)time_dst_src1_src2 / time_dst_src1_src1);
    printf("d[i]=d[i-1]+s2[i]    %7.2f   %7.2f   %7.2f\n",
        time_dst_dst_src2,  1.0E-06 * (double)length / time_dst_dst_src2,
        (double)time_dst_src1_src2 / time_dst_dst_src2);
    printf("d[i]=s1[i-1]+d[i]    %7.2f   %7.2f   %7.2f\n",
        time_dst_src1_dst,  1.0E-06 * (double)length / time_dst_src1_dst,
        (double)time_dst_src1_src2 / time_dst_src1_dst);
    printf("d[i]=d[i-1]+d[i]     %7.2f   %7.2f   %7.2f\n",
        time_dst_dst_dst,   1.0E-06 * (double)length / time_dst_dst_dst,
        (double)time_dst_src1_src2 / time_dst_dst_dst);
    time_at_start = second();
    loop_antidependency(dst, src1, src2, length);
    time_at_end   = second();
    time_dst_src1_src2 = time_at_end - time_at_start;
    time_at_start = second();
    loop_antidependency(dst, src1, src1, length);
    time_at_end   = second();
    time_dst_src1_src1 = time_at_end - time_at_start;
    time_at_start = second();
    loop_antidependency(dst, dst,  src2, length);
    time_at_end   = second();
    time_dst_dst_src2  = time_at_end - time_at_start;
    time_at_start = second();
    loop_antidependency(dst, src1, dst,  length);
    time_at_end   = second();
    time_dst_src1_dst  = time_at_end - time_at_start;
    time_at_start = second();
    loop_antidependency(dst, dst,  dst,  length);
    time_at_end   = second();
    time_dst_dst_dst   = time_at_end - time_at_start;
    printf("d[i]=s1[i+1]+s2[i]   %7.2f   %7.2f   %7.2f\n",
        time_dst_src1_src2, 1.0E-06 * (double)length / time_dst_src1_src2,
        1.0);
    printf("d[i]=s1[i+1]+s1[i]   %7.2f   %7.2f   %7.2f\n",
        time_dst_src1_src1, 1.0E-06 * (double)length / time_dst_src1_src1,
        (double)time_dst_src1_src2 / time_dst_src1_src1);
    printf("d[i]=d[i+1]+s2[i]    %7.2f   %7.2f   %7.2f\n",
        time_dst_dst_src2,  1.0E-06 * (double)length / time_dst_dst_src2,
        (double)time_dst_src1_src2 / time_dst_dst_src2);
    printf("d[i]=s1[i+1]+d[i]    %7.2f   %7.2f   %7.2f\n",
        time_dst_src1_dst,  1.0E-06 * (double)length / time_dst_src1_dst,
        (double)time_dst_src1_src2 / time_dst_src1_dst);
    printf("d[i]=d[i+1]+d[i]     %7.2f   %7.2f   %7.2f\n",
        time_dst_dst_dst,   1.0E-06 * (double)length / time_dst_dst_dst,
        (double)time_dst_src1_src2 / time_dst_dst_dst);
} /* main */

void get_arguments (int argc, char* argv[])
{ /* get_arguments */
    if (argc <= 1) {
        length = DFLT_LENGTH;
    } /* if (argc <= 1) */
    else {
        sscanf(argv[1], "%d", &length);
        if (length < 1) {
            fprintf(stderr, "ERROR: can't work on arrays of length %d\n",
                length);
            exit(-1);
        } /* if (length < 1) */
    } /* if (argc <= 1) */
} /* get_arguments */

void initialize_sources (double* src1, double* src2, int length)
{ /* initialize_sources */
    int index;

    for (index = 0; index < length; index++) {
        src1[index] = 0.01 * index;
        src2[index] = 1.0 / index;
    } /* for index */
} /* initialize_sources */

void loop_dependency (double* dst, double* src1, double* src2, int length)
{ /* loop_dependency */
    int index;

    if ((dst == src1) && (dst == src2)) {
        for (index = 1; index < length; index++) {
            dst[index] = dst[index-1]  + dst[index];
        } /* for index */
    } /* if ((dst == src1) && (dst == src2)) */
    else if (dst == src1) {
        for (index = 1; index < length; index++) {
            dst[index] = dst[index-1]  + src2[index];
        } /* for index */
    } /* if (dst == src1) */
    else if (dst == src2) {
        for (index = 1; index < length; index++) {
            dst[index] = src1[index-1] + dst[index];
        } /* for index */
    } /* if (dst == src2) */
    else if (src1 == src2) {
        for (index = 1; index < length; index++) {
            dst[index] = src1[index-1] + src1[index];
        } /* for index */
    } /* if (dst == src2) */
    else {
        for (index = 1; index < length; index++) {
            dst[index] = src1[index-1] + src2[index];
        } /* for index */
    } /* if (dst == src2)...else */
} /* loop_dependency */

void loop_antidependency (double* dst, double* src1, double* src2, int length)
{ /* loop_antidependency */
    int lengthm1;
    int index;

    lengthm1 = length - 1;
    if ((dst == src1) && (dst == src2)) {
        for (index = 0; index < lengthm1; index++) {
            dst[index] = dst[index+1]  + dst[index];
        } /* for index */
    } /* if ((dst == src1) && (dst == src2)) */
    else if (dst == src1) {
        for (index = 0; index < lengthm1; index++) {
            dst[index] = dst[index+1]  + src2[index];
        } /* for index */
    } /* if (dst == src1) */
    else if (dst == src2) {
        for (index = 0; index < lengthm1; index++) {
            dst[index] = src1[index+1] + dst[index];
        } /* for index */
    } /* if (dst == src2) */
    else if (src1 == src2) {
        for (index = 0; index < lengthm1; index++) {
            dst[index] = src1[index+1] + src1[index];
        } /* for index */
    } /* if (dst == src2) */
    else {
        for (index = 0; index < lengthm1; index++) {
            dst[index] = src1[index+1] + src2[index];
        } /* for index */
    } /* if (dst == src2)...else */
} /* loop_antidependency */


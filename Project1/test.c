#include <stdio.h>

int main(void)
{
    int i = 0;
    while ((i++ < 10) && (i++ % 2 == 0))
        printf("In %d\n", i);
    printf("After %d\n", i);
}
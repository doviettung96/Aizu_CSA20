#include <stdio.h>
#define BUFLEN 12
int main(void)
{
    char buf[BUFLEN];
    while (fgets(buf, BUFLEN, stdin) != NULL)
        printf("%s", buf);
    return 0;
}
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include "fetch_line.h"
static char *trim_line(char *s);
char *fetch_line(char *buf, int buflen, FILE *stream, int *lineno);

static char *trim_line(char *s)
{
    char *start = s;
    while (isspace(*start))
    {
        start++;
    }

    char *ptr = start;
    char *end = start;
    while ((*ptr != '\0') && (*ptr != '#'))
    {
        ptr++;
        if (!isspace(*ptr))
        {
            end = ptr;
        }
    }
    ptr++;

    *end = '\0';
    return start;
}

char *fetch_line(char *buf, int buflen, FILE *stream, int *lineno)
{
    char *s;
    if (fgets(buf, buflen, stream) == NULL)
        return NULL;
    ++*lineno;
    if (buf[strlen(buf) - 1] != '\n')
    {
        fprintf(stderr, "*** reading error: input line %d too long for %s's buf[%d]\n", *lineno, __func__, buflen);
        exit(EXIT_FAILURE);
    }
    s = trim_line(buf);
    if (*s != '\0')
        return s;
    else
        return fetch_line(buf, buflen, stream, lineno);
}

Question 1:
Output:
```
viet@envfs1:~/Projects/HPC/Project1$ ./fetch-line-demo < fgets_demo.c
trimmed line   3: int main(void)

trimmed line   4: {

trimmed line   5: char buf[BUFLEN];

*** reading error: input line 6 too long for fetch_line's buf[40]
```

Question 2:
Output:
```
viet@envfs1:~/Projects/HPC/Project1$ ./fetch-line-demo < data2.txt 
trimmed line   1: 35xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

trimmed line   2: 36xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

*** reading error: input line 3 too long for fetch_line's buf[40]
```

Conclusion: all of the lines with length 39 (+1 of the end character) will be skip. In this case is 37xxx
Other than that, nothing happend since the lines are trimmed already.

Question 3:
Output:
```
viet@envfs1:~/Projects/HPC/Project1$ ./fetch-line-demo < data3.txt 
trimmed line   1: 35xxxxxxxxxxxxxxxxx xxxxxxxxxxxxxxxxxx

trimmed line   2: 35xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 

trimmed line   3: 35xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

trimmed line   5: 35xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
trimmed line   6: 35xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
trimmed line   7: 36xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

*** reading error: input line 8 too long for fetch_line's buf[40]
```

Conclusion: the program works as expected.
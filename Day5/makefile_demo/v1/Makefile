shallow_water_1d: shallow_water_1d.o boundary_conditions.o \
  initial_conditions.o r8mat_write.o r8vec_linspace_new.o r8vec_write.o \
  timestamp.o
	gcc shallow_water_1d.o boundary_conditions.o \
  initial_conditions.o r8mat_write.o r8vec_linspace_new.o r8vec_write.o \
  timestamp.o -o shallow_water_1d

shallow_water_1d.o: shallow_water_1d.c boundary_conditions.h \
  initial_conditions.h r8mat_write.h r8vec_linspace_new.h r8vec_write.h \
  timestamp.h
	gcc -c -Wall -pedantic -std=c99 -O2 shallow_water_1d.c

boundary_conditions.o: boundary_conditions.c boundary_conditions.h
	gcc -c -Wall -pedantic -std=c99 -O2 boundary_conditions.c

initial_conditions.o: initial_conditions.c initial_conditions.h
	gcc -c -Wall -pedantic -std=c99 -O2 initial_conditions.c

r8mat_write.o: r8mat_write.c r8mat_write.h
	gcc -c -Wall -pedantic -std=c99 -O2 r8mat_write.c

r8vec_linspace_new.o: r8vec_linspace_new.c r8vec_linspace_new.h
	gcc -c -Wall -pedantic -std=c99 -O2 r8vec_linspace_new.c

r8vec_write.o: r8vec_write.c r8vec_write.h
	gcc -c -Wall -pedantic -std=c99 -O2 r8vec_write.c

timestamp.o: timestamp.c timestamp.h
	gcc -c -Wall -pedantic -std=c99 -O2 timestamp.c



OFILES = shallow_water_1d.o boundary_conditions.o \
  initial_conditions.o r8mat_write.o r8vec_linspace_new.o r8vec_write.o \
  timestamp.o
TARGET = shallow_water_1d
CC = gcc
CFLAGS = -Wall -pedantic -std=c99 -O2
$(TARGET): $(OFILES)
	$(CC) $(OFILES) -o $@

shallow_water_1d.o: shallow_water_1d.c boundary_conditions.h \
  initial_conditions.h r8mat_write.h r8vec_linspace_new.h r8vec_write.h \
  timestamp.h

boundary_conditions.o: boundary_conditions.c boundary_conditions.h


initial_conditions.o: initial_conditions.c initial_conditions.h


r8mat_write.o: r8mat_write.c r8mat_write.h


r8vec_linspace_new.o: r8vec_linspace_new.c r8vec_linspace_new.h


r8vec_write.o: r8vec_write.c r8vec_write.h


timestamp.o: timestamp.c timestamp.h



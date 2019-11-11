CC=mpicc
FLAGS= -lm -lgsl -fpic -O3 -Wall -Wextra -fopenmp
CTYPES= -shared -o libcMobCal.dll -DTEST_CTYPES=1

cmobcal: cMobCal.c
	$(CC) cMobCal.c $(FLAGS) -o cMobCal.exe

test: cMobCal.c
	$(CC) -DTEST=1 cMobCal.c $(FLAGS) -pg  -o cMobCal.exe

libcmobcal: cMobCal.c
	$(CC)  cMobCal.c $(FLAGS) $(CTYPES)

ptest: cMobCal.c
	mpicc -DTEST=1 cMobCal.c $(FLAGS) -pg  -o cMobCal.exe


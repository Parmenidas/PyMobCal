CC=gcc
FLAGS= -lm -lgsl -pthread -fpic -O3 -Wall -Wextra
CTYPES= -shared -o libcMobCal.dll -DTEST_CTYPES=1

cmobcal: cMobCal.c
	$(CC) cMobCal.c $(FLAGS) -o cMobCal.exe

test: cMobCal.c
	$(CC) cMobCal.c $(FLAGS) -pg -DTEST=1 -o cMobCal.exe

libcmobcal: cMobCal.c
	$(CC)  cMobCal.c $(FLAGS) $(CTYPES)

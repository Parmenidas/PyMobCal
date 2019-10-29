CC=gcc
FLAGS= -lm -lgsl -fpic -O3 -Wall -Wextra
CTYPES= -shared -o libcMobCal.dll

cmobcal: cMobCal.c
	$(CC) cMobCal.c $(FLAGS) -o cMobCal.exe

libcmobcal:
	$(CC)  cMobCal.c $(FLAGS) $(CTYPES)
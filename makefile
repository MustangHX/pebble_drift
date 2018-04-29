
all:

	#module load gcc/5.2.0 

#gcc -E main.c -o main.o
#	gcc -E drift.c -o drift.o
#	gcc -c main.c
#	gcc -c drift.c
#	gcc -o main.o drift.o -o pebble
#	gcc -mcmodel=medium *.c -o pebble -lstdc++
	gcc *.c -o pebble -lstdc++
archive :
	@echo "Creating src_peb.tar"
	@tar cf src_peb.tar *.c
	@tar rf src_peb.tar *.h
	@tar rf src_peb.tar makefile

clean:
	rm -rf *.o pebble

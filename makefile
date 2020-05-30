
all:

# module load gcc/5.2.0

# gcc -E main.c -o main.o
#	gcc -E drift.c -o drift.o
#	gcc -c main.c
#	gcc -c drift.c
#	gcc -o main.o drift.o -o pebble
	gcc -g *.c -o pebble
#	gcc -mcmodel=medium ./*.c -o pebble -lm -lpthread #-lstdc++
#	gcc main.c disk.c drift.c global_var.c group.c Init.c interaction.c opaczhu.c spline.c growth.c output_check.c -o pebble -lstdc++
# gcc -mcmodel=medium main.c disk.c drift.c global_var.c group.c Init.c interaction.c opaczhu.c spline.c -o pebble -lstdc++
archive :
	@echo "Creating src_peb.tar"
	@tar cf src_peb.tar *.c
	@tar rf src_peb.tar *.h
	@tar rf src_peb.tar makefile

clean:
	rm -rf *.o pebble

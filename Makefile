main: main.o ./src/matrix.lib
	gcc -lm -l main.o matrix.lib -o main
main.o: main.c
	gcc main.c -o main.o
matrix.lib: ./src/matrix.c
	gcc -g ./src/matrix.c ./src/complex.c -lm -o ./src/matrix.lib
clean:
	rm ./src/matrix.lib ./main.o ./main
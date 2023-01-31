# compiler used
CC=g++
# flags used
FLAGS=-O3
# further libraries required
LIB=-lm -larmadillo -lstdc++fs

# Makefile
central: *.cpp *.hpp
	$(CC) $(FLAGS) *.cpp *.hpp -o central $(LIB)

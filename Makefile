# compiler used
CC=g++
# flags used
FLAGS=-O3
# further libraries required
LIB=-lm -larmadillo -lstdc++fs
# directory where the include files are located
DIR_INC=include
# directory where the source is located
DIR_SRC=src

# Makefile
central: $(DIR_SRC)/*.cpp $(DIR_INC)/*.hpp
	$(CC) $(FLAGS) $(DIR_SRC)/*.cpp -o central $(LIB)

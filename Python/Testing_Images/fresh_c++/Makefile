
FLAGS = -std=c++11 -stdlib=libc++ -Wall -Wextra -pedantic -DPARALLEL=1
CC_OPT = -O3 -ffast-math -fno-common
PATH_JPEG = /usr/local/Cellar/jpeg/9a/
PATH_BLAS = /usr/local/Cellar/openblas/0.2.16/
LIBJPEG = -L$(PATH_JPEG)lib/ -I$(PATH_JPEG)include/ -ljpeg -DWITH_JPEG=1
LIBOPENBLAS = -L$(PATH_BLAS)lib/ -I$(PATH_BLAS)include/ -lopenblas
ACCELERATE = -framework Accelerate
BIN = fresh
SRC = src/main.cpp

compile_accelerate_nojpeg:
	g++ $(FLAGS) $(CC_OPT) -o $(BIN) $(ACCELERATE) $(SRC)

compile_accelerate_withjpeg:
	g++ $(FLAGS) $(CC_OPT) -o $(BIN) $(ACCELERATE) $(LIBJPEG) $(SRC)

compile_openblas_nojpeg:
	g++ $(FLAGS) $(CC_OPT) -o $(BIN) $(LIBOPENBLAS) $(SRC)

compile_openblas_withjpeg:
	g++ $(FLAGS) $(CC_OPT) -o $(BIN) $(LIBOPENBLAS) $(LIBJPEG) $(SRC)

test: compile_accelerate_nojpeg
	./$(BIN) data/cameraman256.bmp cameraman256_out.bmp 2 1 1 1 0

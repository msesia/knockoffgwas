#compiler
CXX=gcc

CURRENT_DIR=.

# Compiler flags
FLAGS = -O3 -Wno-deprecated-declarations

TARGET = libsqlite3.a

SRC_DIR = ${CURRENT_DIR}/sqlite3
OBJ_DIR = ${CURRENT_DIR}/build
LIB_DIR = ${CURRENT_DIR}/../../lib
SRC_FILES = $(wildcard $(SRC_DIR)/*.c)
OBJ_FILES = $(patsubst $(SRC_DIR)/%.c,$(OBJ_DIR)/%.o,$(SRC_FILES))

all : ${LIB_DIR}/${TARGET}

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	${CXX} $(FLAGS) -o $@ -c $<

${LIB_DIR}/${TARGET} : $(OBJ_FILES)
	ar rvs ${LIB_DIR}/${TARGET} ${OBJ_DIR}/*.o

clean :
	rm -f ${OBJ_DIR}/*
	rm -f ${LIB_DIR}/${TARGET}

# Create output directory
DIRS = ${OBJ_DIR} ${LIB_DIR}
$(shell mkdir -p $(DIRS))

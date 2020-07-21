#compiler
CXX=g++

CURRENT_DIR=.

# Compiler flags
FLAGS = -O3 \
-I${CURRENT_DIR}/genfile/include \
-I${CURRENT_DIR}/include \
-I${CURRENT_DIR}/../sqlite3

TARGET = libdb.a

SRC_DIR = ${CURRENT_DIR}/src
OBJ_DIR = ${CURRENT_DIR}/build
LIB_DIR = ${CURRENT_DIR}/../../lib
SRC_FILES = $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))

all : ${LIB_DIR}/${TARGET}

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	${CXX} $(FLAGS) -o $@ -c $<

${LIB_DIR}/${TARGET} : $(OBJ_FILES)
	ar rvs ${LIB_DIR}/${TARGET} ${OBJ_DIR}/*.o

clean :
	rm -f ${OBJ_DIR}/*
	rm -f ${LIB_DIR}/${TARGET}

# Create output directory
DIRS = ${OBJ_DIR} ${LIB_DIR}
$(shell mkdir -p $(DIRS))

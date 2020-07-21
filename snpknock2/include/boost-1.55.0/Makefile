#compiler
CXX=g++

# Compiler flags
FLAGS = -std=c++11 -Wno-unused-local-typedefs -Wno-deprecated-declarations -I.

TARGET = libboost.a

CURRENT_DIR=.

SRC_DIR = ${CURRENT_DIR}/src
OBJ_DIR = ${CURRENT_DIR}/build
LIB_DIR = ${CURRENT_DIR}/../../lib

SRC_SYSTEM = $(wildcard $(SRC_DIR)/system/src/*.cpp)
OBJ_SYSTEM = $(patsubst $(SRC_DIR)/system/src/%.cpp,$(OBJ_DIR)/system/%.o,$(SRC_SYSTEM))

SRC_THREAD = $(wildcard $(SRC_DIR)/thread/src/*.cpp)
OBJ_THREAD = $(patsubst $(SRC_DIR)/thread/src/%.cpp,$(OBJ_DIR)/thread/%.o,$(SRC_THREAD))

SRC_FILESYSTEM = $(wildcard $(SRC_DIR)/filesystem/src/*.cpp)
OBJ_FILESYSTEM = $(patsubst $(SRC_DIR)/filesystem/src/%.cpp,$(OBJ_DIR)/filesystem/%.o,$(SRC_FILESYSTEM))

SRC_DATETIME = $(wildcard $(SRC_DIR)/datetime/src/posix_time/*.cpp)
OBJ_DATETIME = $(patsubst $(SRC_DIR)/datetime/src/posix_time/%.cpp,$(OBJ_DIR)/datetime/%.o,$(SRC_DATETIME))

SRC_TIMER = $(wildcard $(SRC_DIR)/timer/src/*.cpp)
OBJ_TIMER = $(patsubst $(SRC_DIR)/timer/src/%.cpp,$(OBJ_DIR)/timer/%.o,$(SRC_TIMER))

SRC_CHRONO = $(wildcard $(SRC_DIR)/chrono/src/*.cpp)
OBJ_CHRONO = $(patsubst $(SRC_DIR)/chrono/src/%.cpp,$(OBJ_DIR)/chrono/%.o,$(SRC_CHRONO))

SRC_FILES = ${SRC_SYSTEM} ${SRC_FILESYSTEM} ${SRC_THREAD} ${SRC_DATETIME} ${SRC_TIMER} ${SRC_CHRONO}
OBJ_FILES = ${OBJ_SYSTEM} ${OBJ_FILESYSTEM} ${OBJ_THREAD} ${OBJ_DATETIME} ${OBJ_TIMER} ${OBJ_CHRONO}

all : ${LIB_DIR}/${TARGET}

info:
	@echo "Info"
	@echo "Sources: "${SRC_FILES}
	@echo "Objects: "${OBJ_FILES}

$(OBJ_DIR)/system/%.o: $(SRC_DIR)/system/src/%.cpp
	${CXX} $(FLAGS) -o $@ -c $<

$(OBJ_DIR)/filesystem/%.o: $(SRC_DIR)/filesystem/src/%.cpp
	${CXX} $(FLAGS) -o $@ -c $<

$(OBJ_DIR)/thread/%.o: $(SRC_DIR)/thread/src/%.cpp
	${CXX} $(FLAGS) -o $@ -c $<

$(OBJ_DIR)/datetime/%.o: $(SRC_DIR)/datetime/src/%.cpp
	${CXX} $(FLAGS) -o $@ -c $<

$(OBJ_DIR)/timer/%.o: $(SRC_DIR)/timer/src/%.cpp
	${CXX} $(FLAGS) -o $@ -c $<

$(OBJ_DIR)/chrono/%.o: $(SRC_DIR)/chrono/src/%.cpp
	${CXX} $(FLAGS) -o $@ -c $<

${LIB_DIR}/${TARGET} : $(OBJ_FILES)
	ar rvs ${LIB_DIR}/${TARGET} $(OBJ_FILES)

clean :
	@rm -rf ${OBJ_DIR}/*
	@rm -f ${LIB_DIR}/${TARGET}

# Create output directory
DIRS = ${OBJ_DIR} ${OBJ_DIR}/system ${OBJ_DIR}/filesystem ${OBJ_DIR}/thread ${OBJ_DIR}/datetime ${OBJ_DIR}/timer  ${OBJ_DIR}/chrono ${LIB_DIR}
$(shell mkdir -p $(DIRS))

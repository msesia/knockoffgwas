#compiler
CXX=g++

# Executable file
EXEC = bin/snpknock2

# Static libraries
LIBS_NAMES = alglib.a libboost.a libsqlite3.a libdb.a libzstd.a libbgen.a
STATIC_LIBS = $(addprefix lib/, ${LIBS_NAMES})

# Compiler flags
# -DWITHGPERFTOOLS -g
# -O3
FLAGS = -Ofast -std=c++11 \
-Iinclude/bgen/genfile/include \
-Iinclude/zstd-1.1.0/src \
-Iinclude/alglib-3.15/src \
-Llib

#-L/home/groups/candes/Software/usr/lib \
# -lprofiler

# Linker flags
LINKS = -lpthread -lboost_iostreams -lboost_thread -lboost_filesystem -lboost_system -lz -lzstd -ldb -lsqlite3

# List of objects
OBJECTS = build/utils.o build/bgen_utils.o build/bgen_reader.o build/haps_reader.o build/generator_master.o build/main.o build/chaplotype.o build/filter_reader.o build/arguments.o build/utils_graph.o build/cross_validation.o build/haplotypes.o build/knockoffs.o build/dataset.o build/kinship_computer.o build/kinship.o build/covariates.o build/windows.o build/metadata.o build/ibd.o build/family_bp.o build/family_knockoffs.o build/plink_interface.o ${STATIC_LIBS}

all : bin/snpknock2

bin/snpknock2 : ${OBJECTS}
	${CXX} ${FLAGS} -o bin/snpknock2 ${OBJECTS} ${LINKS}

build/main.o : src/main.cpp src/utils.h src/arguments.h src/generator_master.h src/dataset.h src/kinship.h src/metadata.h src/covariates.h
	${CXX} ${FLAGS} -o build/main.o -c src/main.cpp

build/generator_master.o : src/generator_master.h src/generator_master.cpp src/utils.h src/haplotypes.h src/knockoffs.h src/filter_reader.h src/metadata.h
	${CXX} ${FLAGS} -o build/generator_master.o -c src/generator_master.cpp

build/kinship.o : src/kinship.h src/kinship.cpp src/kinship_computer.h src/utils.h src/metadata.h
	${CXX} ${FLAGS} -o build/kinship.o -c src/kinship.cpp

build/kinship_computer.o : src/kinship_computer.h src/kinship_computer.cpp src/dataset.h src/metadata.h src/covariates.h src/windows.h
	${CXX} ${FLAGS} -o build/kinship_computer.o -c src/kinship_computer.cpp

build/dataset.o : src/dataset.h src/dataset.cpp src/filter_reader.h src/haplotypes.h src/metadata.h
	${CXX} ${FLAGS} -o build/dataset.o -c src/dataset.cpp

build/cross_validation.o : src/cross_validation.h src/cross_validation.cpp src/utils.h
	${CXX} ${FLAGS} -o build/cross_validation.o -c src/cross_validation.cpp

build/knockoffs.o : src/knockoffs.h src/knockoffs.cpp src/filter_reader.h src/haplotypes.h src/utils.h src/plink_interface.h src/cross_validation.h src/metadata.h src/ibd.h src/family_bp.h src/family_knockoffs.h
	${CXX} ${FLAGS} -o build/knockoffs.o -c src/knockoffs.cpp

build/family_bp.o: src/family_bp.h src/family_bp.cpp src/ibd.h
	${CXX} ${FLAGS} -o build/family_bp.o -c src/family_bp.cpp

build/family_knockoffs.o: src/family_knockoffs.h src/family_knockoffs.cpp src/ibd.h
	${CXX} ${FLAGS} -o build/family_knockoffs.o -c src/family_knockoffs.cpp

build/haplotypes.o : src/haplotypes.h src/haplotypes.cpp src/filter_reader.h src/bgen_reader.h src/haps_reader.h src/chaplotype.h src/metadata.h
	${CXX} ${FLAGS} -o build/haplotypes.o -c src/haplotypes.cpp

build/chaplotype.o : src/chaplotype.h src/chaplotype.cpp src/utils.h
	${CXX} ${FLAGS} -o build/chaplotype.o -c src/chaplotype.cpp

build/filter_reader.o : src/filter_reader.h src/filter_reader.cpp src/utils.h
	${CXX} ${FLAGS} -o build/filter_reader.o -c src/filter_reader.cpp

build/arguments.o : src/arguments.h src/arguments.cpp src/utils.h src/metadata.h
	${CXX} ${FLAGS} -o build/arguments.o -c src/arguments.cpp

build/metadata.o : src/metadata.h src/metadata.cpp src/utils.h src/utils_graph.h src/windows.h
	${CXX} ${FLAGS} -o build/metadata.o -c src/metadata.cpp

build/windows.o : src/windows.h src/windows.cpp src/utils.h src/utils_graph.h
	${CXX} ${FLAGS} -o build/windows.o -c src/windows.cpp

build/ibd.o : src/ibd.h src/ibd.cpp src/utils.h src/utils_graph.h
	${CXX} ${FLAGS} -o build/ibd.o -c src/ibd.cpp

build/covariates.o : src/covariates.h src/covariates.cpp src/utils.h src/metadata.h
	${CXX} ${FLAGS} -o build/covariates.o -c src/covariates.cpp

build/plink_interface.o : src/plink_interface.h src/plink_interface.cpp src/utils.h src/metadata.h src/chaplotype.h
	${CXX} ${FLAGS} -o build/plink_interface.o -c src/plink_interface.cpp

build/utils.o : src/utils.h src/utils.cpp
	${CXX} ${FLAGS} -o build/utils.o -c src/utils.cpp

build/utils_graph.o : src/utils_graph.h src/utils_graph.cpp
	${CXX} ${FLAGS} -o build/utils_graph.o -c src/utils_graph.cpp

build/bgen_utils.o : src/bgen_utils.h src/bgen_utils.cpp src/utils.h
	${CXX} ${FLAGS} -o build/bgen_utils.o -c src/bgen_utils.cpp

build/haps_reader.o : src/haps_reader.h src/haps_reader.cpp src/chaplotype.h
	${CXX} ${FLAGS} -o build/haps_reader.o -c src/haps_reader.cpp

build/bgen_reader.o : src/bgen_reader.h src/bgen_reader.cpp src/chaplotype.h src/metadata.h
	${CXX} ${FLAGS} -o build/bgen_reader.o -c src/bgen_reader.cpp


lib/alglib.a :
	cd include/alglib-3.15/; make -f Makefile; cd ../../

lib/libsqlite3.a :
	cd include/sqlite3/; make -f Makefile; cd ../../

lib/libzstd.a :
	cd include/zstd-1.1.0/; make -f Makefile; cd ../../

lib/libdb.a :
	cd include/db/; make -f Makefile; cd ../../

lib/libbgen.a :
	cd include/bgen/; make -f Makefile; cd ../../

lib/libboost.a :
	cd include/boost-1.55.0/; make -f Makefile; cd ../../


clean :
	rm -f build/*
	rm -f lib/*
	rm -f bin/*

# Create build directory
DIRS = build/ bin/
$(shell mkdir -p $(DIRS))

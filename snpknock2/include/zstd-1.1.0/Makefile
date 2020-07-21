# ################################################################
# Copyright (c) 2016-present, Yann Collet, Facebook, Inc.
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree. An additional grant
# of patent rights can be found in the PATENTS file in the same directory.
# ################################################################

CURRENT_DIR=.

SRC_DIR = ${CURRENT_DIR}/src/

# Version numbers
LIBVER_MAJOR_SCRIPT:=`sed -n '/define ZSTD_VERSION_MAJOR/s/.*[[:blank:]]\([0-9][0-9]*\).*/\1/p' < ${SRC_DIR}/zstd.h`
LIBVER_MINOR_SCRIPT:=`sed -n '/define ZSTD_VERSION_MINOR/s/.*[[:blank:]]\([0-9][0-9]*\).*/\1/p' < ${SRC_DIR}/zstd.h`
LIBVER_PATCH_SCRIPT:=`sed -n '/define ZSTD_VERSION_RELEASE/s/.*[[:blank:]]\([0-9][0-9]*\).*/\1/p' < ${SRC_DIR}/zstd.h`
LIBVER_SCRIPT:= $(LIBVER_MAJOR_SCRIPT).$(LIBVER_MINOR_SCRIPT).$(LIBVER_PATCH_SCRIPT)
LIBVER_MAJOR := $(shell echo $(LIBVER_MAJOR_SCRIPT))
LIBVER_MINOR := $(shell echo $(LIBVER_MINOR_SCRIPT))
LIBVER_PATCH := $(shell echo $(LIBVER_PATCH_SCRIPT))
LIBVER := $(shell echo $(LIBVER_SCRIPT))
VERSION?= $(LIBVER)

CPPFLAGS= -I${SRC_DIR} -I${SRC_DIR}/common
CFLAGS ?= -O3
CFLAGS += -Wall -Wextra -Wcast-qual -Wcast-align -Wshadow -Wstrict-aliasing=1 \
          -Wswitch-enum -Wdeclaration-after-statement -Wstrict-prototypes -Wundef
FLAGS   = $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) $(MOREFLAGS)


ZSTD_FILES := ${SRC_DIR}/common/*.c ${SRC_DIR}/compress/*.c ${SRC_DIR}/decompress/*.c ${SRC_DIR}/dictBuilder/*.c

ifeq ($(ZSTD_LEGACY_SUPPORT), 0)
CPPFLAGS  += -DZSTD_LEGACY_SUPPORT=0
else
ZSTD_FILES+= ${SRC_DIR}/legacy/*.c
CPPFLAGS  += -I${SRC_DIR}/legacy -DZSTD_LEGACY_SUPPORT=1
endif


# OS X linker doesn't support -soname, and use different extension
# see : https://developer.apple.com/library/mac/documentation/DeveloperTools/Conceptual/DynamicLibraries/100-Articles/DynamicLibraryDesignGuidelines.html
ifeq ($(shell uname), Darwin)
	SHARED_EXT = dylib
	SHARED_EXT_MAJOR = $(LIBVER_MAJOR).$(SHARED_EXT)
	SHARED_EXT_VER = $(LIBVER).$(SHARED_EXT)
	SONAME_FLAGS = -install_name $(PREFIX)/lib/$@.$(SHARED_EXT_MAJOR) -compatibility_version $(LIBVER_MAJOR) -current_version $(LIBVER)
else
	SONAME_FLAGS = -Wl,-soname=$@.$(SHARED_EXT).$(LIBVER_MAJOR)
	SHARED_EXT = so
	SHARED_EXT_MAJOR = $(SHARED_EXT).$(LIBVER_MAJOR)
	SHARED_EXT_VER = $(SHARED_EXT).$(LIBVER)
endif


.PHONY: default all clean organize

default: clean libzstd organize

all: clean libzstd organize

libzstd: $(ZSTD_FILES)
	@echo compiling static library
	@$(CC) $(FLAGS) -c $^
	@$(AR) rcs $@.a *.o

clean:
	@rm -f *.o *.a ${CURRENT_DIR}/../../lib/libzstd.a
	@echo Cleaning library completed

organize:
	@mkdir -p ${CURRENT_DIR}/../../lib/
	@rm *.o
	@mv *.a ${CURRENT_DIR}/../../lib/

####  Makefile for compilation on Linux  ####

OPT=-O3     # Optimization option by default

CC=gcc
ifeq "$(CC)" "gcc"
    COMPILER=gcc
else ifeq "$(CC)" "clang"
    COMPILER=clang
endif

ARCHITECTURE=_AMD64_
USE_OPT_LEVEL=_GENERIC_
ifeq "$(ARCH)" "x64"
    ARCHITECTURE=_AMD64_
else ifeq "$(ARCH)" "x86"
    ARCHITECTURE=_X86_
else ifeq "$(ARCH)" "ARM"
    ARCHITECTURE=_ARM_
    ARM_SETTING=-lrt
else ifeq "$(ARCH)" "ARM64"
    ARCHITECTURE=_ARM64_
    ARM_SETTING=-lrt
endif

ifeq "$(OPT_LEVEL)" "GENERIC"
    USE_OPT_LEVEL=_GENERIC_
endif

ifeq "$(ARCHITECTURE)" "_AMD64_"
    ifeq "$(USE_OPT_LEVEL)" "_FAST_"
        MULX=-D _MULX_
        ifeq "$(USE_MULX)" "FALSE"
            MULX=
        else
            ADX=-D _ADX_
            ifeq "$(USE_ADX)" "FALSE"
                ADX=
            endif
        endif
    endif
endif

ifeq "$(SET)" "EXTENDED"
    ADDITIONAL_SETTINGS=-fwrapv -fomit-frame-pointer -march=native
endif

AR=ar rcs
RANLIB=ranlib

CFLAGS=$(OPT) $(ADDITIONAL_SETTINGS) -D $(ARCHITECTURE) -D __LINUX__ -D $(USE_OPT_LEVEL) $(MULX) $(ADX)
LDFLAGS=-lm
ifeq "$(USE_OPT_LEVEL)" "_GENERIC_"
    EXTRA_OBJECTS_503=objs503/fp_generic.o
    EXTRA_OBJECTS_751=objs751/fp_generic.o
    EXTRA_OBJECTS_751_uRadix=src/P751_uRadix/P751/generic/fp_generic.o
    EXTRA_OBJECTS_503_uRadix=src/P503_uRadix/P503/generic/fp_generic.o
endif
OBJECTS_503=objs503/P503.o $(EXTRA_OBJECTS_503) objs/random.o objs/fips202.o
OBJECTS_751=objs751/P751.o $(EXTRA_OBJECTS_751) objs/random.o objs/fips202.o
OBJECTS_751_uRadix=src/P751_uRadix/P751/generic/P751.o $(EXTRA_OBJECTS_751_uRadix) objs/random.o objs/fips202.o
OBJECTS_503_uRadix=src/P503_uRadix/P503/generic/P503.o $(EXTRA_OBJECTS_503_uRadix) objs/random.o objs/fips202.o

all: lib503 lib751 tests KATS lib751uRadix lib503uRadix

objs503/%.o: src/P503_Mont/P503/%.c
	@mkdir -p $(@D)
	$(CC) -c $(CFLAGS) $< -o $@

objs751/%.o: src/P751_Mont/P751/%.c
	@mkdir -p $(@D)
	$(CC) -c $(CFLAGS) $< -o $@

# objs751uRadix/%.o: src/P751_uRadix/P751/%.c
# 	@mkdir -p $(@D)
# 	$(CC) -c $(CFLAGS) $< -o $@

# objs503uRadix/%.o: src/P503_uRadix/P503/%.c
# 	@mkdir -p $(@D)
# 	$(CC) -c $(CFLAGS) $< -o $@

ifeq "$(USE_OPT_LEVEL)" "_GENERIC_"		
    objs503/fp_generic.o: src/P503_Mont/P503/generic/fp_generic.c
	    $(CC) -c $(CFLAGS) src/P503_Mont/P503/generic/fp_generic.c -o objs503/fp_generic.o

    objs751/fp_generic.o: src/P751_Mont/P751/generic/fp_generic.c
	    $(CC) -c $(CFLAGS) src/P751_Mont/P751/generic/fp_generic.c -o objs751/fp_generic.o

    # objs751uRadix/fp_generic.o: src/P751_uRadix/P751/generic/fp_generic.c
	#     $(CC) -c $(CFLAGS) src/P751_uRadix/P751/generic/fp_generic.c -o objs751uRadix/fp_generic.o

    # objs503uRadix/fp_generic.o: src/P503_uRadix/P503/generic/fp_generic.c
	#     $(CC) -c $(CFLAGS) src/P503_uRadix/P503/generic/fp_generic.c -o objs503uRadix/fp_generic.o
endif

objs/random.o: src/random/random.c
	@mkdir -p $(@D)
	$(CC) -c $(CFLAGS) src/random/random.c -o objs/random.o

objs/fips202.o: src/sha3/fips202.c
	$(CC) -c $(CFLAGS) src/sha3/fips202.c -o objs/fips202.o

lib503: $(OBJECTS_503)
	rm -rf lib503 sike
	mkdir lib503 sike
	$(AR) lib503/libsike.a $^
	$(RANLIB) lib503/libsike.a

lib751: $(OBJECTS_751)
	rm -rf lib751 sike
	mkdir lib751 sike
	$(AR) lib751/libsike.a $^
	$(RANLIB) lib751/libsike.a

lib751uRadix: $(OBJECTS_751_uRadix)
	rm -rf lib751uRadix sike
	mkdir lib751uRadix sike
	$(AR) lib751uRadix/libsike.a $^
	$(RANLIB) lib751uRadix/libsike.a

lib503uRadix: $(OBJECTS_503_uRadix)
	rm -rf lib503uRadix sike
	mkdir lib503uRadix sike
	$(AR) lib503uRadix/libsike.a $^
	$(RANLIB) lib503uRadix/libsike.a

tests: lib503 lib751 lib751uRadix lib503uRadix
	$(CC) $(CFLAGS) -L./lib503 tests/test_SIKEp503.c tests/test_extras.c -lsike $(LDFLAGS) -o sike/test_KEM_P503_Mont $(ARM_SETTING)
	$(CC) $(CFLAGS) -L./lib751 tests/test_SIKEp751.c tests/test_extras.c -lsike $(LDFLAGS) -o sike/test_KEM_P751_Mont $(ARM_SETTING)
	$(CC) $(CFLAGS) -L./lib751uRadix tests/test_SIKEp751_uRadix.c tests/test_extras.c -lsike $(LDFLAGS) -o sike/test_KEM_P751_uRadix $(ARM_SETTING)
	$(CC) $(CFLAGS) -L./lib503uRadix tests/test_SIKEp503_uRadix.c tests/test_extras.c -lsike $(LDFLAGS) -o sike/test_KEM_P503_uRadix $(ARM_SETTING)

# AES
AES_OBJS=objs/aes.o objs/aes_c.o

objs/%.o: tests/aes/%.c
	@mkdir -p $(@D)
	$(CC) -c $(CFLAGS) $< -o $@

lib503_for_KATs: $(OBJECTS_503) $(AES_OBJS)
	$(AR) lib503/libsike_for_testing.a $^
	$(RANLIB) lib503/libsike_for_testing.a

lib751_for_KATs: $(OBJECTS_751) $(AES_OBJS)
	$(AR) lib751/libsike_for_testing.a $^
	$(RANLIB) lib751/libsike_for_testing.a

lib751uRadix_for_KATs: $(OBJECTS_751_uRadix) $(AES_OBJS)
	$(AR) lib751uRadix/libsike_for_testing.a $^
	$(RANLIB) lib751uRadix/libsike_for_testing.a

lib503uRadix_for_KATs: $(OBJECTS_503_uRadix) $(AES_OBJS)
	$(AR) lib503uRadix/libsike_for_testing.a $^
	$(RANLIB) lib503uRadix/libsike_for_testing.a

KATS: lib503_for_KATs lib751_for_KATs lib751uRadix_for_KATs lib503uRadix_for_KATs
	$(CC) $(CFLAGS) -L./lib503 tests/PQCtestKAT_kem503.c tests/rng/rng.c -lsike_for_testing $(LDFLAGS) -o sike/PQCtestKAT_kem_P503_Mont $(ARM_SETTING)
	$(CC) $(CFLAGS) -L./lib751 tests/PQCtestKAT_kem751.c tests/rng/rng.c -lsike_for_testing $(LDFLAGS) -o sike/PQCtestKAT_kem_P751_Mont $(ARM_SETTING)
	$(CC) $(CFLAGS) -L./lib751uRadix tests/PQCtestKAT_kem751_uRadix.c tests/rng/rng.c -lsike_for_testing $(LDFLAGS) -o sike/PQCtestKAT_kem_P751_uRadix $(ARM_SETTING)
	$(CC) $(CFLAGS) -L./lib503uRadix tests/PQCtestKAT_kem503_uRadix.c tests/rng/rng.c -lsike_for_testing $(LDFLAGS) -o sike/PQCtestKAT_kem_P503_uRadix $(ARM_SETTING)
check: tests

.PHONY: clean

clean:
	rm -rf *.req objs503 objs lib503 sike objs751 lib751 lib751uRadix lib503uRadix

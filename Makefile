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
    EXTRA_OBJECTS_434=objs434/fp_generic.o
    EXTRA_OBJECTS_503=objs503/fp_generic.o
    EXTRA_OBJECTS_610=objs610/fp_generic.o
    EXTRA_OBJECTS_751=objs751/fp_generic.o
    EXTRA_OBJECTS_751_uRadix=objs751uRadix/fp_generic.o
    EXTRA_OBJECTS_503_uRadix=objs503uRadix/fp_generic.o
	EXTRA_OBJECTS_434_uRadix=objs434uRadix/fp_generic.o
	EXTRA_OBJECTS_610_uRadix=objs610uRadix/fp_generic.o
endif
OBJECTS_434=objs434/P434.o $(EXTRA_OBJECTS_434) objs/random.o objs/fips202.o
OBJECTS_503=objs503/P503.o $(EXTRA_OBJECTS_503) objs/random.o objs/fips202.o
OBJECTS_610=objs610/P610.o $(EXTRA_OBJECTS_610) objs/random.o objs/fips202.o
OBJECTS_751=objs751/P751.o $(EXTRA_OBJECTS_751) objs/random.o objs/fips202.o
OBJECTS_751_uRadix=objs751uRadix/P751.o $(EXTRA_OBJECTS_751_uRadix) objs/random.o objs/fips202.o
OBJECTS_503_uRadix=objs503uRadix/P503.o $(EXTRA_OBJECTS_503_uRadix) objs/random.o objs/fips202.o
OBJECTS_434_uRadix=objs434uRadix/P434.o $(EXTRA_OBJECTS_434_uRadix) objs/random.o objs/fips202.o
OBJECTS_610_uRadix=objs610uRadix/P610.o $(EXTRA_OBJECTS_610_uRadix) objs/random.o objs/fips202.o

all: lib434 lib503 lib610 lib751 tests KATS lib751uRadix lib503uRadix lib610uRadix lib434uRadix

objs434/%.o: src/P434_Mont/P434/%.c
	@mkdir -p $(@D)
	$(CC) -c $(CFLAGS) $< -o $@

objs503/%.o: src/P503_Mont/P503/%.c
	@mkdir -p $(@D)
	$(CC) -c $(CFLAGS) $< -o $@
	
objs610/%.o: src/P610_Mont/P610/%.c
	@mkdir -p $(@D)
	$(CC) -c $(CFLAGS) $< -o $@

objs751/%.o: src/P751_Mont/P751/%.c
	@mkdir -p $(@D)
	$(CC) -c $(CFLAGS) $< -o $@

objs751uRadix/%.o: src/P751_uRadix/P751/%.c
	@mkdir -p $(@D)
	$(CC) -c $(CFLAGS) $< -o $@

objs503uRadix/%.o: src/P503_uRadix/P503/%.c
	@mkdir -p $(@D)
	$(CC) -c $(CFLAGS) $< -o $@
	
objs434uRadix/%.o: src/P434_uRadix/P434/%.c
	@mkdir -p $(@D)
	$(CC) -c $(CFLAGS) $< -o $@

objs610uRadix/%.o: src/P610_uRadix/P610/%.c
	@mkdir -p $(@D)
	$(CC) -c $(CFLAGS) $< -o $@

ifeq "$(USE_OPT_LEVEL)" "_GENERIC_"		
    objs434/fp_generic.o: src/P434_Mont/P434/generic/fp_generic.c
	    $(CC) -c $(CFLAGS) src/P434_Mont/P434/generic/fp_generic.c -o objs434/fp_generic.o

    objs503/fp_generic.o: src/P503_Mont/P503/generic/fp_generic.c
	    $(CC) -c $(CFLAGS) src/P503_Mont/P503/generic/fp_generic.c -o objs503/fp_generic.o
		
    objs610/fp_generic.o: src/P610_Mont/P610/generic/fp_generic.c
	    $(CC) -c $(CFLAGS) src/P610_Mont/P610/generic/fp_generic.c -o objs610/fp_generic.o

    objs751/fp_generic.o: src/P751_Mont/P751/generic/fp_generic.c
	    $(CC) -c $(CFLAGS) src/P751_Mont/P751/generic/fp_generic.c -o objs751/fp_generic.o

    objs751uRadix/fp_generic.o: src/P751_uRadix/P751/generic/fp_generic.c
	    $(CC) -c $(CFLAGS) src/P751_uRadix/P751/generic/fp_generic.c -o objs751uRadix/fp_generic.o

    objs503uRadix/fp_generic.o: src/P503_uRadix/P503/generic/fp_generic.c
	    $(CC) -c $(CFLAGS) src/P503_uRadix/P503/generic/fp_generic.c -o objs503uRadix/fp_generic.o
		
    objs434uRadix/fp_generic.o: src/P434_uRadix/P434/generic/fp_generic.c
	    $(CC) -c $(CFLAGS) src/P434_uRadix/P434/generic/fp_generic.c -o objs434uRadix/fp_generic.o

    objs610uRadix/fp_generic.o: src/P610_uRadix/P610/generic/fp_generic.c
	    $(CC) -c $(CFLAGS) src/P610_uRadix/P610/generic/fp_generic.c -o objs610uRadix/fp_generic.o
endif

objs/random.o: src/random/random.c
	@mkdir -p $(@D)
	$(CC) -c $(CFLAGS) src/random/random.c -o objs/random.o

objs/fips202.o: src/sha3/fips202.c
	$(CC) -c $(CFLAGS) src/sha3/fips202.c -o objs/fips202.o

lib434: $(OBJECTS_434)
	rm -rf lib434 sike
	mkdir lib434 sike
	$(AR) lib434/libsike.a $^
	$(RANLIB) lib434/libsike.a
	
lib503: $(OBJECTS_503)
	rm -rf lib503 sike
	mkdir lib503 sike
	$(AR) lib503/libsike.a $^
	$(RANLIB) lib503/libsike.a
	
lib610: $(OBJECTS_610)
	rm -rf lib610 sike
	mkdir lib610 sike
	$(AR) lib610/libsike.a $^
	$(RANLIB) lib610/libsike.a

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
	
lib434uRadix: $(OBJECTS_434_uRadix)
	rm -rf lib434uRadix sike
	mkdir lib434uRadix sike
	$(AR) lib434uRadix/libsike.a $^
	$(RANLIB) lib434uRadix/libsike.a

lib610uRadix: $(OBJECTS_610_uRadix)
	rm -rf lib610uRadix sike
	mkdir lib610uRadix sike
	$(AR) lib610uRadix/libsike.a $^
	$(RANLIB) lib610uRadix/libsike.a

tests: lib434 lib503 lib610 lib751 lib751uRadix lib503uRadix lib434uRadix lib610uRadix
	$(CC) $(CFLAGS) -L./lib434 tests/test_SIKEp434.c tests/test_extras.c -lsike $(LDFLAGS) -o sike/test_KEM_P434_Mont $(ARM_SETTING)
	$(CC) $(CFLAGS) -L./lib503 tests/test_SIKEp503.c tests/test_extras.c -lsike $(LDFLAGS) -o sike/test_KEM_P503_Mont $(ARM_SETTING)
	$(CC) $(CFLAGS) -L./lib610 tests/test_SIKEp610.c tests/test_extras.c -lsike $(LDFLAGS) -o sike/test_KEM_P610_Mont $(ARM_SETTING)
	$(CC) $(CFLAGS) -L./lib751 tests/test_SIKEp751.c tests/test_extras.c -lsike $(LDFLAGS) -o sike/test_KEM_P751_Mont $(ARM_SETTING)
	$(CC) $(CFLAGS) -L./lib751uRadix tests/test_SIKEp751_uRadix.c tests/test_extras.c -lsike $(LDFLAGS) -o sike/test_KEM_P751_uRadix $(ARM_SETTING)
	$(CC) $(CFLAGS) -L./lib503uRadix tests/test_SIKEp503_uRadix.c tests/test_extras.c -lsike $(LDFLAGS) -o sike/test_KEM_P503_uRadix $(ARM_SETTING)
	$(CC) $(CFLAGS) -L./lib434uRadix tests/test_SIKEp434_uRadix.c tests/test_extras.c -lsike $(LDFLAGS) -o sike/test_KEM_P434_uRadix $(ARM_SETTING)
	$(CC) $(CFLAGS) -L./lib610uRadix tests/test_SIKEp610_uRadix.c tests/test_extras.c -lsike $(LDFLAGS) -o sike/test_KEM_P610_uRadix $(ARM_SETTING)

# AES
AES_OBJS=objs/aes.o objs/aes_c.o

objs/%.o: tests/aes/%.c
	@mkdir -p $(@D)
	$(CC) -c $(CFLAGS) $< -o $@
	
lib434_for_KATs: $(OBJECTS_434) $(AES_OBJS)
	$(AR) lib434/libsike_for_testing.a $^
	$(RANLIB) lib434/libsike_for_testing.a

lib503_for_KATs: $(OBJECTS_503) $(AES_OBJS)
	$(AR) lib503/libsike_for_testing.a $^
	$(RANLIB) lib503/libsike_for_testing.a
	
lib610_for_KATs: $(OBJECTS_610) $(AES_OBJS)
	$(AR) lib610/libsike_for_testing.a $^
	$(RANLIB) lib610/libsike_for_testing.a

lib751_for_KATs: $(OBJECTS_751) $(AES_OBJS)
	$(AR) lib751/libsike_for_testing.a $^
	$(RANLIB) lib751/libsike_for_testing.a

lib751uRadix_for_KATs: $(OBJECTS_751_uRadix) $(AES_OBJS)
	$(AR) lib751uRadix/libsike_for_testing.a $^
	$(RANLIB) lib751uRadix/libsike_for_testing.a

lib503uRadix_for_KATs: $(OBJECTS_503_uRadix) $(AES_OBJS)
	$(AR) lib503uRadix/libsike_for_testing.a $^
	$(RANLIB) lib503uRadix/libsike_for_testing.a
	
lib434uRadix_for_KATs: $(OBJECTS_434_uRadix) $(AES_OBJS)
	$(AR) lib434uRadix/libsike_for_testing.a $^
	$(RANLIB) lib434uRadix/libsike_for_testing.a

lib610uRadix_for_KATs: $(OBJECTS_610_uRadix) $(AES_OBJS)
	$(AR) lib610uRadix/libsike_for_testing.a $^
	$(RANLIB) lib610uRadix/libsike_for_testing.a

KATS: lib434_for_KATs lib503_for_KATs lib610_for_KATs lib751_for_KATs lib751uRadix_for_KATs lib503uRadix_for_KATs lib434uRadix_for_KATs lib610uRadix_for_KATs
	$(CC) $(CFLAGS) -L./lib434 tests/PQCtestKAT_kem434.c tests/rng/rng.c -lsike_for_testing $(LDFLAGS) -o sike/PQCtestKAT_kem_P434_Mont $(ARM_SETTING)
	$(CC) $(CFLAGS) -L./lib503 tests/PQCtestKAT_kem503.c tests/rng/rng.c -lsike_for_testing $(LDFLAGS) -o sike/PQCtestKAT_kem_P503_Mont $(ARM_SETTING)
	$(CC) $(CFLAGS) -L./lib610 tests/PQCtestKAT_kem610.c tests/rng/rng.c -lsike_for_testing $(LDFLAGS) -o sike/PQCtestKAT_kem_P610_Mont $(ARM_SETTING)
	$(CC) $(CFLAGS) -L./lib751 tests/PQCtestKAT_kem751.c tests/rng/rng.c -lsike_for_testing $(LDFLAGS) -o sike/PQCtestKAT_kem_P751_Mont $(ARM_SETTING)
	$(CC) $(CFLAGS) -L./lib751uRadix tests/PQCtestKAT_kem751_uRadix.c tests/rng/rng.c -lsike_for_testing $(LDFLAGS) -o sike/PQCtestKAT_kem_P751_uRadix $(ARM_SETTING)
	$(CC) $(CFLAGS) -L./lib503uRadix tests/PQCtestKAT_kem503_uRadix.c tests/rng/rng.c -lsike_for_testing $(LDFLAGS) -o sike/PQCtestKAT_kem_P503_uRadix $(ARM_SETTING)
	$(CC) $(CFLAGS) -L./lib434uRadix tests/PQCtestKAT_kem434_uRadix.c tests/rng/rng.c -lsike_for_testing $(LDFLAGS) -o sike/PQCtestKAT_kem_P434_uRadix $(ARM_SETTING)
	$(CC) $(CFLAGS) -L./lib610uRadix tests/PQCtestKAT_kem610_uRadix.c tests/rng/rng.c -lsike_for_testing $(LDFLAGS) -o sike/PQCtestKAT_kem_P610_uRadix $(ARM_SETTING)
check: tests

.PHONY: clean

clean:
	rm -rf *.req objs503 objs lib503 sike objs751 lib751 lib751uRadix objs751uRadix lib503uRadix objs503uRadix objs434 lib434 objs610 lib610 lib434uRadix objs434uRadix lib610uRadix objs610uRadix

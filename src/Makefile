
# Standard compiler parameters
CFLAGS += -O -I../libmseed -I../evalresp

# Options specific for GCC
GCC = gcc
GCCFLAGS = -O2 -Wall -I../libmseed -I../evalresp

BIN = tracedsp

LDFLAGS = -L../libmseed -L../evalresp/.libs
LDLIBS = -lmseed -lm -levresp

OBJS = $(BIN).o iirfilter.o complex.o getpzfr.o whiten.o convolve.o decimate.o taper.o rotate.o

all: $(BIN)

$(BIN): $(OBJS)
	$(CC) $(CFLAGS) -o ../$@ $(OBJS) $(LDFLAGS) $(LDLIBS)

clean:
	rm -f $(OBJS) ../$(BIN)

cc:
	@$(MAKE) "CC=$(CC)" "CFLAGS=$(CFLAGS)"

gcc:
	@$(MAKE) "CC=$(GCC)" "CFLAGS=$(GCCFLAGS)"

gcc32:
	@$(MAKE) "CC=$(GCC)" "CFLAGS=-m32 $(GCCFLAGS)"

gcc64:
	@$(MAKE) "CC=$(GCC)" "CFLAGS=-m64 $(GCCFLAGS)"

static:
	@$(MAKE) "CC=$(CC)" "CFLAGS=-static $(GCCFLAGS)"

debug:
	$(MAKE) "CFLAGS=-g $(CFLAGS)"

gccdebug:
	$(MAKE) "CC=$(GCC)" "CFLAGS=-g $(GCCFLAGS)"

gcc32debug:
	$(MAKE) "CC=$(GCC)" "CFLAGS=-g -m32 $(GCCFLAGS)"

gcc64debug:
	$(MAKE) "CC=$(GCC)" "CFLAGS=-g -m64 $(GCCFLAGS)"

install:
	@echo
	@echo "No install target, copy the executable(s) yourself"
	@echo

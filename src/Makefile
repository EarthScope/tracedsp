
# Standard compiler parameters
CFLAGS += -I../libmseed -I../evalresp/libsrc

BIN = tracedsp

LDFLAGS = -L../libmseed -L../evalresp/libsrc/evalresp \
          -L ../evalresp/libsrc/evalresp_log \
          -L ../evalresp/libsrc/spline \
          -L ../evalresp/libsrc/mxml
LDLIBS = -lmseed -lm -levalresp -levalresp_log -lspline -lmxmlev

OBJS = $(BIN).o iirfilter.o complex.o getpzfr.o whiten.o convolve.o decimate.o taper.o envelope.o rotate.o

all: $(BIN)

$(BIN): $(OBJS)
	$(CC) $(CFLAGS) -o ../$@ $(OBJS) $(LDFLAGS) $(LDLIBS)

clean:
	rm -f $(OBJS) ../$(BIN)

cc:
	@$(MAKE) "CC=$(CC)" "CFLAGS=$(CFLAGS)"

static:
	@$(MAKE) "CC=$(CC)" "CFLAGS=-static $(GCCFLAGS)"

install:
	@echo
	@echo "No install target, copy the executable(s) to desired location"
	@echo

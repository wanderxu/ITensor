include ../this_dir.mk
include ../options.mk

#Define Flags ----------

TENSOR_HEADERS=$(PREFIX)/itensor/all.h $(PREFIX)/itensor/mps/idmrg.h
CCFLAGS= -I. $(ITENSOR_INCLUDEFLAGS) $(CPPFLAGS) $(OPTIMIZATIONS)
CCGFLAGS= -I. $(ITENSOR_INCLUDEFLAGS) $(DEBUGFLAGS)
LIBFLAGS=-L$(ITENSOR_LIBDIR) $(ITENSOR_LIBFLAGS)
LIBGFLAGS=-L$(ITENSOR_LIBDIR) $(ITENSOR_LIBGFLAGS)

#Rules ------------------

%.o: %.cc $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

.debug_objs/%.o: %.cc $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) -c $(CCGFLAGS) -o $@ $<

#Targets -----------------

build: dmrg_squ

debug: dmrg_squ-g

all: dmrg_squ

dmrg_squ: dmrg_squ.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) dmrg_squ.o -o dmrg_squ $(LIBFLAGS)

dmrg_squ-g: mkdebugdir .debug_objs/dmrg_squ.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCGFLAGS) .debug_objs/dmrg_squ.o -o dmrg_squ-g $(LIBGFLAGS)

mkdebugdir:
	mkdir -p .debug_objs

clean:
	@rm -fr *.o .debug_objs dmrg_squ dmrg_squ-g

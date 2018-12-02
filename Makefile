SDSL_DIR=../sdsl-lite
include $(SDSL_DIR)/Make.helper

PARALLEL_FLAGS=-fopenmp -pthread
LIBS=-L$(LIB_DIR) -lsdsl -ldivsufsort -ldivsufsort64

OTHER_FLAGS=$(PARALLEL_FLAGS)

CXX_FLAGS=$(MY_CXX_FLAGS) $(OTHER_FLAGS) $(MY_CXX_OPT_FLAGS) -I$(INC_DIR) -I ../gbwt/include
LIBOBJS=path.o
SOURCES=$(wildcard *.cpp)
HEADERS=$(wildcard ../gbwt/include/gbwt/*.h)
OBJS=$(SOURCES:.cpp=.o)

LIBRARY=../gbwt/libgbwt.a
PROGRAMS=path

all:$(LIBRARY) $(PROGRAMS)

%.o:%.cpp $(HEADERS)
	$(MY_CXX) $(CXX_FLAGS) -c $<

$(LIBRARY):$(LIBOBJS)
	ar rcs $@ $(LIBOBJS)

path:path.o $(LIBRARY)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)

clean:
	rm -f $(PROGRAMS) $(OBJS)

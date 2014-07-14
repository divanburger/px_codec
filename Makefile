CC=ccache g++

DCFLAGS=-pipe -Wall -std=c++11 -O0 -g
OCFLAGS=-pipe -Wall -std=c++11 -O3 -g
CFLAGS=$(DCFLAGS)
LFLAGS=
LIBS=-lgsl -lgslcblas

BIN=Coder

SRC=$(shell find src -iname *.cpp)
SOURCES=$(SRC:src/%.cpp=%.cpp)
OBJECTS=$(SOURCES:%.cpp=objs/%.o)

DEPS=$(SRC:src/%.cpp=deps/%.d)

all: release

release: main
	
main: $(OBJECTS)
	@echo Making: $(BIN)
	@$(CC) $(CFLAGS) $(LFLAGS) -o $(BIN) $(OBJECTS) $(LIBS)

-include $(DEPS)

objs/%.o: src/%.cpp deps/%.d
	@echo Compiling: $*
	@mkdir -p $(@D)
	@$(CC) $(CFLAGS) -c $< -o $@

depend: $(DEPS)

deps/%.d: src/%.cpp
	@echo Generating dependencies: $*
	@mkdir -p $(@D)
	@$(CC) $(CFLAGS) -MM -MT 'objs/$*.o' $< > deps/$*.d

clean:
	rm -rf deps
	rm -rf objs
	rm -ff $(BIN)

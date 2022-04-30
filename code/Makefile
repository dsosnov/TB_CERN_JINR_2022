CXXFLAGS+=--std=c++17
CXXFLAGS+=$(shell root-config --cflags) -Wall -I. -ggdb3 -fPIC
LDFLAGS+=$(shell root-config --glibs)  -ggdb3 "-Wl,-rpath,./"

LINKD=LinkDef

HEAD=$(filter-out $(LINKD).h, $(wildcard *.h))

CODE=$(wildcard *.C)

OBJS=$(subst ., ., $(CODE:.C=.o))

all: calibration vmm

link.cpp: link.h LinkDef.h
	rootcint -f $@ -c -p $^

libLink.so: link.cpp
	$(CXX) -shared -o $@ $^ $(LDFLAGS) -fPIC

%.o: %.cxx %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@
	@echo ""

calibration: calibration.o libLink.so
	$(CXX) -o $@ $^ $(LDFLAGS) -L. -lLink -pthread -lm -ldl -rdynamic

vmm: vmm.o libLink.so
	$(CXX) -o $@ $^ $(LDFLAGS) -L. -lLink -pthread -lm -ldl -rdynamic

clean:
	rm -rf *.o
	rm -rf link.cpp *_rdict.pcm
	rm -rf vmm calibration
	rm -rf *.so
	@echo $(OBJS)


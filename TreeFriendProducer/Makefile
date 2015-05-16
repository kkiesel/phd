LIBS = $(shell root-config --evelibs)
INCS = -I$(shell root-config --incdir) -I.

WARN = -Wall -Wshadow

OBJ = main.o #TreeFriendProducer.o TreeObjects.o
EXE = TreeFriendProducer



all: $(EXE) Makefile
	@@

%.o:%.cc
	g++ -std=c++11 -o $@ $+ -c -O2 $(INCS) $(WARN)

$(EXE): $(OBJ)
	g++ -std=c++11 -O2 -o $@ $+ $(LIBS) $(WARN)

TreeObjects_h.so: TreeObjects.*
	root -l -b -q -L TreeObjects.h+
	echo Add `pwd` to current LD_LIBRARY_PATH

lib%.so: %.h %_LinkDef.h
	rootcint -f $(subst so,cxx,$@) -c $+
	g++ -O2 -Wall  -shared -fPIC `root-config --ldflags` $(INCS) -o $@ $(subst lib,,$(subst so,cc,$@)) $(subst so,cxx,$@)




clean:
	@rm -f $(OBJ) $(EXE)
	@rm -f *.cxx lib*h # remove rootcints generated source and header files
	@rm -f AutoDict_* *_h.d *_cc.d *_C.d # cint generated files
	@rm -f *pcm
	@rm -f *so


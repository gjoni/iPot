
CXX=g++
CXXFLAGS += -O3 -mtune=native -Wall -pedantic -std=c++0x
CXXLIBS =
INCDIRS =
LIBDIRS =

OBJDIR = obj
SRCDIR = src

DEPS = $(shell find $(SRCDIR) -name '*.h')
SRCS = $(shell find $(SRCDIR) -name '*.cpp')
OBJS = $(patsubst $(SRCDIR)%.cpp, $(OBJDIR)%.o, $(SRCS))

objrrce=obj/rrce.o
objaace=obj/aace.o

all: $(OBJDIR) aace rrce

$(OBJDIR):
	mkdir -p $(OBJDIR)

aace: $(OBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(INCDIRS) $(LIBDIRS) -o aace $(filter-out $(objrrce), $(OBJS)) $(CXXLIBS)

rrce: $(OBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(INCDIRS) $(LIBDIRS) -o rrce $(filter-out $(objaace), $(OBJS)) $(CXXLIBS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCDIRS) -c -o $@ $<

clean:
	rm $(OBJS) aace rrce
	rmdir $(OBJDIR)


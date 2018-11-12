
CXX=g++
CXXFLAGS += -O3 -mtune=native -Wall -pedantic -std=c++0x -Wno-unused-result
CXXLIBS =
INCDIRS =
LIBDIRS =

OBJDIR = obj
SRCDIR = src

DEPS = $(shell find $(SRCDIR) -name '*.h')
SRCS = $(shell find $(SRCDIR) -name '*.cpp')
OBJS = $(patsubst $(SRCDIR)%.cpp, $(OBJDIR)%.o, $(SRCS))

objrrce=$(filter-out obj/aace.o obj/aace18.o obj/aace167.o obj/aace20.o obj/rrce20.o, $(OBJS))
objaace=$(filter-out obj/rrce.o obj/aace18.o obj/aace167.o obj/aace20.o obj/rrce20.o, $(OBJS))
objaace18=$(filter-out obj/rrce.o obj/aace.o obj/aace167.o obj/aace20.o obj/rrce20.o, $(OBJS))
objaace20=$(filter-out obj/rrce.o obj/aace.o obj/aace167.o obj/aace18.o obj/rrce20.o, $(OBJS))
objaace167=$(filter-out obj/rrce.o obj/aace18.o obj/aace.o obj/aace20.o obj/rrce20.o, $(OBJS))
objrrce20=$(filter-out obj/rrce.o obj/aace.o obj/aace167.o obj/aace18.o obj/aace20.o, $(OBJS))

all: $(OBJDIR) aace18 aace20 aace167 rrce20

$(OBJDIR):
	mkdir -p $(OBJDIR)

aace: $(OBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(INCDIRS) $(LIBDIRS) -o aace $(objaace) $(CXXLIBS)

rrce: $(OBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(INCDIRS) $(LIBDIRS) -o rrce $(objrrce) $(CXXLIBS)

rrce20: $(OBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(INCDIRS) $(LIBDIRS) -o rrce20 $(objrrce20) $(CXXLIBS)

aace18: $(OBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(INCDIRS) $(LIBDIRS) -o aace18 $(objaace18) $(CXXLIBS)

aace20: $(OBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(INCDIRS) $(LIBDIRS) -o aace20 $(objaace20) $(CXXLIBS)

aace167: $(OBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(INCDIRS) $(LIBDIRS) -o aace167 $(objaace167) $(CXXLIBS)


$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCDIRS) -c -o $@ $<

clean:
	rm $(OBJS) aace18 aace20 aace167 rrce20
	rmdir $(OBJDIR)
	rm -rf aace || true
	rm -rf rrce || true



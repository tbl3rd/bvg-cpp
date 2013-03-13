PATH := .:$(PATH)

# Need the following tools on PATH.
#
CXX := g++
COMPARE := cmp
ECHO := echo
GREP := grep
MAKE := make
RM := rm
TIME := time
TOUCH := touch
TR := tr
UNCOMPRESS := gzip -c -d
COUNTLINES := wc -l

CXXFLAGS += -O3

PROBABILITY := 20
LARGE := 10000
SMALL := 500

help:
	@$(ECHO) This makefile defines several build targets.
	@$(ECHO) all: Build both the small and large scale executables.
	@$(ECHO) bvg-small: Build the $(SMALL) scale executable.
	@$(ECHO) bvg-large: Build the $(LARGE) scale executable.
	@$(ECHO) small-ok: Validate bvg-small against parent data.
	@$(ECHO) large-ok: Test bvg-large for sanity.
	@$(ECHO) test: Run both the small-ok and large-ok tests.
	@$(ECHO) clean: Remove any generated files.

all: bvg-small bvg-large
CLEAN += bvg-small bvg-large

bvg-small: CXXFLAGS += -DSCALE=$(SMALL)
bvg-small: bvg.cc
	$(CXX) $(CXXFLAGS) -o $@ $?

bvg-large: CXXFLAGS += -DSCALE=$(LARGE)
bvg-large: bvg.cc
	$(CXX) $(CXXFLAGS) -o $@ $?

bitvectors-genes.data: bitvectors-genes.data.gz
	$(UNCOMPRESS) $? > $@
CLEAN += bitvectors-genes.data

bitvectors-genes.data.small: bitvectors-genes.data.small.gz
	$(UNCOMPRESS) $? > $@
CLEAN += bitvectors-genes.data.small

bitvectors-parents.data.small.txt: bitvectors-parents.data.small.txt.gz
	$(UNCOMPRESS) $? > $@
CLEAN += bitvectors-parents.data.small.txt

small-output.txt: bvg-small bitvectors-genes.data.small
	$(TIME) bvg-small $(PROBABILITY) bitvectors-genes.data.small > $@
CLEAN += small-output.txt

small-ok: small-output.txt bitvectors-parents.data.small.txt
	$(COMPARE) $^ && $(TOUCH) small-ok || $(RM) small-ok
CLEAN += small-ok

bitvectors-parents.data: bvg-large bitvectors-genes.data
	$(TIME) bvg-large $(PROBABILITY) bitvectors-genes.data > $@
CLEAN += bitvectors-parents.data

large-ok: bitvectors-parents.data
	$(COUNTLINES) < bitvectors-parents.data | $(TR) -d ' ' > count.tmp
	$(ECHO) $(LARGE) > scale.tmp
	$(COMPARE) count.tmp scale.tmp && $(TOUCH) large-ok
	$(GREP) -e -1 bitvectors-parents.data | wc -l | $(GREP) -e 1
	$(RM) count.tmp scale.tmp
CLEAN += large-ok

test: small-ok large-ok

clean:
	$(RM) -f $(CLEAN)
	$(RM) -rf *.dSYM

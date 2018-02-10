CC=g++
OPT=-Wno-deprecated -O3 -I include -Wall -std=c++03
SRCS=dash.cpp individual.cpp match.cpp pedigree.cpp
OBJS=dash.o individual.o match.o pedigree.o
MAIN=dash_adv

all: clean dash test

dash: $(OBJS)
	$(CC) $(OPT) -o $(MAIN) $(OBJS)
	cp dash_adv ../../bin

$(OBJS): $(SRCS)
	$(CC) $(OPT) -c $*.cpp

test_cluster:
	-@echo -e "---\nRunning Test Case\n---"
	cat test/test.seg | ./dash_adv -fam test/test.fam -win 50 test/generated
	diff -qs test/expected.ped test/generated.ped
	diff -qs test/expected.hcl test/generated.hcl

test: test_cluster

clean:
	-rm -f *.o $(MAIN) ../../bin/$(MAIN) test/generated.*

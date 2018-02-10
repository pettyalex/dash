CC=g++
OPT=-Wno-deprecated -O3 -I include
SRCS=dash.cpp individual.cpp match.cpp pedigree.cpp
OBJS=dash.o individual.o match.o pedigree.o
MAIN=dash

all: clean dash test

dash: $(OBJS)
	$(CC) $(OPT) -o $(MAIN) $(OBJS)

$(OBJS): $(SRCS)
	$(CC) $(OPT) -c $*.cpp

test_cluster:
	-@echo -e "---\nRunning Test Case\n---"
	cat test/test.seg | ./dash -fam test/test.fam -win 50 test/generated
	diff -qs test/expected.ped test/generated.ped
	diff -qs test/expected.hcl test/generated.hcl

test: test_cluster

clean:
	-rm -f *.o $(MAIN) test/generated.*

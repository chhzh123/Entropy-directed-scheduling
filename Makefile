PCC = g++

ALL = main main-multi-r

all: $(ALL)

% : %.cpp
	$(PCC) -std=c++11 $< -o $@

.PHONY: clean
clean:
	-rm -f *.o $(ALL) *.out
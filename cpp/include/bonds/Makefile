.PHONY: all
all: token-gperf-generated.hpp

.PHONY: clean
clean:
		rm -f token-gperf-generated.hpp token.gperf

token-gperf-generated.hpp: token.gperf Makefile
		gperf -L C++ $(<) >$(@)

token.gperf: gperf_input.py
		python3 gperf_input.py

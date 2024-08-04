extra_LDFLAGS = -Wl,--undefined=lookupIdentifier
extra_CXXFLAGS = -std=c++17 -g -O3 -fvisibility=hidden

.PHONY: all
all: token-gperf-generated.hpp

.PHONY: clean
clean:
		rm -f token-gperf-generated.hpp token.gperf

token-gperf-generated.hpp: token.gperf Makefile
		gperf $(<) >$(@)

token.gperf: gperf_input.py
		python3 gperf_input.py

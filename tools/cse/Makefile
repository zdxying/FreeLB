# CSE source-to-source code generator for FreeLB's unrolled-for (.ur.h)
# specializations.
#
# The engine lives in ../../third_party/cse (git submodule of the standalone
# `cse` repository); this directory only wires it into FreeLB and hosts the
# numerical verifiers.
#
#   make          build the csegen driver (builds the engine first)
#   make gen      regenerate generated/lbm/*.ur.h from ../../src/lbm/*.h (@cse marked)
#   make verify   numerically compare generated .ur.h against the pinned
#                 hand-written references in reference/
#   make install  gen + copy generated .ur.h into src/lbm/ (after verify)
#
# All three marked headers (moment, equilibrium, force) are ported.

CSE_ROOT   := ../../third_party/cse
ENGINE_BIN := $(CSE_ROOT)/bin/csegen
CSEGEN     := csegen

LBM_HDR := ../../src/lbm
REFS    := reference
INPUTS  := moment equilibrium force

# Canonical output location (shared with make.mk build system)
GEN_DIR := ../../generated/lbm

# Always let the engine's own make decide what to rebuild, then refresh the
# driver binary (a plain file dependency would go stale after engine updates).
$(CSEGEN): engine
	cp $(ENGINE_BIN) $@

engine:
	$(MAKE) -C $(CSE_ROOT)

all: $(CSEGEN)

gen: $(CSEGEN)
	@mkdir -p $(GEN_DIR)
	@for f in $(INPUTS); do \
		echo "--- generating $$f.ur.h ---"; \
		./$(CSEGEN) $(LBM_HDR)/$$f.h $(GEN_DIR)/$$f.ur.h; \
	done

verify_equilibrium: gen
	python3 verify_equilibrium.py $(REFS)/equilibrium.ur.h $(GEN_DIR)/equilibrium.ur.h

verify_force: gen
	python3 verify_force.py $(REFS)/force.ur.h $(GEN_DIR)/force.ur.h

verify_moment: gen
	python3 verify_moment.py $(REFS)/moment.ur.h $(GEN_DIR)/moment.ur.h

verify: verify_equilibrium verify_force verify_moment
	@echo "=== verification passed ==="

# Snapshot the current src/lbm/*.ur.h as the verification reference. Only
# needed when the hand-written references change.
gen-refs:
	@mkdir -p $(REFS)
	@for f in $(INPUTS); do cp $(LBM_HDR)/$$f.ur.h $(REFS)/$$f.ur.h; done

install: verify
	@for f in $(INPUTS); do \
		cp $(GEN_DIR)/$$f.ur.h $(LBM_HDR)/$$f.ur.h; \
	done

clean:
	rm -f $(CSEGEN)
	rm -rf out

.PHONY: all engine gen gen-refs verify verify_equilibrium verify_force verify_moment install clean
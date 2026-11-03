# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-SnapKitty-Opaque-1.0
# Copyright Ahmad Ali Parr
#
# make sim   build and run the Verilator testbench (ends in PASS), writes sim/fa_dual/cycles.dat
# make figs  copy the measured cycles into docs/img and render the figures to SVG and PNG

# Optional tool env (e.g. a micromamba env with Verilator 5.054); used only if it exists.
TOOLS ?= /workspace/.tools/env/bin
ifneq ($(wildcard $(TOOLS)),)
export PATH := $(TOOLS):$(PATH)
endif

IMG := docs/img

.PHONY: all sim figs clean

all: sim figs

sim:
	sh sim/fa_dual/run.sh

sim/fa_dual/cycles.dat:
	sh sim/fa_dual/run.sh

figs: sim/fa_dual/cycles.dat
	{ echo "# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-SnapKitty-Opaque-1.0"; \
	  echo "# Copyright Ahmad Ali Parr"; cat sim/fa_dual/cycles.dat; } > $(IMG)/cycles.dat
	cd $(IMG) && gnuplot cycles.gp && gnuplot -e "fmt='png'" cycles.gp
	cd $(IMG) && dot -Tsvg dual_dataflow.dot -o dual_dataflow.svg && dot -Tpng -Gdpi=150 dual_dataflow.dot -o dual_dataflow.png

clean:
	rm -rf sim/fa_dual/work sim/fa_dual/work_iv sim/fa_dual/cycles.dat

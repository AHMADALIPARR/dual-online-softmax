#!/bin/sh
# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-SnapKitty-Opaque-1.0
# Copyright Ahmad Ali Parr
#
# Build and run tb_fa_dual with Verilator (default). Works from any directory.
# SIM=icarus only elaborates the RTL (fa_pkg, fa_exp_approx, fa_softmax_merge, fa_dual_online_softmax):
# Icarus 13 has no unpacked structs with real members, which the bound tracker in the TB uses.
set -e
cd "$(dirname "$0")"
R=../../rtl/src
RTL="$R/fa_pkg.sv $R/fa_exp_approx.sv $R/fa_online_softmax.sv $R/fa_softmax_merge.sv $R/fa_dual_online_softmax.sv"
if [ "${SIM:-verilator}" = icarus ]; then
  mkdir -p work_iv
  iverilog -g2012 -o work_iv/fa_dual.vvp -s fa_dual_online_softmax $RTL
  echo "icarus: RTL elaborated OK"
else
  rm -rf work
  verilator --binary --timing --timescale 1ns/1ps -Wall \
    -Wno-DECLFILENAME -Wno-UNUSEDSIGNAL -Wno-UNUSEDPARAM -Wno-PINCONNECTEMPTY \
    -Wno-INITIALDLY -Wno-BLKSEQ -Wno-MULTITOP -Wno-PROCASSINIT \
    --Mdir work -I$R $RTL tb_fa_dual.sv --top-module tb_fa_dual -o tb_fa_dual >/dev/null
  ./work/tb_fa_dual
fi

// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-SnapKitty-Opaque-1.0
// Copyright Ahmad Ali Parr
//
// tb_fa_dual: fa_online_softmax (original) vs fa_dual_online_softmax (ping-pong and least-loaded).
//
//   For random rows (several logit distributions, 1..BR rows, ragged last tile) every unit gets the
//   same tile stream (key tile j0 outer, row inner -- the fa_engine order) and the testbench checks:
//     1. original: final l bit-exact vs an integer model of fa_online_softmax
//     2. dual:     merged l bit-exact vs an integer model of the dual algorithm (lane split as
//                  reported by tile_lane, lazy rescale, one merge)
//     3. both:     |l - L*| <= proven bound, L* = 65536 * sum exp((x - m)/256) in real math
//     4. dual vs original: |l_dual - l_orig| reported (and bounded by the sum of both bounds)
//     5. O (one V column, from the per-lane event streams + merge betas) vs exact attention
//     6. collisions: original alpha/p writebacks in the same cycle (probe), dual snapshot vs fold
//        in the same lane in the same cycle (must never happen); p+fold combined writes counted
//   and reports cycles: original "contract" (one op at a time, as fa_engine drives it), original
//   "pipelined" (fastest collision-free issue), dual ping-pong, dual least-loaded.
//   A directed test then issues the original back-to-back (next tile_start with the last p_start)
//   to show the same-cycle l-write collision.

`timescale 1ns/1ps
/* verilator lint_off WIDTHEXPAND */
/* verilator lint_off WIDTHTRUNC */
/* verilator lint_off SYNCASYNCNET */

module tb_fa_dual;
    import fa_pkg::*;

    localparam int RW    = $clog2(BR);
    localparam int NW    = $clog2(BC + 1);
    localparam int NMAX  = 1024;
    localparam real LOG2E = 1.4426950408889634;

    logic clk;
    logic rst_n;
    initial clk = 1'b0;
    always #5 clk = ~clk;

    int errors;

    // ================================================================== original
    logic               o_init, o_tile_start, o_alpha_valid, o_p_start, o_p_valid;
    logic [RW-1:0]      o_row, o_l_rd_row;
    logic signed [31:0] o_tile_max, o_logit;
    logic [16:0]        o_alpha, o_p;
    logic signed [63:0] o_l_out;

    fa_online_softmax u_orig (
        .clk, .rst_n, .init(o_init),
        .tile_start(o_tile_start), .row(o_row), .tile_max(o_tile_max),
        .alpha_valid(o_alpha_valid), .alpha(o_alpha),
        .p_start(o_p_start), .logit(o_logit), .p_valid(o_p_valid), .p(o_p),
        .l_rd_row(o_l_rd_row), .l_out(o_l_out));

    // ================================================================== dual units (0: ping-pong, 1: least-loaded)
    logic [1:0]                 d_init, d_tile_valid, d_tile_ready, d_tile_lane, d_idle;
    logic [1:0][RW-1:0]         d_tile_row;
    logic [1:0][31:0]           d_tile_max;
    logic [1:0][BC-1:0][31:0]   d_tile_logits;
    logic [1:0][NW-1:0]         d_tile_nkeys;
    logic [1:0][1:0]            d_ev_valid, d_ev_alpha, d_st_comb, d_st_alpha, d_st_skip, d_st_stall;
    logic [1:0][1:0][RW-1:0]    d_ev_row;
    logic [1:0][1:0][16:0]      d_ev_val;
    logic [1:0]                 d_merge_start, d_merge_done, d_o_valid, d_o_out_valid;
    logic [1:0][RW-1:0]         d_merge_row, d_merge_out_row, d_l_rd_row;
    logic [1:0][31:0]           d_merge_m;
    logic [1:0][63:0]           d_merge_l, d_o0, d_o1, d_o_out;
    logic [1:0][16:0]           d_beta0, d_beta1;
    logic [1:0][1:0][63:0]      d_l_lane;
    logic [1:0][1:0][31:0]      d_m_lane;

    for (genvar u = 0; u < 2; u++) begin : g_dut
        fa_dual_online_softmax #(.PINGPONG(u == 0)) u_dual (
            .clk, .rst_n, .init(d_init[u]),
            .tile_valid(d_tile_valid[u]), .tile_ready(d_tile_ready[u]),
            .tile_row(d_tile_row[u]), .tile_max(d_tile_max[u]),
            .tile_logits(d_tile_logits[u]), .tile_nkeys(d_tile_nkeys[u]), .tile_lane(d_tile_lane[u]),
            .ev_valid(d_ev_valid[u]), .ev_alpha(d_ev_alpha[u]), .ev_row(d_ev_row[u]), .ev_val(d_ev_val[u]),
            .idle(d_idle[u]),
            .merge_start(d_merge_start[u]), .merge_row(d_merge_row[u]), .merge_done(d_merge_done[u]),
            .merge_out_row(d_merge_out_row[u]), .merge_m(d_merge_m[u]), .merge_l(d_merge_l[u]),
            .merge_beta0(d_beta0[u]), .merge_beta1(d_beta1[u]),
            .o_valid(d_o_valid[u]), .o0(d_o0[u]), .o1(d_o1[u]),
            .o_out_valid(d_o_out_valid[u]), .o_out(d_o_out[u]),
            .l_rd_row(d_l_rd_row[u]), .l_lane(d_l_lane[u]), .m_lane(d_m_lane[u]),
            .st_combined(d_st_comb[u]), .st_alpha(d_st_alpha[u]), .st_skip(d_st_skip[u]),
            .st_stall(d_st_stall[u]));
    end

    // ================================================================== test data
    int          n_rows, n_keys, n_tiles;
    int          tk = BC;                 // keys per tile actually used (1..BC)
    int          x [BR][NMAX];
    int          v [NMAX];

    // ------------------------------------------------------------------ bit-exact exp (fa_exp_approx)
    function automatic longint exp_hw(input longint delta);
        longint dcl, y8, n, f, t0, t1, base;
        dcl  = (delta > 4095) ? 4095 : delta;
        y8   = (dcl * 369) >> 8;
        n    = y8 >> 8;
        f    = y8 & 255;
        t0   = longint'(exp_lut(5'(f >> 4)));
        t1   = longint'(exp_lut(5'((f >> 4) + 1)));
        base = t0 - (((t0 - t1) * (f & 15)) >> 4);
        return (n >= 16) ? 0 : (base >> n);
    endfunction

    function automatic real exp_true(input longint delta);   // 65536 * exp(-delta/256)
        return 65536.0 * $exp(-real'(delta) / 256.0);
    endfunction

    function automatic int tile_max_of(input int r, input int t);
        int mx;
        mx = x[r][t * tk];
        for (int i = 1; i < tk; i++) if (t * tk + i < n_keys && x[r][t * tk + i] > mx) mx = x[r][t * tk + i];
        return mx;
    endfunction

    function automatic int tile_n(input int t);
        return (n_keys - t * tk < tk) ? n_keys - t * tk : tk;
    endfunction

    // ------------------------------------------------------------------ error model (from exhaustive scan)
    real rho, aabs;

    task automatic scan_exp();
        real e, err, rel;
        rho = 0.0; aabs = 0.0;
        for (int d = 0; d < 4096; d++) begin
            e   = exp_true(d);
            err = (real'(exp_hw(d)) > e) ? real'(exp_hw(d)) - e : e - real'(exp_hw(d));
            if (e >= 1024.0) begin rel = err / e; if (rel > rho) rho = rel; end
        end
        for (int d = 0; d < 4096; d++) begin
            e   = exp_true(d);
            err = (real'(exp_hw(d)) > e) ? real'(exp_hw(d)) - e : e - real'(exp_hw(d));
            if (err - rho * e > aabs) aabs = err - rho * e;
        end
        if (exp_true(4096) > aabs) aabs = exp_true(4096);    // clamp region: hw 0, true < this
        $display("exp_approx error model: |E^ - E| <= rho*E + A,  rho = %.6e  A = %.4f LSB (Q0.16)", rho, aabs);
    endtask

    // bound state of one accumulator: |l^ - L| <= R*L + B ; cnt = number of p's (L <= cnt * 2^16)
    typedef struct { real R; real B; int cnt; } bnd_t;

    function automatic bnd_t bnd_add_p(input bnd_t s);
        bnd_t o;
        o.R = (s.R > rho) ? s.R : rho; o.B = s.B + aabs; o.cnt = s.cnt + 1;
        return o;
    endfunction
    function automatic bnd_t bnd_rescale(input bnd_t s);
        bnd_t o;
        o.R = s.R + rho + s.R * rho;
        o.B = s.B * (1.0 + rho + aabs / 65536.0) + (1.0 + s.R) * real'(s.cnt) * aabs + 1.0;
        o.cnt = s.cnt;
        return o;
    endfunction

    // ------------------------------------------------------------------ integer models
    longint ref_l_orig [BR];  bnd_t ref_b_orig [BR];  int ref_resc_orig [BR];
    longint ref_l_dual [BR];  bnd_t ref_b_dual [BR];
    real    lstar [BR], ostar [BR];

    task automatic model_orig();
        for (int r = 0; r < n_rows; r++) begin
            longint m, l, mn, a;
            bnd_t b;
            m = NEG_INF; l = 0; b = '{0.0, 0.0, 0}; ref_resc_orig[r] = 0;
            for (int t = 0; t < n_tiles; t++) begin
                mn = tile_max_of(r, t); if (m > mn) mn = m;
                a  = exp_hw(mn - m);
                if (l != 0 && a != 65536) begin b = bnd_rescale(b); ref_resc_orig[r]++; end
                l  = (l * a) >> 16;
                m  = mn;
                for (int i = 0; i < tile_n(t); i++) begin
                    l += exp_hw((m >= x[r][t * tk + i]) ? m - x[r][t * tk + i] : 0);
                    b = bnd_add_p(b);
                end
            end
            ref_l_orig[r] = l; ref_b_orig[r] = b;
        end
    endtask

    // lane assignment captured from the DUT: lane_of[u][r][t]
    int lane_of [2][BR][NMAX];

    task automatic model_dual(input int u);
        for (int r = 0; r < n_rows; r++) begin
            longint m [2], l [2], mn, a, mh, ml, lh, ll;
            bnd_t b [2], bh, bl;
            for (int k = 0; k < 2; k++) begin m[k] = NEG_INF; l[k] = 0; b[k] = '{0.0, 0.0, 0}; end
            for (int t = 0; t < n_tiles; t++) begin
                int k;
                k  = lane_of[u][r][t];
                mn = tile_max_of(r, t); if (m[k] > mn) mn = m[k];
                if (m[k] != NEG_INF && mn != m[k]) begin           // lazy: only when the max moves
                    a = exp_hw(mn - m[k]);
                    l[k] = (l[k] * a) >> 16;
                    if (l[k] != 0 || b[k].cnt != 0) b[k] = bnd_rescale(b[k]);
                end
                m[k] = mn;
                for (int i = 0; i < tile_n(t); i++) begin
                    l[k] += exp_hw(m[k] - x[r][t * tk + i]);
                    b[k] = bnd_add_p(b[k]);
                end
            end
            // merge
            if (m[1] > m[0]) begin mh = m[1]; ml = m[0]; lh = l[1]; ll = l[0]; bh = b[1]; bl = b[0]; end
            else             begin mh = m[0]; ml = m[1]; lh = l[0]; ll = l[1]; bh = b[0]; bl = b[1]; end
            a = exp_hw(mh - ml);
            ref_l_dual[r] = ((ll * a) >> 16) + lh;
            if (a != 65536 && bl.cnt != 0) bl = bnd_rescale(bl);
            ref_b_dual[r].R   = (bl.R > bh.R) ? bl.R : bh.R;
            ref_b_dual[r].B   = bl.B + bh.B;
            ref_b_dual[r].cnt = bl.cnt + bh.cnt;
        end
    endtask

    task automatic model_exact();
        for (int r = 0; r < n_rows; r++) begin
            int  mx;
            real s, so, w;
            mx = x[r][0];
            for (int j = 1; j < n_keys; j++) if (x[r][j] > mx) mx = x[r][j];
            s = 0.0; so = 0.0;
            for (int j = 0; j < n_keys; j++) begin
                w = $exp(-real'(mx - x[r][j]) / 256.0);
                s += w; so += w * real'(v[j]);
            end
            lstar[r] = 65536.0 * s;
            ostar[r] = 16.0 * so / s;              // output units 1/256 (v is Q4.4)
        end
    endtask

    // ================================================================== monitors
    longint o_acc_orig [BR];
    longint o_acc_dual [2][2][BR];
    int     orig_ops_q [$];                      // encoded (kind << 24) | (row << 16) | key
    int     dual_tiles_q [2][2][$];              // per unit, per lane: tile index t*BR + r
    int     dual_key [2][2];                     // key index inside the lane's current tile
    longint orig_collisions;
    longint dual_snapfold_conflicts, dual_combined, dual_alpha_issued, dual_skips, dual_stalls;
    bit     orig_track_o;

    always @(posedge clk) if (rst_n) begin
        // original: alpha and p writebacks in the same cycle (lost l update)
        if (u_orig.ea_valid && u_orig.ep_valid) orig_collisions++;
        if (orig_track_o) begin
            if (o_alpha_valid) begin
                int op; op = orig_ops_q.pop_front();
                if ((op >> 24) != 1) begin errors++; $display("FAIL: orig op order (alpha)"); end
                o_acc_orig[(op >> 16) & 255] = (o_acc_orig[(op >> 16) & 255] * longint'(o_alpha)) >>> 16;
            end
            if (o_p_valid) begin
                int op; op = orig_ops_q.pop_front();
                if ((op >> 24) != 2) begin errors++; $display("FAIL: orig op order (p)"); end
                o_acc_orig[(op >> 16) & 255] += longint'(o_p) * longint'(v[op & 65535]);
            end
        end
        for (int u = 0; u < 2; u++) begin
            if (d_tile_valid[u] && d_tile_ready[u]) begin
                int ti; ti = int'(d_tile_row[u]) + BR * (cur_tile[u]);
                lane_of[u][d_tile_row[u]][cur_tile[u]] = int'(d_tile_lane[u]);
                dual_tiles_q[u][d_tile_lane[u]].push_back(ti);
            end
            for (int k = 0; k < 2; k++) begin
                dual_combined     += longint'(d_st_comb[u][k]);
                dual_alpha_issued += longint'(d_st_alpha[u][k]);
                dual_skips        += longint'(d_st_skip[u][k]);
                dual_stalls       += longint'(d_st_stall[u][k]);
                if (d_ev_valid[u][k]) begin
                    int ti, r, t;
                    ti = dual_tiles_q[u][k][0]; r = ti % BR; t = ti / BR;
                    if (int'(d_ev_row[u][k]) != r) begin
                        errors++; $display("FAIL: unit %0d lane %0d event row %0d, want %0d", u, k, d_ev_row[u][k], r);
                    end
                    if (d_ev_alpha[u][k]) begin
                        if (dual_key[u][k] != 0) begin errors++; $display("FAIL: alpha inside a tile"); end
                        o_acc_dual[u][k][r] = (o_acc_dual[u][k][r] * longint'(d_ev_val[u][k])) >>> 16;
                    end else begin
                        o_acc_dual[u][k][r] += longint'(d_ev_val[u][k]) * longint'(v[t * tk + dual_key[u][k]]);
                        dual_key[u][k]++;
                        if (dual_key[u][k] == tile_n(t)) begin
                            dual_key[u][k] = 0;
                            void'(dual_tiles_q[u][k].pop_front());
                        end
                    end
                end
            end
        end
    end

    // snapshot (alpha writeback) and fold in the same lane in the same cycle would put the fold in the
    // wrong frame; the design must never do it
    always @(posedge clk) if (rst_n) begin
        if (g_dut[0].u_dual.req_set[0] && g_dut[0].u_dual.g_lane[0].fold_here) dual_snapfold_conflicts++;
        if (g_dut[0].u_dual.req_set[1] && g_dut[0].u_dual.g_lane[1].fold_here) dual_snapfold_conflicts++;
        if (g_dut[1].u_dual.req_set[0] && g_dut[1].u_dual.g_lane[0].fold_here) dual_snapfold_conflicts++;
        if (g_dut[1].u_dual.req_set[1] && g_dut[1].u_dual.g_lane[1].fold_here) dual_snapfold_conflicts++;
    end

    int cur_tile [2];

    // ================================================================== drivers
    task automatic orig_reset_state();
        @(negedge clk); o_init = 1'b1;
        @(negedge clk); o_init = 1'b0;
        for (int r = 0; r < BR; r++) o_acc_orig[r] = 0;
        orig_ops_q.delete();
    endtask

    // mode 0: contract (issue, wait for valid); mode 1: pipelined (fastest collision-free);
    // mode 2: overlapped (next tile_start together with the last p_start) -- unsafe, demo only
    task automatic run_orig(input int mode, output longint cycles);
        longint c0;
        orig_reset_state();
        c0 = cyc;
        for (int t = 0; t < n_tiles; t++) begin
            for (int r = 0; r < n_rows; r++) begin
                // tile_start (in mode 2 it was already issued together with the previous last p)
                if (!(mode == 2 && (t != 0 || r != 0))) begin
                    o_tile_start = 1'b1; o_row = RW'(r); o_tile_max = tile_max_of(r, t);
                    orig_ops_q.push_back((1 << 24) | (r << 16));
                    @(negedge clk); o_tile_start = 1'b0;
                    if (mode == 0) begin while (!o_alpha_valid) @(negedge clk); end
                end
                for (int i = 0; i < tile_n(t); i++) begin
                    o_p_start = 1'b1; o_logit = x[r][t * tk + i];
                    orig_ops_q.push_back((2 << 24) | (r << 16) | (t * tk + i));
                    if (mode == 2 && i == tile_n(t) - 1 && !(t == n_tiles - 1 && r == n_rows - 1)) begin
                        int nr, nt;
                        nr = (r + 1 == n_rows) ? 0 : r + 1; nt = (r + 1 == n_rows) ? t + 1 : t;
                        o_tile_start = 1'b1; o_row = RW'(nr); o_tile_max = tile_max_of(nr, nt);
                        orig_ops_q.push_back((1 << 24) | (nr << 16));
                    end
                    @(negedge clk); o_p_start = 1'b0; o_tile_start = 1'b0;
                    if (mode == 0) begin while (!o_p_valid) @(negedge clk); end
                end
            end
        end
        // drain: the last p result lands two cycles after its issue
        @(negedge clk); @(negedge clk);
        cycles = cyc - c0;
    endtask

    longint cyc;
    always @(posedge clk) cyc <= cyc + 1;

    task automatic dual_reset_state(input int u);
        @(negedge clk); d_init[u] = 1'b1;
        @(negedge clk); d_init[u] = 1'b0;
        for (int k = 0; k < 2; k++) begin
            for (int r = 0; r < BR; r++) o_acc_dual[u][k][r] = 0;
            dual_tiles_q[u][k].delete(); dual_key[u][k] = 0;
        end
    endtask

    longint d_l_got [2][BR];
    longint d_o_got [2][BR];

    task automatic run_dual(input int u, output longint cycles_tiles, output longint cycles_total);
        longint c0;
        int got;
        dual_reset_state(u);
        c0 = cyc;
        for (int t = 0; t < n_tiles; t++) begin
            for (int r = 0; r < n_rows; r++) begin
                d_tile_valid[u] = 1'b1; d_tile_row[u] = RW'(r); d_tile_max[u] = tile_max_of(r, t);
                d_tile_nkeys[u] = NW'(tile_n(t));
                for (int i = 0; i < BC; i++) d_tile_logits[u][i] = (i < tile_n(t)) ? x[r][t * tk + i] : 0;
                cur_tile[u] = t;
                @(posedge clk);
                while (!d_tile_ready[u]) @(posedge clk);
                @(negedge clk);
            end
        end
        d_tile_valid[u] = 1'b0;
        while (!d_idle[u]) @(negedge clk);
        cycles_tiles = cyc - c0;
        // pipelined merge of all rows (one per cycle)
        got = 0;
        for (int r = 0; r < n_rows; r++) begin
            d_merge_start[u] = 1'b1; d_merge_row[u] = RW'(r);
            @(negedge clk);
            if (d_merge_done[u]) begin d_l_got[u][d_merge_out_row[u]] = d_merge_l[u]; got++; end
        end
        d_merge_start[u] = 1'b0;
        while (got < n_rows) begin
            @(negedge clk);
            if (d_merge_done[u]) begin d_l_got[u][d_merge_out_row[u]] = d_merge_l[u]; got++; end
        end
        cycles_total = cyc - c0;
        // per-row merge of the O accumulators with the held betas (functional check, not timed)
        for (int r = 0; r < n_rows; r++) begin
            d_merge_start[u] = 1'b1; d_merge_row[u] = RW'(r);
            @(negedge clk); d_merge_start[u] = 1'b0;
            while (!d_merge_done[u]) @(negedge clk);
            d_o_valid[u] = 1'b1; d_o0[u] = o_acc_dual[u][0][r]; d_o1[u] = o_acc_dual[u][1][r];
            @(negedge clk); d_o_valid[u] = 1'b0;
            d_o_got[u][r] = d_o_out[u];
        end
    endtask

    // ================================================================== checks and stats
    function automatic longint tdiv(input longint n, input longint d);
        longint q;
        if (d == 0) return 0;
        q = ((n < 0) ? -n : n) / ((d < 0) ? -d : d);
        return ((n < 0) != (d < 0)) ? -q : q;
    endfunction

    real    max_rel_orig, max_rel_dual, max_rel_dd, max_bnd_use_orig, max_bnd_use_dual;
    real    max_oerr_orig, max_oerr_dual;
    longint max_abs_dd;
    longint cyc_sum [4];
    longint n_logits_total;
    int     n_rows_checked;

    task automatic check_case(input string tag);
        longint lo, ld;
        real err, bnd, rel, oe;
        for (int r = 0; r < n_rows; r++) begin
            o_l_rd_row = RW'(r); #1;
            lo = o_l_out;
            if (lo != ref_l_orig[r]) begin
                errors++; $display("FAIL %s: orig row %0d l=%0d model=%0d", tag, r, lo, ref_l_orig[r]);
            end
            err = (real'(lo) > lstar[r]) ? real'(lo) - lstar[r] : lstar[r] - real'(lo);
            bnd = ref_b_orig[r].R * lstar[r] + ref_b_orig[r].B;
            if (err > bnd) begin errors++; $display("FAIL %s: orig row %0d err %f > bound %f", tag, r, err, bnd); end
            rel = err / lstar[r]; if (rel > max_rel_orig) max_rel_orig = rel;
            if (bnd > 0.0 && err / bnd > max_bnd_use_orig) max_bnd_use_orig = err / bnd;
            oe = real'(tdiv(o_acc_orig[r] * 16, lo)) - ostar[r]; if (oe < 0.0) oe = -oe;
            if (oe > max_oerr_orig) max_oerr_orig = oe;

            for (int u = 0; u < 2; u++) begin
                ld = d_l_got[u][r];
                if (ld != ref_l_dual[r] && u == 0) begin
                    errors++; $display("FAIL %s: dual%0d row %0d l=%0d model=%0d", tag, u, r, ld, ref_l_dual[r]);
                end
                err = (real'(ld) > lstar[r]) ? real'(ld) - lstar[r] : lstar[r] - real'(ld);
                bnd = ref_b_dual[r].R * lstar[r] + ref_b_dual[r].B;
                if (err > bnd) begin errors++; $display("FAIL %s: dual%0d row %0d err %f > bound %f", tag, u, r, err, bnd); end
                rel = err / lstar[r]; if (rel > max_rel_dual) max_rel_dual = rel;
                if (bnd > 0.0 && err / bnd > max_bnd_use_dual) max_bnd_use_dual = err / bnd;
                if (((ld > lo) ? ld - lo : lo - ld) > max_abs_dd) max_abs_dd = (ld > lo) ? ld - lo : lo - ld;
                rel = real'((ld > lo) ? ld - lo : lo - ld) / lstar[r]; if (rel > max_rel_dd) max_rel_dd = rel;
                oe = real'(tdiv(d_o_got[u][r] * 16, ld)) - ostar[r]; if (oe < 0.0) oe = -oe;
                if (oe > max_oerr_dual) max_oerr_dual = oe;
                if (oe > 8.0) begin errors++; $display("FAIL %s: dual%0d row %0d O err %f", tag, u, r, oe); end
            end
            n_rows_checked++;
        end
    endtask

    // ------------------------------------------------------------------ data generation
    int seed;
    function automatic int rnd(input int lo, input int hi);
        seed = seed ^ (seed << 13); seed = seed ^ int'(unsigned'(seed) >> 17); seed = seed ^ (seed << 5);   // xorshift32, simulator-independent
        return lo + int'(unsigned'(seed) % unsigned'(hi - lo + 1));
    endfunction

    task automatic gen(input int dsel, input int rows, input int keys);
        n_rows = rows; n_keys = keys; n_tiles = (keys + tk - 1) / tk;
        for (int j = 0; j < keys; j++) v[j] = rnd(-128, 127);
        for (int r = 0; r < rows; r++) begin
            int w;
            w = rnd(-500, 500);
            for (int j = 0; j < keys; j++) begin
                case (dsel)
                    0: x[r][j] = rnd(-1024, 1024);                          // moderate
                    1: x[r][j] = -4000 + j * 24 + rnd(-40, 40);             // rising: max moves every tile
                    2: x[r][j] = 4000 - j * 24 + rnd(-40, 40);              // falling: max never moves
                    3: x[r][j] = rnd(-20000, 20000);                        // wide: peaked, alpha -> 0
                    4: begin w += rnd(-96, 100); x[r][j] = w; end           // random walk
                    5: x[r][j] = 777;                                       // all equal (ties)
                    default: x[r][j] = rnd(-262144, 262144);                // full fa_engine logit range
                endcase
            end
        end
    endtask

    task automatic run_case(input string tag, output longint cyc_c, output longint cyc_p,
                            output longint cyc_d, output longint cyc_dl);
        longint dummy;
        model_orig(); model_exact();
        orig_track_o = 1'b1;
        run_orig(0, cyc_c);
        orig_track_o = 1'b0;
        // re-run pipelined, check l again (O tracked from the contract run)
        begin
            longint keep [BR];
            for (int r = 0; r < BR; r++) keep[r] = o_acc_orig[r];
            run_orig(1, cyc_p);
            for (int r = 0; r < BR; r++) o_acc_orig[r] = keep[r];
        end
        run_dual(0, dummy, cyc_d);
        run_dual(1, dummy, cyc_dl);
        model_dual(0);
        check_case(tag);
        // the least-loaded unit split the tiles differently: check it against its own model
        model_dual(1);
        for (int r = 0; r < n_rows; r++)
            if (d_l_got[1][r] != ref_l_dual[r]) begin
                errors++; $display("FAIL %s: dual-LL row %0d l=%0d model=%0d", tag, r, d_l_got[1][r], ref_l_dual[r]);
            end
        n_logits_total += longint'(n_rows) * longint'(n_keys);
    endtask

    // ================================================================== main
    initial begin
        int fd;
        longint cc, cp, cd, cdl;
        errors = 0; seed = 20261102; cyc = 0;
        o_init = 0; o_tile_start = 0; o_p_start = 0; o_row = '0; o_tile_max = 0; o_logit = 0; o_l_rd_row = '0;
        d_init = '0; d_tile_valid = '0; d_tile_row = '0; d_tile_max = '0; d_tile_logits = '0; d_tile_nkeys = '0;
        d_merge_start = '0; d_merge_row = '0; d_o_valid = '0; d_o0 = '0; d_o1 = '0; d_l_rd_row = '0;
        orig_collisions = 0; dual_snapfold_conflicts = 0; dual_combined = 0; dual_alpha_issued = 0;
        dual_skips = 0; dual_stalls = 0; orig_track_o = 0;
        max_rel_orig = 0; max_rel_dual = 0; max_rel_dd = 0; max_abs_dd = 0; max_bnd_use_orig = 0;
        max_bnd_use_dual = 0; max_oerr_orig = 0; max_oerr_dual = 0; n_logits_total = 0; n_rows_checked = 0;
        for (int i = 0; i < 4; i++) cyc_sum[i] = 0;
        rst_n = 0;
        repeat (3) @(posedge clk);
        rst_n = 1;
        repeat (2) @(posedge clk);

        scan_exp();

        // ---------------- random regression
        for (int tc = 0; tc < 700; tc++) begin
            int dsel;
            dsel = tc % 7;
            tk   = (tc % 3 == 2) ? rnd(1, BC) : BC;          // partial tiles exercise back-to-back alphas / stalls
            gen(dsel, rnd(1, BR), (tc % 5 == 0) ? rnd(1, 9) : rnd(1, 160));
            run_case($sformatf("case %0d dist %0d rows %0d keys %0d", tc, dsel, n_rows, n_keys), cc, cp, cd, cdl);
            cyc_sum[0] += cc; cyc_sum[1] += cp; cyc_sum[2] += cd; cyc_sum[3] += cdl;
        end
        $display("random regression: %0d cases, %0d rows, %0d logits, errors = %0d", 700, n_rows_checked, n_logits_total, errors);
        $display("  collisions: original (alpha+p writeback same cycle, safe issue) = %0d", orig_collisions);
        $display("  collisions: dual (snapshot+fold same lane same cycle)            = %0d", dual_snapfold_conflicts);
        $display("  dual: p+fold retired in one combined write = %0d, phantom alphas = %0d, lazy skips = %0d, stall cycles = %0d",
                 dual_combined, dual_alpha_issued, dual_skips, dual_stalls);
        $display("  max |l - L*|/L*: original %.3e  dual %.3e ;  max |l_dual - l_orig| = %0d LSB (%.3e rel)",
                 max_rel_orig, max_rel_dual, max_abs_dd, max_rel_dd);
        $display("  max err / proven bound: original %.3f  dual %.3f", max_bnd_use_orig, max_bnd_use_dual);
        $display("  max |O - O*| (1/256 units): original %.3f  dual %.3f", max_oerr_orig, max_oerr_dual);
        $display("  cycles (sum over cases): orig-contract %0d  orig-pipelined %0d  dual-pingpong %0d  dual-leastloaded %0d",
                 cyc_sum[0], cyc_sum[1], cyc_sum[2], cyc_sum[3]);
        if (orig_collisions != 0) begin errors++; $display("FAIL: original collided under safe issue"); end
        if (dual_snapfold_conflicts != 0) begin errors++; $display("FAIL: dual snapshot/fold conflict"); end

        tk = BC;
        // ---------------- throughput sweep (written for gnuplot)
        fd = $fopen("cycles.dat", "w");
        $fdisplay(fd, "# rows keys dist orig_contract orig_pipelined dual_pingpong dual_leastloaded");
        $display("throughput sweep (cycles incl. drain and merge):");
        $display("  rows keys dist | orig-contract orig-pipelined dual-pp dual-ll | speedup vs pipelined (pp)");
        foreach (sweep_keys[i]) begin
            for (int rr = 0; rr < 2; rr++) begin
                for (int dd = 0; dd < 3; dd++) begin
                    int rows, dsel;
                    rows = (rr == 0) ? BR : 1;
                    dsel = (dd == 0) ? 0 : (dd == 1) ? 1 : 4;
                    gen(dsel, rows, sweep_keys[i]);
                    run_case($sformatf("sweep rows %0d keys %0d dist %0d", rows, sweep_keys[i], dsel), cc, cp, cd, cdl);
                    $fdisplay(fd, "%0d %0d %0d %0d %0d %0d %0d", rows, sweep_keys[i], dsel, cc, cp, cd, cdl);
                    $display("  %4d %4d %4d | %13d %14d %7d %7d | %.2fx", rows, sweep_keys[i], dsel, cc, cp, cd, cdl,
                             real'(cp) / real'(cd));
                end
            end
        end
        $fclose(fd);

        // ---------------- directed: original with overlapped issue (next tile_start + last p_start)
        begin
            longint c, coll0;
            int bad;
            gen(1, 2, 16); model_orig();
            coll0 = orig_collisions;
            run_orig(2, c);
            bad = 0;
            for (int r = 0; r < n_rows; r++) begin
                o_l_rd_row = RW'(r); #1;
                if (o_l_out != ref_l_orig[r]) bad++;
                $display("  overlap demo row %0d: original l = %0d, correct = %0d", r, o_l_out, ref_l_orig[r]);
            end
            $display("directed overlap test (original): %0d same-cycle alpha/p writebacks, %0d of %0d rows wrong",
                     orig_collisions - coll0, bad, n_rows);
        end

        // ---------------- directed: one-key tiles, rising max -> every same-lane tile needs an alpha while one is pending
        begin
            longint s0, cc2, cp2, cd2, cdl2;
            s0 = dual_stalls; tk = 1;
            gen(1, 1, 64);
            run_case("directed stall", cc2, cp2, cd2, cdl2);
            $display("directed stall test (tk=1, rising): stall cycles = %0d, cycles orig-pipelined %0d dual-pp %0d dual-ll %0d",
                     dual_stalls - s0, cp2, cd2, cdl2);
            if (dual_stalls == s0) begin errors++; $display("FAIL: stall path not exercised"); end
            tk = BC;
        end

        $display("tb_fa_dual: errors = %0d", errors);
        if (errors == 0) $display("PASS"); else $fatal(1, "FAIL");
        $finish;
    end

    int sweep_keys [5] = '{16, 64, 256, 512, 1024};

    initial begin #2_000_000_000; $display("timeout"); $fatal(1); end
endmodule

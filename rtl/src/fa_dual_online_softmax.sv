// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-SnapKitty-Opaque-1.0
// Copyright Ahmad Ali Parr
//
// DUAL_ONLINE_SOFTMAX: two-lane (dual-processor) online softmax with phantom-logit alpha, lazy
// rescale, frame-split deferred rescale on one shared multiplier, and an exact final merge.
// See docs/dual-online-softmax.md for the algorithm, the proof and the error bound.
//
//   tile in (row, tile_max, BC logits) --dispatch--> lane 0 | lane 1      (ping-pong per row)
//
//   lane k (k = 0, 1): own m_k[row], l_k[row], one EXP_APPROX, one issue slot per cycle.
//     begin tile: m_new = max(m_k, tile_max);  m_k <= m_new
//                 if m_k != NEG_INF && m_new != m_k   (max really moved: LAZY rescale)
//                     issue the PHANTOM LOGIT m_k:  alpha = exp(-(m_new - m_k))      (same exp unit)
//                 then issue one logit per cycle:   p = exp(-(m_new - logit))
//     writeback (one exp result per lane per cycle, never two):
//       alpha : FRAME SPLIT -- l_prev = l_k[row] goes to the shared rescaler, l_k[row] <= 0
//       p     : l_k[row] <= l_k[row] + p (+ fold, if the rescaler returns for this row this cycle)
//   RESCALER (shared by both lanes, round-robin): fold = (l_prev * alpha) >> 16, added back to
//     l_k[row] DEFERRED, through the same single write (one write per register per cycle).
//     A lane holds at most one outstanding rescale; a second max change in that lane waits.
//
//   Per lane the integer result is bit-identical to fa_online_softmax fed with that lane's tiles:
//     floor(l_prev * alpha / 2^16) + sum(p) is the same integer whatever the order of the adds,
//     and a skipped rescale is a multiply by exp(0) = 65536, i.e. the identity.
//
//   merge_start(row) -> fa_softmax_merge -> merge_l = exact (m, l) of the whole row.
//   ev_* : per-lane event stream for the output accumulators (alpha before the tile's p's,
//          exactly the order fa_pv_accum expects); merge_beta0/1 then rescale o0/o1.
//
// tile_nkeys must be 1..BC. init and merge_start only while idle.

module fa_dual_online_softmax
    import fa_pkg::*;
#(
    parameter bit PINGPONG   = 1'b1,     // 1: per-row ping-pong; 0: least-loaded lane (ping-pong on ties)
    parameter int FIFO_DEPTH = 4         // tiles queued per lane (power of two)
) (
    input  logic                       clk,
    input  logic                       rst_n,

    input  logic                       init,           // m <= -inf, l <= 0 for all rows, both lanes

    input  logic                       tile_valid,
    output logic                       tile_ready,
    input  logic [$clog2(BR)-1:0]      tile_row,
    input  logic signed [31:0]         tile_max,
    input  logic [BC-1:0][31:0]        tile_logits,    // signed, logit units (1/256)
    input  logic [$clog2(BC+1)-1:0]    tile_nkeys,     // 1..BC
    output logic                       tile_lane,      // lane that takes the tile (valid with ready)

    output logic [1:0]                 ev_valid,       // per lane, one pulse per exp result
    output logic [1:0]                 ev_alpha,       // 1: rescale event, 0: p event
    output logic [1:0][$clog2(BR)-1:0] ev_row,
    output logic [1:0][16:0]           ev_val,

    output logic                       idle,

    input  logic                       merge_start,
    input  logic [$clog2(BR)-1:0]      merge_row,
    output logic                       merge_done,
    output logic [$clog2(BR)-1:0]      merge_out_row,
    output logic signed [31:0]         merge_m,
    output logic signed [63:0]         merge_l,
    output logic [16:0]                merge_beta0,
    output logic [16:0]                merge_beta1,

    input  logic                       o_valid,
    input  logic signed [63:0]         o0,
    input  logic signed [63:0]         o1,
    output logic                       o_out_valid,
    output logic signed [63:0]         o_out,

    input  logic [$clog2(BR)-1:0]      l_rd_row,
    output logic [1:0][63:0]           l_lane,         // per-lane partial l (debug / tests)
    output logic [1:0][31:0]           m_lane,

    output logic [1:0]                 st_combined,    // p and fold retired into one write
    output logic [1:0]                 st_alpha,       // phantom logit issued (max moved)
    output logic [1:0]                 st_skip,        // tile began with no rescale (lazy)
    output logic [1:0]                 st_stall        // tile held: lane's rescale still outstanding
);

    localparam int RW = $clog2(BR);
    localparam int NW = $clog2(BC + 1);
    localparam int IW = $clog2(BC);
    localparam int PW = $clog2(FIFO_DEPTH);

    typedef struct packed {
        logic [RW-1:0]        row;
        logic signed [31:0]   tmax;
        logic [BC-1:0][31:0]  logits;
        logic [NW-1:0]        nkeys;
    } tile_t;

    function automatic logic [31:0] pdelta(input logic signed [31:0] mc, input logic signed [31:0] x);
        return (mc >= x) ? 32'(mc - x) : 32'd0;
    endfunction

    // ------------------------------------------------------------------ dispatch
    logic [BR-1:0]  toggle;                      // next lane for each row
    logic [1:0]     push, full;
    logic [1:0][PW:0] load;                      // queued tiles + tile in progress
    logic           sel;

    always_comb begin
        if (PINGPONG || load[0] == load[1]) sel = toggle[tile_row];
        else                                sel = (load[1] < load[0]);
    end
    assign tile_lane  = sel;
    assign tile_ready = !full[sel] && !init;
    assign push[0]    = tile_valid && tile_ready && !sel;
    assign push[1]    = tile_valid && tile_ready &&  sel;

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            for (int r = 0; r < BR; r++) toggle[r] <= r[0];
        end else if (init) begin
            for (int r = 0; r < BR; r++) toggle[r] <= r[0];
        end else if (tile_valid && tile_ready) begin
            toggle[tile_row] <= !sel;
        end
    end

    tile_t tile_in;
    always_comb begin
        tile_in.row    = tile_row;
        tile_in.tmax   = tile_max;
        tile_in.logits = tile_logits;
        tile_in.nkeys  = tile_nkeys;
    end

    // ------------------------------------------------------------------ shared rescaler
    logic [1:0]                 req_valid;
    logic [1:0][RW-1:0]         req_row;
    logic [1:0][63:0]           req_l;
    logic [1:0][16:0]           req_a;
    logic [1:0]                 req_set;          // from lanes: snapshot this cycle
    logic [1:0][RW-1:0]         req_set_row;
    logic [1:0][63:0]           req_set_l;
    logic [1:0][16:0]           req_set_a;

    logic       rr, grant_v, grant;
    logic       fold_v, fold_lane;
    logic [RW-1:0] fold_row;
    logic [63:0]   fold_val;

    always_comb begin
        grant_v = |req_valid;
        if (req_valid[0] && req_valid[1]) grant = rr;
        else                              grant = req_valid[1];
    end

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            req_valid <= '0; req_row <= '0; req_l <= '0; req_a <= '0; rr <= 1'b0;
            fold_v <= 1'b0; fold_lane <= 1'b0; fold_row <= '0; fold_val <= '0;
        end else if (init) begin
            req_valid <= '0; fold_v <= 1'b0;
        end else begin
            fold_v <= grant_v;
            if (grant_v) begin
                fold_lane <= grant;
                fold_row  <= req_row[grant];
                fold_val  <= 64'((128'(req_l[grant]) * 128'(req_a[grant])) >> 16);
                rr        <= !grant;
            end
            for (int k = 0; k < 2; k++) begin
                if (req_set[k]) begin
                    req_valid[k] <= 1'b1;
                    req_row[k]   <= req_set_row[k];
                    req_l[k]     <= req_set_l[k];
                    req_a[k]     <= req_set_a[k];
                end else if (grant_v && grant == 1'(k)) begin
                    req_valid[k] <= 1'b0;
                end
            end
        end
    end

    // ------------------------------------------------------------------ lanes
    logic [1:0][63:0] mrg_l;
    logic [1:0][31:0] mrg_m;
    logic [1:0]       lane_idle;

    for (genvar k = 0; k < 2; k++) begin : g_lane
        logic signed [31:0] m [BR];
        logic signed [63:0] l [BR];

        // tile FIFO
        tile_t           fifo [FIFO_DEPTH];
        logic [PW:0]     wr_ptr, rd_ptr;
        logic            empty, pop;
        tile_t           head;
        assign empty   = (wr_ptr == rd_ptr);
        assign full[k] = (wr_ptr[PW] != rd_ptr[PW]) && (wr_ptr[PW-1:0] == rd_ptr[PW-1:0]);
        assign head    = fifo[rd_ptr[PW-1:0]];

        // sequencer
        logic               busy, pend;
        logic [RW-1:0]      cur_row;
        logic signed [31:0] cur_m;
        logic [BC-1:0][31:0] cur_logits;
        logic [NW-1:0]      cur_n;
        logic [IW:0]        idx;

        logic signed [31:0] h_m_old, h_m_new;
        logic               need_alpha, begin_ok;
        logic               iss_valid, iss_alpha;
        logic [RW-1:0]      iss_row;
        logic [31:0]        iss_delta;

        always_comb begin
            h_m_old    = m[head.row];
            h_m_new    = (head.tmax > h_m_old) ? head.tmax : h_m_old;
            need_alpha = (h_m_old != NEG_INF) && (h_m_new != h_m_old);
            begin_ok   = !busy && !empty && !(need_alpha && pend) && !init;
            pop        = begin_ok;

            iss_valid = 1'b0; iss_alpha = 1'b0; iss_row = cur_row; iss_delta = '0;
            if (busy) begin
                iss_valid = 1'b1;
                iss_delta = pdelta(cur_m, $signed(cur_logits[idx[IW-1:0]]));
            end else if (begin_ok) begin
                iss_valid = 1'b1;
                iss_row   = head.row;
                iss_alpha = need_alpha;
                iss_delta = need_alpha ? 32'(h_m_new - h_m_old) : pdelta(h_m_new, $signed(head.logits[0]));
            end
        end

        assign st_alpha[k] = begin_ok && need_alpha;
        assign st_skip[k]  = begin_ok && !need_alpha;
        assign st_stall[k] = !busy && !empty && need_alpha && pend && !init;
        assign load[k]     = (wr_ptr - rd_ptr) + (PW+1)'(busy);

        always_ff @(posedge clk or negedge rst_n) begin
            if (!rst_n) begin
                wr_ptr <= '0; rd_ptr <= '0;
                busy <= 1'b0; cur_row <= '0; cur_m <= NEG_INF; cur_logits <= '0; cur_n <= '0; idx <= '0;
            end else if (init) begin
                wr_ptr <= '0; rd_ptr <= '0; busy <= 1'b0;
            end else begin
                if (push[k]) begin
                    fifo[wr_ptr[PW-1:0]] <= tile_in;
                    wr_ptr <= wr_ptr + 1'b1;
                end
                if (busy) begin
                    idx  <= idx + 1'b1;
                    busy <= (32'(idx) + 32'd1 < 32'(cur_n));
                end else if (begin_ok) begin
                    rd_ptr     <= rd_ptr + 1'b1;
                    cur_row    <= head.row;
                    cur_m      <= h_m_new;
                    cur_logits <= head.logits;
                    cur_n      <= head.nkeys;
                    if (need_alpha) begin idx <= '0; busy <= 1'b1; end
                    else            begin idx <= 1;  busy <= (head.nkeys > 1); end
                end
            end
        end

        // exp unit (phantom logits and real logits share it)
        logic          wb_valid, wb_alpha;
        logic [RW-1:0] wb_row;
        logic [16:0]   wb_res;
        fa_exp_approx u_exp (
            .clk, .rst_n,
            .in_valid(iss_valid),
            .delta(iss_delta),
            .out_valid(wb_valid),
            .result(wb_res)
        );
        always_ff @(posedge clk or negedge rst_n) begin
            if (!rst_n) begin wb_alpha <= 1'b0; wb_row <= '0; end
            else if (iss_valid) begin wb_alpha <= iss_alpha; wb_row <= iss_row; end
        end

        assign ev_valid[k] = wb_valid;
        assign ev_alpha[k] = wb_alpha;
        assign ev_row[k]   = wb_row;
        assign ev_val[k]   = wb_res;

        // snapshot to the rescaler (frame split)
        assign req_set[k]     = wb_valid && wb_alpha;
        assign req_set_row[k] = wb_row;
        assign req_set_l[k]   = l[wb_row];
        assign req_set_a[k]   = wb_res;

        logic fold_here;
        assign fold_here      = fold_v && (fold_lane == 1'(k));
        assign st_combined[k] = fold_here && wb_valid && !wb_alpha && (wb_row == fold_row);

        // the only writer of m/l: one write per register per cycle, all same-cycle terms summed
        always_ff @(posedge clk or negedge rst_n) begin
            if (!rst_n) begin
                for (int r = 0; r < BR; r++) begin m[r] <= NEG_INF; l[r] <= '0; end
                pend <= 1'b0;
            end else if (init) begin
                for (int r = 0; r < BR; r++) begin m[r] <= NEG_INF; l[r] <= '0; end
                pend <= 1'b0;
            end else begin
                if (begin_ok) m[head.row] <= h_m_new;
                if (begin_ok && need_alpha) pend <= 1'b1;
                else if (fold_here)         pend <= 1'b0;
                for (int r = 0; r < BR; r++) begin
                    logic hit_p, hit_f;
                    hit_p = wb_valid && !wb_alpha && (wb_row == RW'(r));
                    hit_f = fold_here && (fold_row == RW'(r));
                    if (wb_valid && wb_alpha && wb_row == RW'(r))
                        l[r] <= '0;                                      // new frame
                    else if (hit_p || hit_f)
                        l[r] <= l[r] + (hit_p ? 64'(wb_res) : 64'd0) + (hit_f ? fold_val : 64'd0);
                end
            end
        end

        assign lane_idle[k] = empty && !busy && !wb_valid && !pend && !req_valid[k];
        assign mrg_l[k]     = l[merge_row];
        assign mrg_m[k]     = m[merge_row];
        assign l_lane[k]    = l[l_rd_row];
        assign m_lane[k]    = m[l_rd_row];
    end

    assign idle = &lane_idle && !fold_v;

    // ------------------------------------------------------------------ exact merge
    fa_softmax_merge u_merge (
        .clk, .rst_n,
        .start(merge_start), .row(merge_row),
        .m0($signed(mrg_m[0])), .l0($signed(mrg_l[0])),
        .m1($signed(mrg_m[1])), .l1($signed(mrg_l[1])),
        .done(merge_done), .out_row(merge_out_row), .m_out(merge_m), .l_out(merge_l),
        .beta0(merge_beta0), .beta1(merge_beta1),
        .o_valid, .o0, .o1, .o_out_valid, .o_out
    );

endmodule : fa_dual_online_softmax

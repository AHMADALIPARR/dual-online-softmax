// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-SnapKitty-Opaque-1.0
// Copyright Ahmad Ali Parr
//
// SOFTMAX_MERGE: exact merge of two partial online-softmax states (m0, l0), (m1, l1) of one row.
//
//   m    = max(m0, m1)
//   beta = exp(-(m - min(m0, m1)))          (one EXP_APPROX evaluation: the "phantom logit" min(m0, m1))
//   l    = ((l_lo * beta) >> 16) + l_hi      (the lane holding the larger max is never rescaled)
//   beta0 / beta1 = per-lane output-accumulator factors (Q0.16); the larger-max lane gets 65536 exactly.
//
// A lane that saw no tile has m == NEG_INF and l == 0; it merges to beta = 0 and contributes nothing.
// Fully pipelined: start may pulse every cycle; done/l/m/beta follow two cycles later with out_row.
//
// After done, the held beta0/beta1 rescale the output accumulators column by column (one exp per row,
// d multiplies -- the batched part of the rescale):
//   o_valid, o0, o1 -> o_out_valid, o_out = ((o0 * beta0) >>> 16) + ((o1 * beta1) >>> 16)   (one cycle)

module fa_softmax_merge
    import fa_pkg::*;
(
    input  logic                       clk,
    input  logic                       rst_n,

    input  logic                       start,
    input  logic [$clog2(BR)-1:0]      row,
    input  logic signed [31:0]         m0,
    input  logic signed [63:0]         l0,
    input  logic signed [31:0]         m1,
    input  logic signed [63:0]         l1,

    output logic                       done,           // pulse
    output logic [$clog2(BR)-1:0]      out_row,
    output logic signed [31:0]         m_out,
    output logic signed [63:0]         l_out,
    output logic [16:0]                beta0,
    output logic [16:0]                beta1,

    input  logic                       o_valid,
    input  logic signed [63:0]         o0,
    input  logic signed [63:0]         o1,
    output logic                       o_out_valid,    // pulse
    output logic signed [63:0]         o_out
);

    // stage 0: pick the larger max, launch the single exp
    logic               hi_is_1;
    logic signed [31:0] m_hi, m_lo;
    always_comb begin
        hi_is_1 = (m1 > m0);
        m_hi    = hi_is_1 ? m1 : m0;
        m_lo    = hi_is_1 ? m0 : m1;
    end

    logic        eb_valid;
    logic [16:0] eb_res;
    fa_exp_approx u_exp_beta (
        .clk, .rst_n,
        .in_valid(start),
        .delta(32'(m_hi - m_lo)),
        .out_valid(eb_valid),
        .result(eb_res)
    );

    logic                  s1_hi_is_1;
    logic [$clog2(BR)-1:0] s1_row;
    logic signed [31:0]    s1_m;
    logic signed [63:0]    s1_l_hi, s1_l_lo;

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            s1_hi_is_1 <= 1'b0; s1_row <= '0; s1_m <= NEG_INF; s1_l_hi <= '0; s1_l_lo <= '0;
        end else if (start) begin
            s1_hi_is_1 <= hi_is_1;
            s1_row     <= row;
            s1_m       <= m_hi;
            s1_l_hi    <= hi_is_1 ? l1 : l0;
            s1_l_lo    <= hi_is_1 ? l0 : l1;
        end
    end

    // stage 1: one rescale (the lower-max lane only) + add
    logic signed [63:0] lo_scaled;
    always_comb lo_scaled = 64'((128'(s1_l_lo) * 128'(eb_res)) >> 16);

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            done <= 1'b0; out_row <= '0; m_out <= NEG_INF; l_out <= '0;
            beta0 <= 17'd65536; beta1 <= 17'd65536;
        end else begin
            done <= eb_valid;
            if (eb_valid) begin
                out_row <= s1_row;
                m_out   <= s1_m;
                l_out   <= lo_scaled + s1_l_hi;
                beta0   <= s1_hi_is_1 ? eb_res : 17'd65536;
                beta1   <= s1_hi_is_1 ? 17'd65536 : eb_res;
            end
        end
    end

    // output-accumulator merge with the held betas
    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            o_out_valid <= 1'b0; o_out <= '0;
        end else begin
            o_out_valid <= o_valid;
            if (o_valid)
                o_out <= 64'((128'(o0) * $signed({111'b0, beta0})) >>> 16)
                       + 64'((128'(o1) * $signed({111'b0, beta1})) >>> 16);
        end
    end

endmodule : fa_softmax_merge

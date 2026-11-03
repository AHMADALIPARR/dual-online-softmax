// ONLINE_SOFTMAX: running max m[row] and running denominator l[row] for BR rows.
//
//   tile_start(row, tile_max): m_new = max(m, tile_max); alpha = exp(-(m_new - m));
//                              l <= l * alpha            -> alpha_valid, alpha
//   p_start(logit)           : p = exp(-(m_cur - logit)); l <= l + p   -> p_valid, p
//
// One operation at a time. alpha/p are held stable until the next operation.

module fa_online_softmax
    import fa_pkg::*;
(
    input  logic                       clk,
    input  logic                       rst_n,

    input  logic                       init,           // m <= -inf, l <= 0 for all rows

    input  logic                       tile_start,
    input  logic [$clog2(BR)-1:0]      row,
    input  logic signed [31:0]         tile_max,
    output logic                       alpha_valid,    // pulse
    output logic [16:0]                alpha,

    input  logic                       p_start,
    input  logic signed [31:0]         logit,
    output logic                       p_valid,        // pulse
    output logic [16:0]                p,

    input  logic [$clog2(BR)-1:0]      l_rd_row,
    output logic signed [63:0]         l_out
);

    logic signed [31:0] m [BR];
    logic signed [63:0] l [BR];

    logic [$clog2(BR)-1:0] row_q;
    logic signed [31:0]    cur_m;

    // alpha exponent unit
    logic signed [31:0] m_old, m_new;
    logic               ea_valid;
    logic [16:0]        ea_res;
    always_comb begin
        m_old = m[row];
        m_new = (tile_max > m_old) ? tile_max : m_old;
    end
    fa_exp_approx u_exp_alpha (
        .clk, .rst_n,
        .in_valid(tile_start),
        .delta(32'(m_new - m_old)),
        .out_valid(ea_valid),
        .result(ea_res)
    );

    // probability exponent unit
    logic signed [31:0] pdelta;
    logic               ep_valid;
    logic [16:0]        ep_res;
    always_comb pdelta = (cur_m >= logit) ? (cur_m - logit) : 32'sd0;
    fa_exp_approx u_exp_p (
        .clk, .rst_n,
        .in_valid(p_start),
        .delta(32'(pdelta)),
        .out_valid(ep_valid),
        .result(ep_res)
    );

    logic alpha_valid_q, p_valid_q;
    logic [16:0] alpha_q, p_q;

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            for (int i = 0; i < BR; i++) begin m[i] <= NEG_INF; l[i] <= '0; end
            row_q <= '0;
            cur_m <= NEG_INF;
            alpha_valid_q <= 1'b0; p_valid_q <= 1'b0;
            alpha_q <= '0; p_q <= '0;
        end else begin
            alpha_valid_q <= 1'b0;
            p_valid_q     <= 1'b0;

            if (init) begin
                for (int i = 0; i < BR; i++) begin m[i] <= NEG_INF; l[i] <= '0; end
            end

            if (tile_start) begin
                row_q   <= row;
                cur_m   <= m_new;
                m[row]  <= m_new;
            end

            if (ea_valid) begin
                l[row_q]      <= 64'((128'(l[row_q]) * 128'(ea_res)) >> 16);
                alpha_q       <= ea_res;
                alpha_valid_q <= 1'b1;
            end

            if (ep_valid) begin
                l[row_q]    <= l[row_q] + 64'(ep_res);
                p_q         <= ep_res;
                p_valid_q   <= 1'b1;
            end
        end
    end

    assign alpha_valid = alpha_valid_q;
    assign alpha       = alpha_q;
    assign p_valid     = p_valid_q;
    assign p           = p_q;
    assign l_out       = l[l_rd_row];

endmodule : fa_online_softmax

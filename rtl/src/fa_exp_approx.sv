// EXP_APPROX: exp(-delta/256) in Q0.16 (65536 == 1.0) for delta >= 0.
//   y = delta * log2(e);  n = integer part, f = fraction
//   2^-f  ~ 16-segment LUT with linear interpolation;  result = 2^-f >> n
// One cycle of latency.

module fa_exp_approx
    import fa_pkg::*;
(
    input  logic               clk,
    input  logic               rst_n,

    input  logic               in_valid,
    input  logic [31:0]        delta,       // non-negative, logit units (1/256)

    output logic               out_valid,
    output logic [16:0]        result
);

    logic [11:0] dcl;
    logic [23:0] y8;
    logic [4:0]  n;
    logic [7:0]  f;
    logic [16:0] t0, t1, diff, base;
    logic [20:0] interp;
    logic [16:0] res_c;

    always_comb begin
        dcl  = (delta > 32'd4095) ? 12'd4095 : delta[11:0];
        y8   = 24'(dcl) * 24'd369 >> 8;
        n    = y8[12:8];
        f    = y8[7:0];
        t0   = exp_lut({1'b0, f[7:4]});
        t1   = exp_lut({1'b0, f[7:4]} + 5'd1);
        diff   = t0 - t1;
        interp = diff * f[3:0];
        base   = t0 - 17'(interp >> 4);
        res_c = (y8[23:8] >= 16'd16) ? 17'd0 : (base >> n);
    end

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            out_valid <= 1'b0;
            result    <= '0;
        end else begin
            out_valid <= in_valid;
            result    <= res_c;
        end
    end

endmodule : fa_exp_approx

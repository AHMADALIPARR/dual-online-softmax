// Shared constants and lookup tables for the FlashAttention engine.
// Numeric formats are defined bit-exactly by software/fa_int_model.py.

package fa_pkg;

    localparam int D_MAX   = 16;            // max head dimension
    localparam int BR      = 4;             // query rows per tile
    localparam int BC      = 4;             // keys per tile (part of the numeric definition)
    localparam int NEG_INF = -(1 << 30);    // initial running max

    localparam int FA_TRAP = -2;            // SUBLEQ 'a' operand that selects the FA opcode

    // exp(-i/16 in log2 units) table, Q0.16: round(65536 * 2^(-i/16)), i = 0..16
    function automatic logic [16:0] exp_lut(input logic [4:0] i);
        case (i)
            5'd0:  exp_lut = 17'd65536;
            5'd1:  exp_lut = 17'd62757;
            5'd2:  exp_lut = 17'd60097;
            5'd3:  exp_lut = 17'd57549;
            5'd4:  exp_lut = 17'd55109;
            5'd5:  exp_lut = 17'd52773;
            5'd6:  exp_lut = 17'd50535;
            5'd7:  exp_lut = 17'd48393;
            5'd8:  exp_lut = 17'd46341;
            5'd9:  exp_lut = 17'd44376;
            5'd10: exp_lut = 17'd42495;
            5'd11: exp_lut = 17'd40693;
            5'd12: exp_lut = 17'd38968;
            5'd13: exp_lut = 17'd37316;
            5'd14: exp_lut = 17'd35734;
            5'd15: exp_lut = 17'd34219;
            default: exp_lut = 17'd32768;
        endcase
    endfunction

    // round(256 / sqrt(d)), d = 1..16
    function automatic logic [8:0] scale_q(input logic [4:0] d);
        case (d)
            5'd1:  scale_q = 9'd256;
            5'd2:  scale_q = 9'd181;
            5'd3:  scale_q = 9'd148;
            5'd4:  scale_q = 9'd128;
            5'd5:  scale_q = 9'd114;
            5'd6:  scale_q = 9'd105;
            5'd7:  scale_q = 9'd97;
            5'd8:  scale_q = 9'd91;
            5'd9:  scale_q = 9'd85;
            5'd10: scale_q = 9'd81;
            5'd11: scale_q = 9'd77;
            5'd12: scale_q = 9'd74;
            5'd13: scale_q = 9'd71;
            5'd14: scale_q = 9'd68;
            5'd15: scale_q = 9'd66;
            default: scale_q = 9'd64;
        endcase
    endfunction

endpackage : fa_pkg

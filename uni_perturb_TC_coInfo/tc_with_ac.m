function output = tc_with_ac(l)

tc_bit_num = 2 * l * l;
full_bit_num = tc_bit_num + 1;
% (tc_bit -- ref. bit)

output = zeros(2, full_bit_num, full_bit_num, 'logical');

for i = 1: l
    for j = 1: l
        plq_idx = i + (j - 1) * l;
        p1 = i + (j - 1) * 2 * l;
        p2 = p1 + l;

        p3 = p1 + 2 * l;
        if p3 > tc_bit_num
            p3 = p3 - tc_bit_num;
        end
          
        p4 = p2 + 1;
        if mod(p2, l) == 0
            p4 = p2 + 1 - l;
        end
        p_list = [p1, p2, p3, p4];


        if i == l && j == l
            zz_loop = [(1 + l): (2 * l): (2 * l^2), full_bit_num];
            output(2, plq_idx, zz_loop) = 1;
        else
            output(2, plq_idx, p_list) = 1;
        end
    end
end


for i = 1: l
    for j = 1: l
        vtx_idx = l^2 + i + (j - 1) * l;

        v1 = i + (j - 1) * 2 * l;
        v2 = v1 + l;
        v3 = v1 - 1;
        if v3 < (j - 1) * 2 * l + 1
            v3 = v1 - 1 + l;
        end
        v4 = v1 - l;
        if v4 <= 0
            v4 = v4 + tc_bit_num;
        end

        v_list = [v1, v2, v3, v4];

        if i == l && j == l
            xx_loop = [l + (1: l), full_bit_num];
            output(1, vtx_idx, xx_loop) = 1;
        else
            output(1, vtx_idx, v_list) = 1;
        end
    end
end

z_loop = 1: l;
output(2, full_bit_num, z_loop) = 1;


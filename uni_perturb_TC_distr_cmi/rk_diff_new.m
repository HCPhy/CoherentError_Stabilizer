function output = rk_diff_new(l, pn)
    % Precompute the toric code state
    [state, lx_idx, lz_idx] = tc_state(l);

    % Apply ECC procedure: add noise and evaluate
    state = uni_noise(l, state, pn);

    % Precompute plaquette and vertex index lists once for commutation evaluation
    [plaq_lists, vtx_lists] = precompute_indices(l);

    com_mat = commute_eval(l, state, plaq_lists, vtx_lists);

    % Find null space over GF(2)
    sol = null_gf2(com_mat');

    % Compute the operational conditional mutual information
    output = op_cond_mi(l, sol);
end

%% Calculating the opCMI
function output = op_cond_mi(l, sol)
    % Precompute domain sets once
    domain_shifts = (0:(l^2-1)) + 1;
    half_size = l^2 / 2;

    domainA_sets = cell(l,1);
    domainB_sets = cell(l,1);

    for d = 1:l
        disp_idx = d * l;
        domainAp = mod(domain_shifts + disp_idx - 1, l^2) + 1;
        domainBp = mod(domainAp + half_size - 1, l^2) + 1;

        domainA_sets{d} = [domainAp, domainAp + l^2];
        domainB_sets{d} = [domainBp, domainBp + l^2];
    end

    output_val = 0;
    for d = 1:l
        domainA = domainA_sets{d};
        domainB = domainB_sets{d};

        matA = sol(domainA, :);
        matB = sol(domainB, :);
        matAB = sol([domainA, domainB], :);

        output_val = output_val + bitrank_mex(matA) + bitrank_mex(matB) - bitrank_mex(matAB);
    end

    output = output_val / l;
end

%% Generate stabilizer tableau of the toric code state
function [output, lx_idx, lz_idx] = tc_state(l)
    tc_bit_num = 2 * l * l;
    full_bit_num = tc_bit_num;

    output = zeros(2, full_bit_num, full_bit_num, 'logical');

    % Z-plaquette stabilizers
    for i = 1:l
        for j = 1:l
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
                loop = (1 + l): (2 * l): (2 * l^2);
                output(2, plq_idx, loop) = 1;
                lz_idx = plq_idx;
            else
                output(2, plq_idx, p_list) = 1;
            end
        end
    end

    % X-vertex stabilizers
    for i = 1:l
        for j = 1:l
            vtx_idx = l^2 + i + (j - 1) * l;

            v1 = i + (j - 1) * 2 * l;
            v2 = v1 + l;
            v3 = v1 - 1;
            if v3 < (j - 1) * 2 * l + 1
                v3 = v3 + l;
            end
            v4 = v1 - l;
            if v4 <= 0
                v4 = v4 + tc_bit_num;
            end

            v_list = [v1, v2, v3, v4];

            if i == l && j == l
                dual_loop = 1: (2 * l) : (2 * l^2);
                output(1, vtx_idx, dual_loop) = 1;
                lx_idx = vtx_idx;
            else
                output(1, vtx_idx, v_list) = 1;
            end
        end
    end
end

%% Apply 4-body random Clifford noise as "uni_noise"
function output = uni_noise(l, state, pn)
    tc_bit_num = 2 * l * l;

    xs = (1:4)*2;
    zs = xs - 1;

    [I, J] = ndgrid(1:l, 1:l);
    idx_mask = ~(I == l & J == l) & (rand(l) <= pn);
    idx_list = find(idx_mask);

    % Get the number of symplectic transformations without storing them
    nsym4 = numberofsymplectic(4);

    function apply_random_clifford(p_list)
        x_stab = squeeze(state(1, :, p_list));
        z_stab = squeeze(state(2, :, p_list));

        temp = zeros(8, tc_bit_num);
        temp(xs, :) = x_stab.';
        temp(zs, :) = z_stab.';

        % Instead of precomputing all, just compute one random G on-demand
        cliff_idx = randi(nsym4);
        g = double(symplectic(cliff_idx,4)); 

        temp = mod(g * temp, 2);

        state(1, :, p_list) = temp(xs, :)';
        state(2, :, p_list) = temp(zs, :)';
    end

    for idx = idx_list'
        i = I(idx);
        j = J(idx);

        p1 = i + (j - 1)*2*l;
        p2 = p1 + l;
        p3 = mod(p1 + 2*l - 1, tc_bit_num) + 1;
        p4 = p2 + 1;
        if mod(p2, l) == 0
            p4 = p2 + 1 - l;
        end
        p_list = [p1, p2, p3, p4];

        apply_random_clifford(p_list);
        apply_random_clifford(p_list);  % As per original code
    end

    output = state;
end
%% Precompute Plaquette and Vertex Indices
function [plaq_lists, vtx_lists] = precompute_indices(l)
    tc_bit_num = 2*l*l;
    plaq_lists = zeros(l*l,4);
    vtx_lists = zeros(l*l,4);

    for i = 1:l
        for j = 1:l
            idx = i + (j-1)*l;
            % Plaquette
            p1 = i + (j - 1)*2*l;
            p2 = p1 + l;
            p3 = p1 + 2*l;
            if p3 > tc_bit_num, p3 = p3 - tc_bit_num; end
            p4 = p2 + 1;
            if mod(p2, l) == 0
                p4 = p2 + 1 - l;
            end
            plaq_lists(idx,:) = [p1, p2, p3, p4];

            % Vertex
            v1 = i + (j - 1)*2*l;
            v2 = v1 + l;
            v3 = v1 - 1;
            if v3 < (j - 1)*2*l + 1, v3 = v3 + l; end
            v4 = v1 - l;
            if v4 <= 0, v4 = v4 + tc_bit_num; end
            vtx_lists(idx,:) = [v1, v2, v3, v4];
        end
    end
end

%% Evaluation of the commutator matrix
function output = commute_eval(l, state, plaq_lists, vtx_lists)
    stab_num = size(state, 2);

    % Extract X and Z planes
    xstates = squeeze(state(1,:,:)); % stab_num x (2*l*l)
    zstates = squeeze(state(2,:,:));

    check_num = 2*l*l;
    output = zeros(check_num, stab_num);

    % Plaquette checks (Z-type)
    % Skip the last (l,l) plaquette condition as original code does
    for idx = 1:(l*l)
        if idx == l*l && (true) % original code: if i==l && j==l continue
            continue;
        end
        ps = plaq_lists(idx,:);
        output(idx,:) = mod(sum(xstates(:,ps),2)', 2);
    end

    % Vertex checks (X-type)
    for idx = 1:l*l
        vs = vtx_lists(idx,:);
        output(l*l+idx,:) = mod(sum(zstates(:,vs),2)', 2);
    end
end

%% Check measurement (not extensively optimized, but can be improved if needed)
function output = check_measure(l, state)
    tc_bit_num = 2*l*l;

    % Convert to uint8 for faster XOR
    output = uint8(state);

    for i = 1:l
        for j = 1:l
            if i == l && j == l
                continue
            end

            % Z-plaquette
            p1 = i + (j - 1)*2*l;
            p2 = p1 + l;
            p3 = mod(p1 + 2*l - 1, tc_bit_num) + 1;
            p4 = mod(p2 - 1 + 1, l) + floor((p2 - 1)/l)*l + 1;
            if mod(p2, l) == 0
                p4 = p2 + 1 - l;
            end
            p_list = [p1, p2, p3, p4];

            sigs = mod(sum(output(1, :, p_list), 3), 2);
            sig_ids = find(sigs);
            if ~isempty(sig_ids)
                piv = sig_ids(1);
                piv_slice = output(:, piv, :);
                if numel(sig_ids) > 1
                    all_ids = sig_ids(2:end);
                    rep_piv = repmat(piv_slice, [1, numel(all_ids), 1]);
                    output(:, all_ids, :) = bitxor(output(:, all_ids, :), rep_piv);
                end
                output(:, piv, :) = 0;
                output(2, piv, p_list) = 1;
            end

            % X-vertex
            v1 = i + (j - 1)*2*l;
            v2 = v1 + l;
            v3 = v1 - 1 + l*(v1 - 1 < (j - 1)*2*l + 1);
            v4 = v1 - l + tc_bit_num*(v1 - l <= 0);
            v_list = [v1, v2, v3, v4];

            sigs = mod(sum(output(2, :, v_list), 3), 2);
            sig_ids = find(sigs);

            if ~isempty(sig_ids)
                piv = sig_ids(1);
                piv_slice = output(:, piv, :);
                if numel(sig_ids) > 1
                    all_ids = sig_ids(2:end);
                    rep_piv = repmat(piv_slice, [1, numel(all_ids), 1]);
                    output(:, all_ids, :) = bitxor(output(:, all_ids, :), rep_piv);
                end
                output(:, piv, :) = 0;
                output(1, piv, v_list) = 1;
            end
        end
    end

    % Convert back to logical
    output = logical(output);
end

%% Direct sum of two square matrices
function output = directsum(m1, m2)
[n1, ~] = size(m1);
[n2, ~] = size(m2);
output = zeros(n1 + n2, n1 + n2, 'int8');
output(1:n1, 1:n1) = m1;
output(n1+1:n1+n2, n1+1:n1+n2) = m2;
end

%% Other necessary functions (symplectic, numberofsymplectic, etc.)
% These are unchanged for optimization since they may be external or complex.

function x = numberofcosets(n)
x = 2^(2*n - 1) * (2^(2*n) - 1);
end

function x = numberofsymplectic(n)
x = 1;
for j = 1:n
    x = x * numberofcosets(j);
end
end

function output = symplectic(i, n)
% Same as original code
nn = 2 * n;
s = 2^nn - 1;
k = mod(i, s) + 1;
i = floor(i / s);
f1 = int2bits(k, nn);
e1 = zeros(1, nn, 'int8');
e1(1) = 1;
T = findtransvection(e1, f1);
bits = int2bits(mod(i, 2^(nn-1)), nn - 1);
eprime = e1; eprime(2:nn) = bits;
h0 = transvection(T(1, :), eprime);
h0 = transvection(T(2, :), h0);
if bits(1) == 1
    f1 = zeros(1, nn, 'int8');
end
if n ~= 1
    g_rest = symplectic(floor(i / 2^(nn - 1)), n - 1);
    g = directsum(int8([1 0;0 1]), g_rest);
else
    g = int8([1 0;0 1]);
end
for j = 1:nn
    g(j, :) = transvection(T(1, :), g(j, :));
    g(j, :) = transvection(T(2, :), g(j, :));
    g(j, :) = transvection(h0, g(j, :));
    g(j, :) = transvection(f1, g(j, :));
end
output = g;
end

function output = int2bits(i, n)
output = zeros(1, n, 'int8');
for j = 1:n
    output(j) = bitget(i, j);
end
end

function output = findtransvection(x, y)
n = length(x);
output = zeros(2, n, 'int8');
if isequal(x, y)
    return;
end
if inner(x, y) == 1
    output(1, :) = mod(x + y, 2);
    return;
end
z = zeros(1, n, 'int8');
for i = 1:(n/2)
    ii = 2*(i-1)+1;
    if ((x(ii)+x(ii+1))~=0)&&((y(ii)+y(ii+1))~=0)
        z(ii) = mod(x(ii)+y(ii),2);
        z(ii+1)=mod(x(ii+1)+y(ii+1),2);
        if (z(ii)+z(ii+1))==0
            z(ii+1)=1;
            if x(ii)~=x(ii+1)
                z(ii)=1;
            end
        end
        output(1,:)=mod(x+z,2);
        output(2,:)=mod(y+z,2);
        return;
    end
end
for i=1:(n/2)
    ii=2*(i-1)+1;
    if ((x(ii)+x(ii+1))~=0)&&((y(ii)+y(ii+1))==0)
        if x(ii)==x(ii+1)
            z(ii+1)=1;
        else
            z(ii+1)=x(ii);
            z(ii)=x(ii+1);
        end
        break;
    end
end
for i=1:(n/2)
    ii=2*(i-1)+1;
    if ((x(ii)+x(ii+1))==0)&&((y(ii)+y(ii+1))~=0)
        if y(ii)==y(ii+1)
            z(ii+1)=1;
        else
            z(ii+1)=y(ii);
            z(ii)=y(ii+1);
        end
        break;
    end
end
output(1,:)=mod(x+z,2);
output(2,:)=mod(y+z,2);
end

function vout = transvection(k, v)
vout = mod(v + inner(k, v)*k, 2);
end

function t = inner(v, w)
n = length(v);
if mod(n, 2) ~= 0
    error('Vectors must have even length');
end
indices = 1:2:n-1;
t = sum(v(indices).*w(indices+1) + w(indices).*v(indices+1));
t = mod(t, 2);
end

function output = bits2int(b, nn)
output = 0;
for j=1:nn
    output = output + b(j)*2^(j-1);
end
end

% The null_gf2, bitrank_mex functions are assumed to be provided externally.
% End of optimized code.
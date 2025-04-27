function output = rk_diff(l, pn)


state = tc_state_with_ac(l);
%Commutation Check
% 
% x = squeeze(state(1, :, :));
% z = squeeze(state(2, :, :));
% 
% find(mod(x * z' + z * x',2))





% ECC procedure
state = uni_noise(l, state, pn);
state = check_measure(l, state);

output = co_info(l, state);

% %post ECC check
% xc = squeeze(state(1, :, :));
% zc = squeeze(state(2, :, :));
% xz_matc = logical([xc,zc]);
% 
% output = bitrank_mex(xz_matc);





end





function output = co_info(l, state)

bit_list =[1 + 2*l^2, 2 + 2*l^2];

sub_x = squeeze(state(1, :, bit_list));
sub_z = squeeze(state(2, :, bit_list));
sub_mat = logical([sub_x, sub_z]);

output = bitrank_mex(sub_mat) - 2;

end




%% generate stabilizer tableau of the toric code state
function output = tc_state_with_ac(l)


tc_bit_num = 2 * l * l;
full_bit_num = tc_bit_num + 2;

output = zeros(2, full_bit_num, full_bit_num);

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
            loop_p = 1: l;
            loop_v = (1 + l): (2 * l): (2 * l^2);

            
            
            %ZpZ1
            output(2, plq_idx, [loop_p, tc_bit_num + 1]) = 1;
            %ZvZ2
            output(2, tc_bit_num + 1, [loop_v, tc_bit_num + 2]) = 1;


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
            
            dual_loop_v = 1: (2 * l): (2*l^2 - 2*l + 1);
            dual_loop_p = (1 : l) + l;
            
            %XvX1
            output(1, vtx_idx, [dual_loop_v, tc_bit_num + 1]) = 1;
            %XpX2
            output(1, tc_bit_num + 2, [dual_loop_p, tc_bit_num + 2]) = 1;
        else
            output(1, vtx_idx, v_list) = 1;
        end
    end
end



end


%% 4-body-random cliff_uni as noise
function output = uni_noise(l, state, pn)
    tc_bit_num = 2 * l * l;
    xs = (1:4) * 2;
    zs = xs - 1;

    [I, J] = ndgrid(1:l, 1:l);
    idx_mask = ~(I == l & J == l) & rand(l) <= pn;
    idx_list = find(idx_mask);
    
    
    for idx = idx_list'
        i = I(idx);
        j = J(idx);

        if i + j == 1
            continue
        end

        p1 = i + (j - 1) * 2 * l;
        p2 = p1 + l;
        p3 = mod(p1 + 2 * l - 1, tc_bit_num) + 1;
        p4 = mod(p2, l) + floor((p2 - 1) / l) * l + 1;
        p_list = [p1, p2, p3, p4];

        

        x_stab = squeeze(state(1, :, p_list));
        z_stab = squeeze(state(2, :, p_list));

        temp = zeros(8, tc_bit_num + 2);
        temp(xs, :) = x_stab(:,:)';
        temp(zs, :) = z_stab(:,:)';

        cliff_idx = randi(numberofsymplectic(4));
        g = double(symplectic(cliff_idx, 4));

        temp = mod(g * temp, 2);

        state(1, :, p_list) = temp(xs, :)';
        state(2, :, p_list) = temp(zs, :)';
    end

    for idx = idx_list'
        i = I(idx);
        j = J(idx);

        if i + j == 0
            continue
        end

        p1 = i + (j - 1) * 2 * l;
        p2 = p1 + l;
        p3 = mod(p1 + 2 * l - 1, tc_bit_num) + 1;
        p4 = mod(p2, l) + floor((p2 - 1) / l) * l + 1;
        p_list = [p1, p2, p3, p4];

        

        x_stab = squeeze(state(1, :, p_list));
        z_stab = squeeze(state(2, :, p_list));

        temp = zeros(8, tc_bit_num + 2);
        temp(xs, :) = x_stab(:,:)';
        temp(zs, :) = z_stab(:,:)';

        cliff_idx = randi(numberofsymplectic(4));
        g = double(symplectic(cliff_idx, 4));

        temp = mod(g * temp, 2);

        state(1, :, p_list) = temp(xs, :)';
        state(2, :, p_list) = temp(zs, :)';
    end

    output = state;
end

function output = check_measure(l, state)
    tc_bit_num = 2 * l * l;

    % Convert to uint8 for faster XOR operations
    output = uint8(state);

    for i = 1:l
        for j = 1:l
            if i == l && j == l
                continue
            end

            % For Z-plaq
            p1 = i + (j - 1) * 2 * l;
            p2 = p1 + l;
            p3 = mod(p1 + 2 * l - 1, tc_bit_num) + 1;
            p4 = mod(p2, l) + floor((p2 - 1) / l) * l + 1;
            p_list = [p1, p2, p3, p4];

            sigs = mod(sum(output(1, :, p_list), 3), 2);
            sig_ids = find(sigs);

            if ~isempty(sig_ids)
                piv = sig_ids(1);
                piv_slice = output(:, piv, :);

                if numel(sig_ids) > 1
                    all_ids = sig_ids(2:end);
                    % Replicate the pivot slice to match the number of sig_ids
                    rep_piv = repmat(piv_slice, [1, numel(all_ids), 1]);
                    output(:, all_ids, :) = bitxor(output(:, all_ids, :), rep_piv);
                end

                output(:, piv, :) = 0;
                output(2, piv, p_list) = 1;
            end

            % For X-vertex
            v1 = i + (j - 1) * 2 * l;
            v2 = v1 + l;
            v3 = v1 - 1 + l * (v1 - 1 < (j - 1) * 2 * l + 1);
            v4 = v1 - l + tc_bit_num * (v1 - l <= 0);
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

%% The k-bit random_cliffords.

function output = directsum(m1, m2)
% Computes the direct sum of two matrices
[n1, n1_col] = size(m1);
[n2, n2_col] = size(m2);
if n1 ~= n1_col || n2 ~= n2_col
    error('Matrices must be square');
end
output = zeros(n1 + n2, n1 + n2, 'int8');
output(1:n1, 1:n1) = m1;
output(n1+1:n1+n2, n1+1:n1+n2) = m2;
end

% function t = inner(v, w)
% % Symplectic inner product
% n = length(v);
% if mod(n, 2) ~= 0
%     error('Vectors must have even length');
% end
% indices = 1:2:n-1; % indices 1,3,5,...
% t = sum(v(indices) .* w(indices + 1) + w(indices) .* v(indices + 1));
% t = mod(t, 2);
% end

function vout = transvection(k, v)
% Applies transvection Z_k to v
vout = mod(v + inner(k, v) * k, 2);
end

function output = int2bits(i, n)
% Converts integer i to a length n array of bits (LSB first)
output = zeros(1, n, 'int8'); % row vector
for j = 1:n
    output(j) = bitget(i, j);
end
end

function output = findtransvection(x, y)
% Finds h1,h2 such that y = Z_h1 Z_h2 x
n = length(x);
output = zeros(2, n, 'int8');
if isequal(x, y)
    return;
end
if inner(x, y) == 1
    output(1, :) = mod(x + y, 2);
    return;
end
% Find a pair where they are both not 00
z = zeros(1, n, 'int8');
for i = 1:(n/2)
    ii = 2*(i-1) + 1; % MATLAB indices
    if ((x(ii) + x(ii+1)) ~= 0) && ((y(ii) + y(ii+1)) ~= 0)
        % Found the pair
        z(ii) = mod(x(ii) + y(ii), 2);
        z(ii+1) = mod(x(ii+1) + y(ii+1), 2);
        if (z(ii) + z(ii+1)) == 0 % They were the same so they added to 00
            z(ii+1) = 1;
            if x(ii) ~= x(ii+1)
                z(ii) = 1;
            end
        end


        output(1, :) = mod(x + z, 2);
        output(2, :) = mod(y + z, 2);
        return;
    end
end
% Didn't find a pair, so look for specific cases
% First y==00 and x doesn't
for i = 1:(n/2)
    ii = 2*(i-1) + 1;
    if ((x(ii) + x(ii+1)) ~= 0) && ((y(ii) + y(ii+1)) == 0)
        % Found the pair
        if x(ii) == x(ii+1)
            z(ii+1) = 1;
        else
            z(ii+1) = x(ii);
            z(ii) = x(ii+1);
        end
        break;
    end
end
% Finally x==00 and y doesn't
for i = 1:(n/2)
    ii = 2*(i-1) + 1;
    if ((x(ii) + x(ii+1)) == 0) && ((y(ii) + y(ii+1)) ~= 0)
        % Found the pair
        if y(ii) == y(ii+1)
            z(ii+1) = 1;
        else
            z(ii+1) = y(ii);
            z(ii) = y(ii+1);
        end
        break;
    end
end


output(1, :) = mod(x + z, 2);
output(2, :) = mod(y + z, 2);
end

function g = symplectic(i, n)
% Output symplectic canonical matrix i of size 2n x 2n
nn = 2 * n;   % This is convenient to have
% Step 1
s = 2^nn - 1;
k = mod(i, s) + 1;
i = floor(i / s);
% Step 2
f1 = int2bits(k, nn);
% Step 3
e1 = zeros(1, nn, 'int8'); % Define first basis vector
e1(1) = 1;
T = findtransvection(e1, f1); % Use Lemma 2 to compute T
% Step 4
bits = int2bits(mod(i, 2^(nn-1)), nn - 1);
% Step 5
eprime = e1;
eprime(2:nn) = bits(1:nn - 1);
h0 = transvection(T(1, :), eprime);
h0 = transvection(T(2, :), h0);
% Step 6
if bits(1) == 1
    f1 = zeros(1, nn, 'int8'); % f1 *= 0
end
% Step 7
% Define the 2x2 identity matrix
id2 = zeros(2, 2, 'int8');
id2(1,1) = 1;
id2(2,2) = 1;
if n ~= 1
    g_rest = symplectic(floor(i / 2^(nn - 1)), n - 1);
    g = directsum(id2, g_rest);
else
    g = id2;
end
for j = 1:nn
    g(j, :) = transvection(T(1, :), g(j, :));
    g(j, :) = transvection(T(2, :), g(j, :));
    g(j, :) = transvection(h0, g(j, :));
    g(j, :) = transvection(f1, g(j, :));
end
end

function output = bits2int(b, nn)
% Converts an nn-bit string b to an integer between 0 and 2^n - 1
output = 0;
for j = 1:nn
    output = output + b(j) * 2^(j - 1);
end
end

function x = numberofcosets(n)
% Returns the number of different cosets
x = 2^(2*n - 1) * (2^(2*n) - 1);
end

function x = numberofsymplectic(n)
% Returns the number of symplectic group elements
x = 1;
for j = 1:n
    x = x * numberofcosets(j);
end
end

function gnidx = symplecticinverse(n, gn)

% Produce an index associated with group element gn
nn = 2 * n;   % This is convenient to have
% Step 1
v = gn(1, :);
w = gn(2, :);
% Step 2
e1 = zeros(1, nn, 'int8'); % Define first basis vector
e1(1) = 1;
T = findtransvection(v, e1); % Use Lemma 2 to compute T
% Step 3
tw = w;
tw = transvection(T(1, :), tw);
tw = transvection(T(2, :), tw);
b = tw(1);
h0 = zeros(1, nn, 'int8');
h0(1) = 1;
h0(2) = 0;
for j = 3:nn
    h0(j) = tw(j);
end
% Step 4
bb = zeros(1, nn - 1, 'int8');
bb(1) = b;
for j = 3:nn
    bb(j - 1) = tw(j);
end
zv = bits2int(v, nn) - 1; % Number between 0...2^(2n)-2
zw = bits2int(bb, nn - 1);  % Number between 0..2^(2n-1)-1
cvw = zw * (2^(2*n) - 1) + zv;
% cvw is a number indexing the unique coset specified by (v,w)
% Step 5
if n == 1
    gnidx = cvw;
    return;
end
% Step 6
gprime = gn;
if b == 0
    for j = 1:nn
        gprime(j, :) = transvection(T(2, :), transvection(T(1, :), gn(j, :)));
        gprime(j, :) = transvection(h0, gprime(j, :));
        gprime(j, :) = transvection(e1, gprime(j, :));
    end
else
    for j = 1:nn
        gprime(j, :) = transvection(T(2, :), transvection(T(1, :), gn(j, :)));
        gprime(j, :) = transvection(h0, gprime(j, :));
    end
end
% Step 7
gnew = gprime(3:nn, 3:nn); % Take submatrix
gnidx = symplecticinverse(n - 1, gnew) * numberofcosets(n) + cvw;
end


%% verify sympmat
function symp_vf(n)
    n_sym = numberofsymplectic(n);

    symp_metric = logical(kron(eye(n), [0, 1; 1, 0]));

    for i = 1:n_sym
        g = logical(symplectic(i, n));

        % Correct the symplectic condition
        lhs = symp_metric;
        rhs = mod(g' * symp_metric * g, 2);

        if isequal(lhs, rhs)
            disp(i / n_sym);
            continue;
        else
            disp('gg');
            break;
        end
    end
end

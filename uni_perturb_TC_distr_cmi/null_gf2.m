function N = null_gf2(H)
    H = logical(mod(H, 2)); % Convert to logical for GF(2) operations
    [~, n] = size(H);

    % Perform RREF over GF(2)
    [RREF, pivots] = gf2_rref(H);

    free_vars = setdiff(1:n, pivots);
    num_free_vars = length(free_vars);
    N = false(n, num_free_vars);

    % Back-substitution for each free variable
    for i = 1:num_free_vars
        v = false(n, 1);
        v(free_vars(i)) = true;

        % Back-substitution
        for j = length(pivots):-1:1
            col = pivots(j);
            v(col) = mod(sum(RREF(j, :) & v'), 2);
        end
        N(:, i) = v;
    end
end

function [RREF, pivots] = gf2_rref(H)
    H = logical(mod(H, 2)); % Convert to logical for GF(2) operations
    [m, n] = size(H);
    RREF = H;
    pivots = [];

    row = 1;
    for col = 1:n
        % Find pivot in column col at or below row
        pivot_row = find(RREF(row:end, col), 1) + row - 1;
        if ~isempty(pivot_row)
            % Swap rows if necessary
            if pivot_row ~= row
                RREF([row, pivot_row], :) = RREF([pivot_row, row], :);
            end
            pivots(end + 1) = col;

            % Eliminate other entries in the column
            idx = find(RREF(:, col));
            idx(idx == row) = [];
            RREF(idx, :) = xor(RREF(idx, :), RREF(row, :));

            row = row + 1;
            if row > m
                break;
            end
        end
    end
end


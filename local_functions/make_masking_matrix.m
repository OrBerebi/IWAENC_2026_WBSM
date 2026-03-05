function W = make_masking_matrix(omega_V, omega_filtered,tol)
    % omega_V: [directions x 3]
    % omega_filtered: [des_directions x 3]
    % Compare only elevation and azimuth (1st and 2nd cols), rounded to 4 decimals

    % Extract and round elevation and azimuth
    omega_V_2D = round(omega_V(:, 1:2), 4);                % [directions x 2]
    omega_filtered_2D = round(omega_filtered(:, 1:2), 4);  % [des_directions x 2]

    % Initialize mask with low value
    directions = size(omega_V_2D, 1);
    W = ones(directions, 1) * tol;
    %W = zeros(directions, 1);

    % Find matches based on rounded [elevation, azimuth]
    for i = 1:size(omega_filtered_2D, 1)
        match = ismember(omega_V_2D, omega_filtered_2D(i, :), 'rows');
        W(match) = 1;
    end

    % Check how many matches were found
    num_ones = sum(W == 1);
    expected = size(omega_filtered_2D, 1);

    if num_ones == expected
        fprintf('[PASS] Found all %d matching directions.\n', num_ones);
    else
        fprintf('[FAIL] Only %d out of %d directions matched.\n', num_ones, expected);
    end

end

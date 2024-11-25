function [b, a] = lpc2filter(LP_coeffs)
    % Ensure LP_coeffs is a row vector
    if iscolumn(LP_coeffs)
        LP_coeffs = LP_coeffs';  % Convert to row vector if it's a column vector
    end
    
    % Ensure that LP_coeffs has at least 2 coefficients (p + 1 coefficients)
    if length(LP_coeffs) < 2
        error('LP coefficients must have at least 2 elements');
    end
    
    % The first coefficient is 1, followed by the negated LP coefficients (a(2:end))
    b = [1, -LP_coeffs(2:end)];
    a = 1;  % The denominator is just 1 for the LP filter
end

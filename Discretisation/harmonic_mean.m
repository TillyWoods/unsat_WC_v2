function hm = harmonic_mean(b)
    % 26th october 2022
    % Approximate values on the cell edges by the harmonic mean of the
    % cell centre values. b is a N column vector of cell centre values. hm
    % is a N+1 column vector of the harmonic means at the cell edges.
    b1 = [b;b(end)];
    b2 = [b(1);b];
    hm = 2*b1.*b2./(b1+b2);         % harmonic mean
    hm(or(b1==0,b2==0)) = 0;        % the harmonic mean should be zero if both arguments are zero (but the above function gives NaN)
end
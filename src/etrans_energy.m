function r=etrans_energy(eq,s)
%ETRANS_ENERGY Calculates energy terms averged over samples
%
%   R = ETRANS_ENERGY(EQ, S) calculates energy terms averged over samples.
%
%   EQ is a structure with parameters:
%       EQ.omegaM2 - square of classical frequence
%       EQ.chi_c - coupled constant in classical equation
%       EQ.chi_q - coupled constant in quantum equation
%       EQ.d(N,1) - diagonal of quantum matrix
%       EQ.s(N-1,1) - subdiagonal of quantum matrix
%
%   This structure is stored in etrans output file
%
%   S(N,4,M) is array of M samples all phase variables for N sites
%
%   RETURN:
%   R - structure with energy terms:
%       R.MEv - mean kinetic energy of classical sites
%       R.SEv - mean square kinetic energy of classical sites
%       R.MEu - mean potential energy of classical sites
%       R.SEu - mean square potential energy of classical sites
%       R.MEq - mean charge energy
%       R.SEq - mean square charge energy
%       R.MEb - mean charge-chain interaction energy
%       R.SEb - mean square charge-chain interaction energy
[n, ~, nsamp] = size(s);

if numel(eq.d) ~= n
    error('Number of sites does not match');
end

x = reshape(s(:,1,:), [n, nsamp]);
y = reshape(s(:,2,:), [n, nsamp]);
u = reshape(s(:,3,:), [n, nsamp]);
v = reshape(s(:,4,:), [n, nsamp]);

a = 0.5 * eq.omegaM2 * u.*u;
r.SEu = mean(a.*a,2);
r.MEu = mean(a,2);

a = 0.5 * v.*v;
r.SEv = mean(a.*a,2);
r.MEv = mean(a,2);

b2 = (x.*x + y.*y);
a = (eq.chi_q * eq.chi_c / eq.omegaM2) * u .* b2 .* b2;
r.SEb = mean(a.*a,2);
r.MEb = mean(a,2);


a = (eq.s * ones(1,nsamp)) .* (x(1:end-1,:).*x(2:end,:)) + y(1:end-1,:).*y(2:end,:);
a = (eq.d * ones(1,nsamp)) .* b2 + [zeros(1,nsamp); a] + [a; zeros(1,nsamp)];

r.SEq = mean(a.*a,2);
r.MEq = mean(a,2);
end
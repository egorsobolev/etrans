function x=etrans_read_fp(fn, nsamp)
%ETRANS_READ_FP Read etrans final point file
%
%   X = ETRANS_READ_FP(FN) reads FN as final point file of etrans.
%
%   X = ETRANS_READ_FP(FN, N) reads N samples from final point file FN
%   of etrans program.
%
%   RETURN:
%   X(N,4,M) array of M samples all phase variables for N sites
if nargin < 1
  error('Specify filename');
end

f = fopen(fn, 'r');
r.n = fread(f, 1, 'int');
r.n0 = fread(f, 1, 'int');

if nargin > 1
    x = fread(f, 4 * r.n * nsamp, 'double');
else
    x = fread(f, inf, 'double');
    nsamp = numel(x) / (4 * r.n);
end
    x = reshape(x, [r.n, 4, nsamp]);

fclose(f);

end

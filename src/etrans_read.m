function r=etrans_read(fn)
%ETRANS_READ Read etrans output file
%
%   R = ETRANS_READ(FN) reads FN as output file of etrans. Etrans is
%   a program for integration of the semiclassical model of charge transfer
%   in DNA.
%
%   R = 
%           n: length of chain
%       nstep: number of large steps with output
%       nskip: number of simple steps without output
%           h: time of one step
%       chain: chain type:
%               OSC - oscillators (linear model)
%               PBD - Peyrard-Bishop-Dauxois model
%          eq: equation parameters (depens on chain type)
%       nsamp: number of averaged samples
%       nfunc: number of functionals stored in output file
%     M<name>: average over samples of functional <name>
%    hM<name>: corresponding time step
%    nM<name>: corresponding last trajectory point
%           ...
%     S<name>: average over sumples sum of squares of functional <name>
%    hS<name>: corresponding time step
%    nS<name>: corresponding last trajectory point
%           ...
%         cmd: command used for run etrans
%
%   If R.chain == 'OSC', then R.eq is
%
%   R.eq =
%          kt: Temperature
%     omegaB2: Elastic
%     omegaM2: Dispersion
%       Gamma: Friction
%       chi_c: Coupled constant in classical equation
%      sigmaF: Standard deviation of random force
%          n0: Origin
%       chi_q: Coupled constant in quantum equation
%      lambda: Electric field
%           d: Diagonal of quantum matrix
%           s: Side element of quantum matrix
%
%   If R.chain == 'PBD', then R.eq is
%
%   R.eq =
%          ...
%          n0: Origin
%       chi_q: Coupled constant in quantum equation
%      lambda: Electric field
%           d: Diagonal of quantum matrix
%           s: Side element of quantum matrix
%
%   Example
%      r=etrans_read('g301_300K.dat');

%   Copyright 2015  IMPB RAS.

if nargin < 1
  error('Specify filename');
end

f = fopen(fn, 'r');
r.n = fread(f, 1, 'int');
r.nstep = fread(f, 1, 'int64');
r.nskip = fread(f, 1, 'int');
r.h = fread(f, 1, 'double');

r.chain = fread(f, 3, '*char')';
if strcmp(r.chain, 'OSC')
    p = fread(f, 6, 'double');
    r.eq.kt = p(1);
    r.eq.omegaB2 = p(2);
    r.eq.omegaM2 = p(3);
    r.eq.Gamma = p(4);
    r.eq.chi_c = p(5);
    r.eq.sigmaF = p(6);
elseif strcmp(r.chain, 'PBD')
    p = fread(f, 10, 'double');
    r.eq.omegaM2 = p(1);
    r.eq.sigma = p(2);
    r.eq.omegaB2 = p(3);
    r.eq.epsilon = p(4);
    r.eq.rho = p(5);
    r.eq.theta0 = p(6);
    r.eq.theta = p(7);
    r.eq.sigmaF = p(8);
    r.eq.Gamma = p(9);
    r.eq.chi_c = p(10);
else
    error('Unknown chain model');
end

r.eq.n0 = fread(f, 1, 'int');
r.eq.chi_q = fread(f, 1, 'double');
r.eq.lambda = fread(f, 1, 'double');
r.eq.d = fread(f, r.n, 'double');
r.eq.s = fread(f, r.n-1, 'double');

r.nsamp = fread(f, 1, 'int');
r.nfunc = fread(f, 1, 'int');

name = cell(r.nfunc,1);
flag = zeros(r.nfunc,1);
m = zeros(r.nfunc,1);
s = zeros(r.nfunc,1);
for i = 1:r.nfunc
    len = fread(f, 1, 'unsigned char');
    name{i} = fread(f, len, '*char')';
    flag(i) = fread(f, 1, 'unsigned char');
    if bitand(flag(i), 1)
        m(i) = fread(f, 1, 'int64');
    end
    if bitand(flag(i), 2)
        s(i) = fread(f, 1, 'int64');
    end
end

n = fread(f, r.nfunc, 'int');

for i = 1:r.nfunc
    if bitand(flag(i), 1)
        ni = r.nstep / m(i);
        r.(['M' name{i}]) = fread(f, [n(i) ni+1], 'double');
        r.(['hM' name{i}]) = r.h * m(i) * r.nskip;
        r.(['nM' name{i}]) = ni;
    end
    if bitand(flag(i), 2)
        ni = r.nstep / s(i);
        r.(['S' name{i}]) = fread(f, [n(i) ni+1], 'double');
        r.(['hS' name{i}]) = r.h * s(i) * r.nskip;
        r.(['nS' name{i}]) = ni;
    end
end

argc = fread(f, 1, 'int');
n = fread(f, 1, 'int');
r.cmd = char(fread(f, n, 'char')');
for i=2:argc
  n = fread(f, 1, 'int');
  r.cmd = [r.cmd ' ' char(fread(f, n, 'char')')];
end
fclose(f);

return

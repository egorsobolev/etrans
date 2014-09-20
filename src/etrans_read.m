function r=etrans_read(fn)

if nargin < 1
  error('Specify filename');
end

f = fopen(fn, 'r');
r.n = fread(f, 1, 'int');
r.nsamp = fread(f, 1, 'int');
r.nstep = fread(f, 1, 'int');
r.h = fread(f, 1, 'single');
r.pos = fread(f, 1, 'int');
r.uos = fread(f, 1, 'int');

if r.pos
    r.o_nstep = r.nstep / r.uos;
else
    r.o_nstep = 0;
end
r.o_h = r.h * r.uos;
if r.uos
    r.q_nstep = r.nstep / r.pos;
else
    r.q_nstep = 0;
end
r.q_h = r.h * r.pos;

r.p = fread(f, [r.n r.q_nstep+1], 'double');
r.c = fread(f, [r.nstep+1 3], 'double');
r.u = fread(f, [r.n r.o_nstep+1], 'double');
r.p2 = fread(f, [r.n r.q_nstep+1], 'double');
r.c2 = fread(f, [r.nstep+1 3], 'double');
r.u2 = fread(f, [r.n r.o_nstep+1], 'double');

argc = fread(f, 1, 'int');
n = fread(f, 1, 'int');
r.cmd = char(fread(f, n, 'char')');
for i=2:argc
  n = fread(f, 1, 'int');
  r.cmd = [r.cmd ' ' char(fread(f, n, 'char')')];
end
fclose(f);

return

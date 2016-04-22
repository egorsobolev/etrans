function etrans_plot(s, f)

if nargin < 2
    f = 'MEq';
    g = @(a) a.MEq+a.MEb+a.MEu+a.MEv;
else
    g = @(a) a.(f);
end

t0 = 0;
hold all
for a = s
    plot(t0 + (0:a.(['n' f])) * a.(['h' f]), g(a));
    t0 = t0 + a.(['n' f])*a.(['h' f]);
end
hold off

end

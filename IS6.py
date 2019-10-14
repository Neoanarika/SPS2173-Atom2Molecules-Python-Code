import sympy as sp

A, x, L = sp.symbols('A, x, L')
p = A*x*(L-x)
p.subs(A, 112).subs(L, 10)

sp.plot(p.subs(A, 1).subs(L, 1), (x, 0, 1), line_color='g')
sp.integrate(p, x)
rst = sp.integrate(p, (x, 0, L))
rst

M = sp.symbols('M')
sp.solve(rst-M, A)

x, e, t = sp.symbols('x, e, t', real=True)
h_bar, L, m = sp.symbols('\\hbar, L, m', real=True, positive=True)
n = sp.symbols('n', integer=True, positive=True)
En = (n*sp.pi*h_bar/L)**2/(2*m)
haha = sp.sqrt(2/L)*sp.sin(((n*sp.pi)/L)*x)
wowo = sp.exp(-sp.I*En*t/h_bar)
p = abs(wowo*haha)**2
p
sp.integrate(p, (x, 0, L))

for i in range(1, 10):
    sp.plot(p.subs(L, 1).subs(n, i).subs(t, 1).subs(m, 1), (x, 0, 1))

sp.re(wowo*haha).subs(L, 1).subs(n, 1).subs(t, 1).subs(m, 1)
for i in range(1, 100):
    sp.plot(sp.re(wowo*haha).subs(L, 1).subs(n, 1).subs(t, 0.01*i).subs(m, 1).subs(h_bar, 1), (x, 0, 1), ylim=(-2, 2))

lower = 1/8*L
upper = 1/4*L

for i in range(1, 6):
    rst = sp.integrate(p.subs(L, 1).subs(n, i).subs(t, 1).subs(m, 1),
                       (x,
                       lower.subs(L, 1),
                       upper.subs(L, 1)))
    prob = sp.N(rst)
    print(f"{prob}")

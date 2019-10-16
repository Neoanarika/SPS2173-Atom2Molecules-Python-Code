import sympy as sp

'''
Expectated value  we expect it to be at the peak in a symmetric distribution
KIV : expectation vs mode

Wavefunctions as a superposition of other wavefunctions in another basis
Example when in particle in a box the basis of energy in particle box repersents the state
, but we are measuring something like ppsition we change to the position basis
when we do the measurment. Key idea: Measurment changes basis
[thats where the orthonormality comes in]

A operator is a "recipe" not really like a function

Quantum mechanics
'''

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
    sp.plot(sp.re(wowo*haha).subs(L, 1).subs(n, 3).subs(t, 0.01*i).subs(m, 1).subs(h_bar, 1), (x, 0, 1), ylim=(-2, 2))


x_hat = lambda y: x*y
x_rslt = sp.integrate(sp.conjugate(wowo*haha)*x_hat(wowo*haha), (x, 0, L))
sp.simplify(x_rslt.subs(n, 2))
xx_rslt = sp.integrate(sp.conjugate(wowo*haha)*x_hat(x_hat(wowo*haha)), (x, 0, L))

'''
Its L/2 for everything because the graph is symmetric
'''
for i in range(1, 10):
    print(sp.simplify(rslt.subs(n, i)))

p_hat = lambda y: (-sp.I * h_bar)*sp.diff(y, x)
wavefunction = wowo*haha
p_rslt = sp.integrate(sp.conjugate(wavefunction)*p_hat(wavefunction).subs(m, 1), (x, 0, L))
sp.simplify(p_rslt.subs(n, 1))
pp_rslt = sp.integrate((sp.conjugate(wavefunction)*p_hat(p_hat(wavefunction))).subs(m, 1), (x, 0, L))
pp_rslt
T_hat = lambda y: (-(h_bar**2)/(2*m))*sp.diff(sp.diff(y, x), x)
wavefunction = (wowo*haha).subs(h_bar, 1)
rslt = sp.integrate(sp.conjugate(wavefunction)*T_hat(wavefunction), (x, 0, L))
sp.simplify(rslt)
sp.simplify(rslt.subs(n, 1))

varience_of_x = sp.simplify(sp.sqrt(xx_rslt - x_rslt**2))
varience_of_p = sp.simplify(sp.sqrt(pp_rslt - p_rslt**2))
sp.simplify(varience_of_x*varience_of_p.subs(n, 1))

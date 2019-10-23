import sympy as sp
import matplotlib.pyplot as plt

'''
We need to put real and positive if not the expression below will become
very complex
'''
x = sp.symbols('x', real=True)
L = sp.symbols('L', real=True, positive=True)
'''
I want to check the boundary conditions and whether the function is
normalised
'''
f = sp.sqrt(30/(L**5))*x*(L-x)

'''
You do not need to add the abs, when u squared it. In reality it is not
a absolute function. Furthermore, sympy requires you to use sp.Abs not abs function
unsuprisingly
'''
sp.integrate(f**2, (x, 0, L))
sp.plot((f**2).subs(L, 1),(x, 0, 1))
plt.show()

'''
Expectation of (E, x, p)

Expectation of a quantity = integral(conjugate(wavefunction)*operator(wavefunction))
'''

h_bar, L, m = sp.symbols('\\hbar, L, m', real=True, positive=True)

'''
Not suprising from the graph above, since the probability distribution is
symmetric so it will on average most likely be ard L/2. Take not the peakness
argument will not work for n=2 onwards !!
'''
x_hat = lambda y: x*y
expectation_x = sp.simplify(sp.integrate(sp.conjugate(f)*x_hat(f), (x, 0, L)))
sp.simplify(expectation_x)


'''
Expectation of p is also not suprising because of the symmetry of the box,
since the particle is moving left or right with equal probability so the total
probability is ard 0. Interestingly enough, even tho momenumtum does not follow
Classical the expectation of momentum does !!
'''
p_hat = lambda y: (-sp.I * h_bar)*sp.diff(y, x)
expectation_p = sp.integrate(sp.conjugate(f)*p_hat(f).subs(m, 1), (x, 0, L))
sp.simplify(expectation_p)

T_hat = lambda y: (-(h_bar**2)/(2*m))*sp.diff(sp.diff(y, x), x)
expectation_T = sp.integrate(sp.conjugate(f)*T_hat(f), (x, 0, L))
sp.simplify(expectation_T)

'''
T|S> = E|S>, the energy is just the eigenvalue.

Solving schrodinger equation is basically solving the eigenvalue of applying
the energy operator on the wavefunction (eigenfunction) !!
'''

'''
The wavefunction are eigen-function of a particular operator, so you can get
multiple eigen-functions which is orthogonal to each other. Why because?
You do expect to be in state 2 when you know you are in state 1.

Any function can be decompose into a linar combination of the eigenfunction
i.e. any wavefunction is an element of the span of the eigenfunctions.

In otherwords, the eigenfunctions form a basis for a vector space that
contains all the valid wavefunctions.

The dot product here is with the intergal !!
'''

n = sp.symbols('n', real=True, positive=True)
psi_n = sp.sqrt(2/L)*sp.sin(((n*sp.pi)/L)*x)
c_n = sp.integrate(sp.conjugate(psi_n)*f, (x,0,L)) # Dot product
for i in range(1, 6):
    print(sp.N(c_n.subs(n, i)))

'''
the reason why all the even values of n are 0 is because the area of
graph negtaive and positive cancels out when u integrate. 
'''

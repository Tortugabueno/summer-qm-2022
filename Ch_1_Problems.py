import numpy as np
import sympy as sy
from sympy.plotting import plot

oo = sy.oo

pi = sy.pi

i = sy.I

# ======== 1.1(a) ======== (tortugabueno)
print('')
print('1.1(a)')
print('')

def mean(list):
  return sum(list)/len(list)

ages = [14, 15, 16, 16, 16, 22, 22, 24, 24, 25, 25, 25, 25, 25]
n = len(ages)

def square(list):
  temp_list = list.copy()
  for i in range(len(list)):
    temp_list[i] = temp_list[i]**2
  return temp_list

print('< j**2 > = ' + str(mean(square(ages))) + ' years.')
print('')
print('< j >**2 = ' + str(mean(ages)**2) + ' years.')


# ======== 1.1(b) ======== (tortugabueno)
print('')
print('1.1(b)')
print('')

def deviations(list):
  m = mean(list)
  temp_list = list.copy()
  for i in range(len(list)):
    temp_list[i] = list[i] - m
  return temp_list

for dev in deviations(ages):
  print(dev)

def variance(list):
  return mean(square(deviations(ages)))

def stddev(list):
  return variance(ages)**(1/2)

print('')
print('The standard deviation is ' + str(stddev(ages)) + ' years.' )
print('')

print('1.1(c)')
print('')

sigma = (mean(square(ages)) - mean(ages)**2)**(1/2)
  
print('Using Eqn 1.12, the standard deviation is ' + str(sigma) + ' years.')
print('')


# ======== 1.2(a) ======== (tortugabueno)
print('')
print('1.2(a)')
print('')

h = sy.symbols('h', positive = True)

def rho(var):
  return (1/2)*(h*var)**(-1/2)

def function_expectation(func, density, var):
  return sy.integrate(func(var)*density(var),var)

def id(var):
  return var

def sq(var):
  return var**2

x = sy.symbols('x', real=True)
p1 = plot(x*x, show=False)
p2 = plot(x, show=False)

def density_mean(density, var):
  return function_expectation(id, density, var)

def density_variance(density, var):
  return function_expectation(sq, density, var) - density_mean(density, var)**2

F = density_variance(rho, x)

mean_rho = sy.nsimplify(density_mean(rho,x).subs(x,h) - density_mean(rho,x).subs(x,0))

var_rho = sy.nsimplify(F.subs(x,h) - F.subs(x,0))

stddev_rho = sy.simplify(var_rho**(1/2))

print('The variance of the distribution with density given by rho is ' + str(sy.nsimplify(F.subs(x,h) - F.subs(x,0))) + ' so the standard deviation is ' + str(stddev_rho) + '.')

# ======== 1.2(b) ======== (tortugabueno)
print('')
print('1.2(b)')
print('')
print('We want to find the probability that a randomly chosen photograph will show a position x more than 1 sd from the average. Since the mean of the distribution is mean_rho = ' + str(mean_rho) +' and standard deviation of the distribution is stddev_rho = ' + str(stddev_rho) + ', we seek the probability that a photograph shows the position less than ' + str(sy.simplify(mean_rho - .298142396999972*h)) + ' or greater than ' + str(sy.simplify(mean_rho + .298142396999972*h)) + '. This is the complement of the probability that the photograph shows the position between ' + str(sy.simplify(mean_rho - .298142396999972*h)) + ' and ' + str(sy.simplify(mean_rho + .298142396999972*h)) + '.') 
print('')

a = mean_rho - .298142396999972*h
b = mean_rho + .298142396999972*h

t = sy.symbols('t')

P = sy.integrate(rho(t), (t, a, b))

print('Therefore, the probability we seek is ' + str(sy.simplify(sy.expand(1 - P))) + '.')
print('')

# ======== 1.3(a) ======== (tortugabueno)
print('1.3(a)')
print('')

A = sy.symbols('A', )
l = sy.symbols('l', positive = True)
a = sy.symbols('a')
b = sy.symbols('b')

def rho(var):
  return A*sy.exp(-l*(var - a)**2)

integral_rho = sy.simplify(sy.integrate(rho(x),(x,-oo,oo)))

A_sol = sy.solve(sy.simplify(sy.integrate(rho(x),(x,-oo,oo))) - 1,A)

print(
  'Since integrating rho(x) over (-oo, oo) gives ' + str(integral_rho) + ' = 1, it follows that A = ' + str(A_sol) + '.' )

# ======== 1.3(b) ======== (tortugabueno)
print('')
print('1.3(b)')
print('')

rho = sy.sqrt(l/pi)*sy.exp(-l*(x - a)**2)
  
mean_rho = sy. simplify(sy.integrate( x*rho , (x, -oo, oo)))

print('< x > = ' + str(mean_rho))

sq_rho = sy.simplify(sy.integrate( x**2*rho, (x, -oo, oo)))

print('< x**2 > = ' + str(sq_rho))

variance_rho = sy.simplify(sq_rho - mean_rho**2)

stddev_rho = sy.simplify(variance_rho**(1/2))

print('stddev_rho = ' + str(stddev_rho) + '.')

# ======== 1.3(b) ======== (tortugabueno)
print('')
print('1.3(c)')
print('')
  
print('Figure 1.3(c) is a plot of rho(x) with lambda = 1, a = 0 ')

p = plot(rho.subs(l,1).subs(a,0), (x,-5,5), title = 'Figure 1.3(b)')

del globals()['A']
  
# ======== 1.4(a) ======== (tortugabueno)
print('')
print('1.4(a)')
print('')

A = sy.symbols('A')

Psi_1 = A * x / a
Psi_2 = A * (b - x) / (b - a)

S = sy.simplify(sy.integrate(Psi_1**2, (x, 0, a)) + sy.integrate(Psi_2**2, (x, a, b)))

A_sol = sy.solve(S - 1, A)
print(str(A_sol))
print('')

Psi_1 = sy.sqrt(3) / sy.sqrt(b) * x / a
Psi_2 = sy.sqrt(3) / sy.sqrt(b) * (b - x) / (b - a)

print('Psi_1 is ' + str(Psi_1))
print('')

print('Since the integral S of |Psi|**2 over (-oo,oo) is ' + str(S) +
      ', solving S = 1 gives A = ' + str(A_sol) + '.')

# ======== 1.4(b) ======== (tortugabueno)
print('')
print('1.4(b)')
print('')

print('Figure 1.4(b) is a plot of Psi for a = 2, b=3')
print('')

Psi_plot = sy.Piecewise((0, x < 0), (Psi_1.subs(a, 2).subs(b, 3), x < 2), (Psi_2.subs(a, 2).subs(b, 3), x<3), (0,True))

p1 = plot(Psi_plot, (x,0,3), show = False , title = 'Figure 1.4(a)')
p1.show()

# ======== 1.4(c) ======== (tortugabueno)
print('1.4(c) ')
print('')

print('Since the maximum of |Psi(x,0)|**2 occurs at x = a, this is the most likely position of the particle.')

# ======== 1.4(d) ======== (tortugabueno)
print('')
print('1.4(d) ')
print('')

P = sy.integrate(Psi_1**2, (x,0,a))

print('The probability that the particle is observed between 0 and a when t= 0 is the integral of |Psi(x,0)|**2 over (0,a), which is ' + str(P) +'.')
print('')


print('When a = b, the probability that x is in (0, a) is ' + str(P.subs(b,a)) +'.')
print('')

print('When 2a = b, the probability that x is in (0, a) is ' + str(P.subs(b,2*a)) +'.')

print('')
print('1.4(e) ')
print('')

mean = sy.simplify(sy.integrate(x*Psi_1**2, (x, 0, a)) + sy.integrate(x*Psi_2**2, (x, a, b)))

print('< x > = ' + str(mean) + '.')

print('')
print('1.5(a) ')
print('')

A = sy.symbols('A', positive = True, real = True)
t = sy.symbols('t', positive = True, real = True)
w = sy.symbols('w', positive = True, real = True)
k = sy.symbols('k', positive = True, real = True)
i = sy.I

def abs(z):
  if (sy.im(z) == 0):
    return sy.Abs(z)
  return (sy.re(z)**2 + sy.im(z)**2)**(1/2)

Psi = A*sy.exp(-k*abs(x))*sy.exp(-i*w*t)

Psi_0 = Psi.subs(t,0)

print('First, we substitute t = 0 to obtain Psi(x,0) = ' + str(sy.simplify(Psi_0)) + '. Notice that Psi(x,0) is an even real-valued function. It follows that the integral of Psi(x,0) over (-oo,oo) is equal to twice the integral of Psi(x,0) over (0,oo).')
print('')

Psi_pos = A*sy.exp(-k*x)

S = sy.integrate(2*(Psi_pos)**2, (x,0,oo))

A_sol = sy.solve(S-1,A)

print('Then we have that the integral of 2*Psi(x,0)**2 over (-oo,oo) is ' + str(S) + ', so A = ' + str(A_sol) + '.' )
print('')

Psi = sy.sqrt(k)*sy.exp(-k*sy.Abs(x))*sy.exp(-i*w*t)
Psi_0 = Psi.subs(t,0)

print('1.5(b) ')
print('')

S = sy.integrate(x*(Psi_0)**2, (x,-oo,oo))

print('< x > = ' + str(S))
print('')

S_2 = sy.integrate(x**2*(Psi_0)**2, (x,-oo,oo))

print('< x**2 > = ' + str(S_2))

print('')
print('1.5(c) ')
print('')

sigma = sy.sqrt(S_2 - S**2)
print('sigma = ' + str(sigma))
print('')

P = Psi.subs(x,sigma).subs(t,0)**2

print('Psi(sigma,0) = ' + str(P))
print('')

print('Figure 1.5(c) is a plot of Psi for k = 1')
print('')

Psi_1_0 = Psi.subs(k,1).subs(t,0)

x_min = -3.5
x_max = 3.5

L = 0 - sigma.subs(k,1)
R = 0 + sigma.subs(k,1)

plot_1_5_a = plot((Psi_1_0)**2, (x,x_min,x_max), show = False, title = 'Figure 1.5(c)')

def shadebeyond(plot_name, func, var, L, R, x_min, x_max, title_str):
    y_func = sy.Piecewise( (func, var < L), (func, var > R), (0,True) )
    x_array = np.linspace(x_min, x_max, 1000)
    y_array = sy.lambdify(x, y_func)(x_array)
    p1 = plot(0,(var,0,0), fill={'x': x_array,'y1':y_array,'color':'blue'}, show = False, title = title_str)
    p1.extend(plot_name)
    return p1
    
plot_1_5_a = shadebeyond(plot_1_5_a, (Psi_1_0)**2, x, L, R, x_min, x_max, 'Figure 1.5(c)')
plot_1_5_a.show()

def shadebetween(plot_name, func, var, L, R, x_min, x_max):
    y_func = sy.Piecewise( (0, x < L), (func, x < R), (0,True) )
    x_array = np.linspace(x_min, x_max, 1000)
    y_array = sy.lambdify(x, y_func)(x_array)
    p1 = plot(0,(var,0,0), fill={'x': x_array,'y1':y_array,'color':'blue'}, show = False)
    p1.extend(plot_name)
    return p1

S = sy.integrate(2*(Psi_0)**2, (x, sigma, oo)) 

print('The probability that x is observed more than one standard deviation from < x > is ' + str(S) + '.')

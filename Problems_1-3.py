import numpy as np
from sympy import *
from sympy.plotting import plot

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

h = symbols('h', positive = true)

def rho(var):
  return (1/2)*(h*var)**(-1/2)

def function_expectation(func, density, var):
  return integrate(func(var)*density(var),var)

def id(var):
  return var

def sq(var):
  return var**2

x = symbols('x')
p1 = plot(x*x, show=False)
p2 = plot(x, show=False)

def density_mean(density, var):
  return function_expectation(id, density, var)

def density_variance(density, var):
  return function_expectation(sq, density, var) - density_mean(density, var)**2

F = density_variance(rho, x)

mean_rho = nsimplify(density_mean(rho,x).subs(x,h) - density_mean(rho,x).subs(x,0))

var_rho = nsimplify(F.subs(x,h) - F.subs(x,0))

stddev_rho = simplify(var_rho**(1/2))

print('The variance of the distribution with density given by rho is ' + str(nsimplify(F.subs(x,h) - F.subs(x,0))) + ' so the standard deviation is ' + str(stddev_rho) + '.')

# ======== 1.2(b) ======== (tortugabueno)
print('')
print('1.2(b)')
print('')
print('We want to find the probability that a randomly chosen photograph will show a position x more than 1 sd from the average. Since the mean of the distribution is mean_rho = ' + str(mean_rho) +' and standard deviation of the distribution is stddev_rho = ' + str(stddev_rho) + ', we seek the probability that a photograph shows the position less than ' + str(simplify(mean_rho - .298142396999972*h - x)) + ' or greater than ' + str(simplify(mean_rho + .298142396999972*h + x)) + '. This is the complement of the probability that the photograph shows the position between ' + str(simplify(mean_rho - .298142396999972*h - x)) + ' and ' + str(simplify(mean_rho + .298142396999972*h + x)) + '.') 

print('')

a = mean_rho - .298142396999972*h
b = mean_rho + .298142396999972*h

t = symbols('t')

P = integrate(rho(t), (t, a, b))

print('Therefore, the probability we seek is ' + str(simplify(expand(1 - P))) + '.')
print('')

# ======== 1.3(a) ======== (tortugabueno)
print('1.3(a)')
print('')

A = symbols('A', )
l = symbols('l', positive = true)
a = symbols('a')
b = symbols('b')

def rho(var):
  return A*exp(-l*(var-a)**2)

integral_rho = simplify(integrate(rho(x),(x,-oo,oo)))

A = solve(simplify(integrate(rho(x),(x,-oo,oo))) - 1,A)

print(
  'Since integrating rho(x) over (-oo, oo) gives ' + str(integral_rho) + ' = 1, it follows that A = ' + str(A) + '.' )

# ======== 1.3(b) ======== (tortugabueno)
print('')
print('1.3(b)')
print('')

def rho(var):
  return sqrt(l/pi)*exp(-l*(var-a)**2)
  
mean_rho = simplify(integrate( x*rho(x), (x, -oo, oo)))

print('< x > = ' + str(mean_rho))

sq_rho = simplify(integrate( x**2*rho(x), (x, -oo, oo)))

print('< x**2 > = ' + str(sq_rho))

variance_rho = simplify(sq_rho - mean_rho**2)

stddev_rho = simplify(variance_rho**(1/2))

print('stddev_rho = ' + str(stddev_rho) + '.')

# ======== 1.3(b) ======== (tortugabueno)
print('')
print('1.3(c)')
print('')

def rho(var):
  return sqrt(1/pi)*exp(-1*(var-0)**2)
  
print('Below is a plot of rho(x) with lambda = 1, a = 0 ')

p = plot(rho(x))



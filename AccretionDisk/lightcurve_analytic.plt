set xlabel "Time (day)"
set ylabel "Luminosity (erg/s)"

set logscale y
unset logscale x

p = 1
q = 2./3.

G = 6.67259e-8
c = 2.99792458e10
alpha = 0.1
k_b = 1.380658e-16
thompson = 6.6524e-25
m_p = 1.6726231e-24
k_es = thompson / m_p
a = 5.67051e-5;

M_0 = 1e-8 * 1.988e33
M_c = 10 * 1.988e33
J = 1.64488e+44 


k = (q/((4*q - 2*p + 4) * (5 * q - 2*p + 4)))**(1/q)
gamma_1 = 3.6851e-4
gamma_2 = 6.2626e-4

r = 6 * G * M_c / (c**2)
r_0 = (J / M_0)**2 / (G * M_c) * (gamma_2 / gamma_1)**2

E_0 = M_0 / (4*pi*(r_0**2)*gamma_2)

V_C = ((27./32.) * (0.1**4 * k_b**4 * k_es) / ((0.63**4)* (m_p**4) * a * G * M_c))**(1./3.)
V_0 = V_C  * r_0**p * E_0**q
t_0 = 4 * r_0**2 / (3 * V_0)

print V_C
print (19./16. -1)* M_0 / t_0* 0.1 * G * M_c / (2*r)
f(x) =  (19./16. -1)* M_0 / t_0 *(1 + (x * 8.64e4)/t_0)**( -19./16.) * 0.1 * G * M_c / (2*r)

plot[:][1e33:] "lightcurve_analytic_0.033000.txt" with lines, "lightcurve_analytic_0.233000.txt" with lines, "lightcurve_analytic_0.363000.txt" with lines, "lightcurve_analytic_1.023000.txt" with lines
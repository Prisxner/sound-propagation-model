import matplotlib.colors
import mpmath
import numpy as np
import pylab as p
import scipy
from scipy import integrate
import matplotlib.pyplot as plt
import matplotlib.transforms as trns
import sympy as sp
import mpmath as mp
from mpmath import cos, quad, sin, exp, pi, sqrt, tan

T = 7
L = 1
h_1, h_2, h_3, l_3 = 0.3*L, 0.3*L, 0.3*L, 0.3*L
alfa, beta = float(45/180*pi), float(135/180*pi)
k = float(2*pi/L)

if alfa > 0 and alfa < float(pi/2):
    l_1 = float(sin(alfa) * (h_2 + (h_1 / cos(alfa) - h_2)) - (h_1 / cos(alfa) - h_2) / sin(alfa))
    l_2 = float(((h_1 / cos(alfa)) - h_2) / (tan(alfa)))
elif alfa > float(pi/2) and alfa < float(pi):
    l_1 = float(((h_2 / sin(alfa - (pi / 2))) + h_1) / (tan(pi - alfa)))
    l_2 = float((((h_2 / sin(alfa - (pi / 2))) + h_1) / (sin(pi - alfa))) - (h_2 / tan(alfa - (pi / 2))))

if beta > 0 and beta < float(pi/2):
    l_4 = float(((h_2 / sin(pi / 2 - beta)) + h_3) / (tan(beta)))
    l_5 = float(- h_2 * (1 / tan(pi / 2 - beta)) + ((h_2 / sin(pi / 2 - beta)) + h_3) / sin(beta))
elif beta > float(pi/2) and beta < float(pi):
    l_4 = float(((h_2 / cos(pi - beta)) - h_3) / (tan(pi - beta)))
    l_5 = float(h_2 * tan(pi - beta) - ((h_2 / cos(pi - beta)) - h_3) / sin(pi - beta))

print(l_1, l_2, l_5,l_4)
def delta(n):
    if n == 0:
        return 1
    else:
        return 0

def ksi(n,k,h):
    if k > (np.pi*n/h):
        return np.sqrt(k**2 - (np.pi*n/h)**2)
    elif k < (pi*n/h):
        return 1j * np.sqrt((np.pi*n/h)**2 - k**2)

def eps(n):
    if n == 0:
        return 1
    else:
        return 0.5

'''
    Перейдемо до пошуку коефіцієнтів у першому рівнянні системи
'''

Am, Bm, Cm, Dm, Fm, Lm, Km, Em, Gm, Hm, Rm, Xm = ([[0 for i in range(T)] for j in range(T)] for q in range(12))

for m in range(T):
    for n in range(T):
        if m == n:
            Am[m][n] = complex(- (1 + delta(m))*(h_1/2))

for m in range(T):
    for n in range(T):
        if m == n:
            Fm[m][n] = complex((1 + delta(m))*(h_1/2))

for m in range(T):
    for n in range(T):
        Lm[m][n] = complex(quad(lambda y: cos((pi*n*(l_1*sin(alfa) + y*cos(alfa)))/(h_2)) * cos(m*pi*y/h_1) * exp(1j*ksi(n, k, h_2)*(l_1*cos(alfa) - y*sin(alfa) + l_2)),[0, h_1]))

for m in range(T):
    for n in range(T):
        Km[m][n] = complex(quad(lambda y: cos(- pi*n)*cos(m*pi*y/h_1)*exp(1j*ksi(n, k, l_1)*y), [0, h_1]))

for m in range(T):
    for n in range(T):
        Em[m][n] = complex(quad(lambda y: cos(pi*n*(-l_1*cos(alfa) + y*sin(alfa))/l_2) * cos(pi*m*y/h_1)*exp(1j*ksi(n, k, l_2)*(l_1*sin(alfa) + y*cos(alfa))), [0, h_1]))

systemone = np.hstack([Am, Bm, Cm, Dm, Fm, Lm, Km, Em, Gm, Hm, Rm, Xm])

'''
    Перейдемо до другого рівняння
'''

Am_1, Bm_1, Cm_1, Dm_1, Fm_1, Lm_1, Km_1, Em_1, Gm_1, Hm_1, Rm_1, Xm_1 = ([[0 for i in range(T)] for j in range(T)] for q in range(12))

for m in range(T):
    for n in range(T):
        if m == n:
            Am_1[m][n] = complex(- (1 + delta(m))*(h_1/2)*(-1j * ksi(m, k, h_1)))

for m in range(T):
    for n in range(T):
        if m == n:
            Fm_1[m][n] = complex((1 + delta(m))*(h_1/2)*(1j * ksi(m, k, h_1)))

'''for m in range(T):
    for n in range(T):
        Km_1[m][n] = complex(quad(lambda y: (pi*n/l_1)*(-sin(pi*n))*exp(1j*ksi(n, k, l_1)*y)*cos(pi*m*y/h_1), [0, h_1]))'''

for m in range(T):
    for n in range(T):
        Lm_1[m][n] = complex(quad(lambda y: exp(-1j*ksi(n, k, h_2)*(-l_1*cos(alfa)+y*sin(alfa)-l_2)) * ((pi*n*sin(alfa)/h_2)*sin(pi*n*(l_1*sin(alfa)+y*cos(alfa))/h_2) + cos(pi*n*(l_1*sin(alfa)+y*cos(alfa))/h_2)*(-1j*cos(alfa)*ksi(n, k, h_2))) * cos(m*pi*y/h_1), [0, h_1]))

for m in range(T):
    for n in range(T):
        Em_1[m][n] = complex(quad(lambda y: exp(1j*ksi(n, k, l_2)*(l_1*sin(alfa)+y*cos(alfa))) * ((pi*n*cos(alfa)/l_2)*(-sin(pi*n*(-l_1*cos(alfa)+y*sin(alfa))/l_2)) + cos(pi*n*(-l_1*cos(alfa)+y*sin(alfa))/l_2)*(-1j*sin(alfa)*ksi(n, k, l_2))) * cos(m*pi*y/h_1), [0, h_1]))

systemtwo = np.hstack([Am_1, Bm_1, Cm_1, Dm_1, Fm_1, Lm_1, Km_1, Em_1, Gm_1, Hm_1, Rm_1, Xm_1])

'''
    Перейдемо до третього рівняння
'''

Am_2, Bm_2, Cm_2, Dm_2, Fm_2, Lm_2, Km_2, Em_2, Gm_2, Hm_2, Rm_2, Xm_2 = ([[0 for i in range(T)] for j in range(T)] for q in range(12))

for m in range(T):
    for n in range(T):
        Fm_2[m][n] = complex(quad(lambda y: cos(pi*n*(l_2*sin(alfa) + y*cos(alfa))/h_1) * exp(1j*ksi(n,k,h_1)*(l_2*cos(alfa) - y*sin(alfa) + l_1)) * cos(m*pi*y/h_2), [0, h_2]))

for m in range(T):
    for n in range(T):
        if m == n:
            Lm_2[m][n] = complex((1+delta(m))*(h_2/2))

for m in range(T):
    for n in range(T):
        Km_2[m][n] = complex(quad(lambda y: cos(pi*n*(l_2*cos(alfa) - y*sin(alfa))/l_1) * exp(1j*ksi(n,k,l_1)*(l_2*sin(alfa) + y*cos(alfa))) * cos(pi*m*y/h_2), [0, h_2]))

for m in range(T):
    for n in range(T):
        Em_2[m][n] = complex(quad(lambda y: cos(pi*n) * cos(pi*m*y/h_2) * exp(1j*ksi(n,k,l_2)*y), [0, h_2]))

for m in range(T):
    for n in range(T):
        if m == n:
            Cm_2[m][n] = complex(- (1+delta(m))*(h_2/2))

for m in range(T):
    for n in range(T):
        if m == n:
            Dm_2[m][n] = complex(- (1+delta(m))*(h_2/2)*exp(1j*ksi(n, k, h_2)*l_3))

systemthree = np.hstack([Am_2, Bm_2, Cm_2, Dm_2, Fm_2, Lm_2, Km_2, Em_2, Gm_2, Hm_2, Rm_2, Xm_2])

'''
    Перейдемо до рівняння 4
'''

Am_3, Bm_3, Cm_3, Dm_3, Fm_3, Lm_3, Km_3, Em_3, Gm_3, Hm_3, Rm_3, Xm_3 = ([[0 for i in range(T)] for j in range(T)] for q in range(12))

for m in range(T):
    for n in range(T):
        Fm_3[m][n] = complex(quad(lambda y: exp(1j*ksi(n,k,h_1)*(l_2*cos(alfa) - y*sin(alfa) + l_1)) * (
            (pi*n*sin(alfa)/h_1)*(-sin(pi*n*(l_2*sin(alfa) + y*cos(alfa))/h_1)) + cos(pi*n*(l_2*sin(alfa) + y*cos(alfa))/h_1)*(1j*ksi(n,k,h_1)*cos(alfa))
        ) * cos(pi*m*y/h_2), [0, h_2]))

for m in range(T):
    for n in range(T):
        if m == n:
            Lm_3[m][n] = complex((1+delta(n))*(h_2/2)*(-1j*ksi(n,k,h_2)))

for m in range(T):
    for n in range(T):
        Km_3[m][n] = complex(quad(lambda y: exp(1j*ksi(n, k, l_1)*(l_2*sin(alfa) + y*cos(alfa))) * (
            (pi*n*cos(alfa)/l_1)*(-sin(pi*n*(l_2*cos(alfa) - y*sin(alfa))/l_1)) + cos(pi*n*(l_2*cos(alfa) - y*sin(alfa))/l_1)*(1j*ksi(n,k,l_1)*sin(alfa))
        ) * cos(pi*m*y/h_2), [0, h_2]))

for m in range(T):
    for n in range(T):
        if m == n:
            Cm_3[m][n] = complex(-(1+delta(m))*(h_2/2)*(1j*ksi(n,k,h_2)))

for m in range(T):
    for n in range(T):
        if m == n:
            Dm_3[m][n] = complex(-(1+delta(m))*(h_2/2)*(-1j*ksi(n,k,h_2))*exp(1j*ksi(n,k,h_2)*l_3))

systemfour = np.hstack([Am_3, Bm_3, Cm_3, Dm_3, Fm_3, Lm_3, Km_3, Em_3, Gm_3, Hm_3, Rm_3, Xm_3])

'''
    Перейдемо до рівняння 5
'''

Am_4, Bm_4, Cm_4, Dm_4, Fm_4, Lm_4, Km_4, Em_4, Gm_4, Hm_4, Rm_4, Xm_4 = ([[0 for i in range(T)] for j in range(T)] for q in range(12))

for m in range(T):
    for n in range(T):
        Gm_4[m][n] = complex(quad(lambda y: cos(pi*n*(-y*cos(beta))/h_3) * cos(pi*m*y/h_2) * exp(1j*ksi(n,k,h_3)*y*sin(beta)), [0, h_2]))

for m in range(T):
    for n in range(T):
        Rm_4[m][n] = complex(quad(lambda y: cos(pi*m*y/h_2)*exp(-1j*ksi(n,k,l_5)*(y-h_2)), [0, h_2]))

for m in range(T):
    for n in range(T):
        Xm_4[m][n] = complex(quad(lambda y: cos(pi*n*(-y*sin(beta))/l_4) * cos(pi*m*y/h_2) * exp(1j*ksi(n,k,l_4)*(y*cos(beta) + h_3)), [0,h_2]))

for m in range(T):
    for n in range(T):
        if m == n:
            Cm_4[m][n] = complex(-(1+delta(m))*(h_2/2)*exp(1j*ksi(n,k,h_2)*l_3))

for m in range(T):
    for n in range(T):
        if m == n:
            Dm_4[m][n] = complex(-(1+delta(m))*(h_2/2))

for m in range(T):
    for n in range(T):
        if m == n:
            Hm_4[m][n] = complex((1+delta(m))*(h_2/2))

systemfive = np.hstack([Am_4, Bm_4, Cm_4, Dm_4, Fm_4, Lm_4, Km_4, Em_4, Gm_4, Hm_4, Rm_4, Xm_4])

'''
    Перейдемо до рівняння 6
'''

Am_5, Bm_5, Cm_5, Dm_5, Fm_5, Lm_5, Km_5, Em_5, Gm_5, Hm_5, Rm_5, Xm_5 = ([[0 for i in range(T)] for j in range(T)] for q in range(12))

for m in range(T):
    for n in range(T):
        Gm_5[m][n] = complex(quad(lambda y: exp(1j*ksi(n,k,h_3)*y*sin(beta)) * (
            (pi*n*sin(beta)/h_3)*(-sin(pi*n*(-y*cos(beta))/h_3)) + cos(pi*n*(-y*cos(beta))/h_3)*(1j*ksi(n,k,h_3)*cos(beta))
        ) * (cos(pi*m*y/h_2)), [0, h_2]))

for m in range(T):
    for n in range(T):
        Xm_5[m][n] = complex(quad(lambda y: exp(1j*ksi(n,k,l_4)*(y*cos(beta) + h_3)) * (
            (pi*n*cos(beta)/l_4)*sin(pi*n*(-y*sin(beta))/l_4) + cos(pi*n*(-y*sin(beta))/l_4)*(-1j*ksi(n,k,l_4)*sin(beta))
        ) * (cos(pi*m*y/h_2)), [0, h_2]))

for m in range(T):
    for n in range(T):
        if m == n:
            Cm_5[m][n] = complex(-(1+delta(m))*(h_2/2)*(1j*ksi(n,k,h_2))*exp(1j*ksi(n,k,h_2)*l_3))

for m in range(T):
    for n in range(T):
        if m == n:
            Dm_5[m][n] = complex(-(1+delta(m))*(h_2/2))*(-1j*ksi(n,k,h_2))

for m in range(T):
    for n in range(T):
        if m == n:
            Hm_5[m][n] = complex((1+delta(m))*(h_2/2)*(1j*ksi(n,k,h_2)))

systemsix = np.hstack([Am_5, Bm_5, Cm_5, Dm_5, Fm_5, Lm_5, Km_5, Em_5, Gm_5, Hm_5, Rm_5, Xm_5])

'''
    Перейдемо до рівняння 7
'''

Am_6, Bm_6, Cm_6, Dm_6, Fm_6, Lm_6, Km_6, Em_6, Gm_6, Hm_6, Rm_6, Xm_6 = ([[0 for i in range(T)] for j in range(T)] for q in range(12))

for m in range(T):
    for n in range(T):
        if m == n:
            Gm_6[m][n] = complex((1+delta(m))*(h_3/2))

for m in range(T):
    for n in range(T):
        Xm_6[m][n] = complex(quad(lambda y: cos(pi*m*y/h_3)*exp(-1j*ksi(n,k,l_4)*(y-h_3)), [0, h_3]))

for m in range(T):
    for n in range(T):
        if m == n:
            Bm_6[m][n] = complex(-(1+delta(m))*(h_3/2))

for m in range(T):
    for n in range(T):
        Hm_6[m][n] = complex(quad(lambda y: cos(pi*n*(-y*cos(beta))/h_2) * cos(pi*m*y/h_3) * exp(1j*ksi(n,k,h_2)*y*sin(beta)), [0, h_3]))

for m in range(T):
    for n in range(T):
        Rm_6[m][n] = complex(quad(lambda y: cos(pi*n*y*sin(beta)/l_5) * cos(pi*m*y/h_3) * exp(1j*ksi(n,k,l_5)*(y*cos(beta) + h_2)), [0, h_3]))

systemseven = np.hstack([Am_6, Bm_6, Cm_6, Dm_6, Fm_6, Lm_6, Km_6, Em_6, Gm_6, Hm_6, Rm_6, Xm_6])

'''
    Перейдемо до рівняння 8
'''

Am_7, Bm_7, Cm_7, Dm_7, Fm_7, Lm_7, Km_7, Em_7, Gm_7, Hm_7, Rm_7, Xm_7 = ([[0 for i in range(T)] for j in range(T)] for q in range(12))

for m in range(T):
    for n in range(T):
        if m == n:
            Gm_7[m][n] = complex((1+delta(m))*(h_3/2)*(-1j*ksi(n,k,h_3)))

for m in range(T):
    for n in range(T):
        Hm_7[m][n] = complex(quad(lambda y: exp(1j*ksi(n,k,h_2)*y*sin(beta)) * (
            (pi*n*sin(beta)/h_2)*sin(pi*n*(-y*cos(beta))/h_2) + cos(pi*n*(-y*cos(beta))/h_2)*(-1j*ksi(n,k,h_2)*cos(beta))
        ) * cos(pi*m*y/h_3), [0, h_3]))

for m in range(T):
    for n in range(T):
        if m == n:
            Bm_7[m][n] = complex(-(1+delta(m))*(h_3/2)*(1j*ksi(n,k,h_3)))

for m in range(T):
    for n in range(T):
        Rm_7[m][n] = complex(quad(lambda y: exp(1j*ksi(n,k,l_5)*(y*cos(beta) + h_2)) * (
            (pi*n*cos(beta)/l_5)*sin(pi*n*y*sin(beta)/l_5) + cos(pi*n*y*sin(beta)/l_5)*(1j*ksi(n,k,l_5)*sin(beta))
        ) * cos(pi*m*y/h_3), [0, h_3]))

systemeight = np.hstack([Am_7, Bm_7, Cm_7, Dm_7, Fm_7, Lm_7, Km_7, Em_7, Gm_7, Hm_7, Rm_7, Xm_7])

'''
    Перейдемо до рівняння 9
'''

Am_8, Bm_8, Cm_8, Dm_8, Fm_8, Lm_8, Km_8, Em_8, Gm_8, Hm_8, Rm_8, Xm_8 = ([[0 for i in range(T)] for j in range(T)] for q in range(12))

for m in range(T):
    for n in range(T):
        Lm_8[m][n] = complex(quad(lambda x: exp(-1j*ksi(n,k,h_2)*(x*cos(alfa) - l_2)) * (
            (pi*n*cos(alfa)/h_2)*(-sin(pi*n*(-x*sin(alfa))/h_2)) + cos(pi*n*(-x*sin(alfa))/h_2)*(-1j*ksi(n,k,h_2)*sin(alfa))
        ) * cos(pi*m*x/l_1), [(-l_1), 0]))

for m in range(T):
    for n in range(T):
        Em_8[m][n] = complex(quad(lambda x: exp(1j*ksi(n,k,l_2)*(-x*sin(alfa))) * (
            (pi*n*sin(alfa)/l_2)*(-sin(pi*n*x*cos(alfa)/l_2)) + cos(pi*n*x*cos(alfa)/l_2)*(1j*ksi(n,k,l_2)*cos(alfa))
        ) * cos(pi*m*x/l_1), [(-l_1), 0]))

for m in range(T):
    for n in range(T):
        if m == n:
            Km_8[m][n] = complex((1+delta(m))*(l_1/2)*(1j*ksi(n,k,l_1)))

systemnine = np.hstack([Am_8, Bm_8, Cm_8, Dm_8, Fm_8, Lm_8, Km_8, Em_8, Gm_8, Hm_8, Rm_8, Xm_8])

'''
    Перейдемо до рівняння 10
'''

Am_9, Bm_9, Cm_9, Dm_9, Fm_9, Lm_9, Km_9, Em_9, Gm_9, Hm_9, Rm_9, Xm_9 = ([[0 for i in range(T)] for j in range(T)] for q in range(12))

for m in range(T):
    for n in range(T):
        Fm_9[m][n] = complex(quad(lambda x: exp(1j*ksi(n,k,h_1)*(x*cos(alfa) + l_1)) * (
            (pi*n*cos(alfa)/h_1)*(-sin(pi*n*x*sin(alfa)/h_1)) + cos(pi*n*x*sin(alfa)/h_1)*(-1j*ksi(n,k,h_1)*sin(alfa))
        ) * cos(pi*m*x/l_2), [0, l_2]))

for m in range(T):
    for n in range(T):
        Km_9[m][n] = complex(quad(lambda x: exp(1j*ksi(n,k,l_1)*x*sin(alfa)) * (
            (pi*n*sin(alfa)/l_1)*sin(pi*n*x*cos(alfa)/l_1) + cos(pi*n*x*cos(alfa)/l_1)*(1j*ksi(n,k,l_1)*cos(alfa))
        ) * cos(pi*m*x/l_2), [0, l_2]))

for m in range(T):
    for n in range(T):
        if m == n:
            Em_9[m][n] = complex((1+delta(m))*(l_2/2)*(1j*ksi(n,k,l_2)))

systemten = np.hstack([Am_9, Bm_9, Cm_9, Dm_9, Fm_9, Lm_9, Km_9, Em_9, Gm_9, Hm_9, Rm_9, Xm_9])

'''
    Перейдемо до рівняння 11
'''

Am_10, Bm_10, Cm_10, Dm_10, Fm_10, Lm_10, Km_10, Em_10, Gm_10, Hm_10, Rm_10, Xm_10 = ([[0 for i in range(T)] for j in range(T)] for q in range(12))

for m in range(T):
    for n in range(T):
        Gm_10[m][n] = complex(quad(lambda x: exp(1j*ksi(n,k,h_3)*((x - l_2 - l_3)*cos(beta) + h_2*sin(beta))) * (
            (pi*n*cos(beta)/h_3)*sin(pi*n*((x - l_2 - l_3)*sin(beta) - h_2*cos(beta))/h_3) + cos(pi*n*((x - l_2 - l_3)*sin(beta) - h_2*cos(beta))/h_3)*(1j*ksi(n,k,h_3)*sin(beta))
        ) * cos(pi*m*(x - l_2 - l_3)/l_5), [l_2+l_3, l_2+l_3+l_5]))

for m in range(T):
    for n in range(T):
        Xm_10[m][n] = complex(quad(lambda x: exp(-1j*ksi(n,k,l_4)*((x - l_2 - l_3)*sin(beta) - h_2*cos(beta) - h_3)) * (
            (pi*n*sin(beta)/l_4)*(sin(pi*n*(-(x - l_2 - l_3)*cos(beta) - h_2*sin(beta))/l_4)) + cos(pi*n*(-(x - l_2 - l_3)*cos(beta) - h_2*sin(beta))/l_4)*(1j*ksi(n,k,l_4)*cos(beta))
        ) * cos(pi*m*(x - l_2 - l_3)/l_5), [l_2+l_3, l_2+l_3+l_5]))

for m in range(T):
    for n in range(T):
        if m == n:
            Rm_10[m][n] = complex((1+delta(m))*(l_5/2)*(-1j*ksi(n,k,l_5)))

systemeleven = np.hstack([Am_10, Bm_10, Cm_10, Dm_10, Fm_10, Lm_10, Km_10, Em_10, Gm_10, Hm_10, Rm_10, Xm_10])

'''
    Перейдемо до рівняння 12
'''

Am_11, Bm_11, Cm_11, Dm_11, Fm_11, Lm_11, Km_11, Em_11, Gm_11, Hm_11, Rm_11, Xm_11 = ([[0 for i in range(T)] for j in range(T)] for q in range(12))

for m in range(T):
    for n in range(T):
        Hm_11[m][n] = complex(quad(lambda x: exp(1j*ksi(n,k,h_2)*(-x*cos(beta) + h_3*sin(beta))) * (
            (pi*n*cos(beta)/h_2)*sin(pi*n*(-x*sin(beta) - h_3*cos(beta))/h_2) + cos(pi*n*(-x*sin(beta) - h_3*cos(beta))/h_2)*(1j*ksi(n,k,h_2)*sin(beta))
        ) * cos(pi*m*x/l_4), [-l_4,0]))

for m in range(T):
    for n in range(T):
        Rm_11[m][n] = complex(quad(lambda x: exp(-1j*ksi(n,k,l_5)*(-x*sin(beta) - h_3*cos(beta) - h_2)) * (
            (pi*n*sin(beta)/l_5)*(-sin(pi*n*(-x*cos(beta) + h_3*sin(beta))/l_5)) + cos(pi*n*(-x*cos(beta) + h_3*sin(beta))/l_5)*(1j*ksi(n,k,l_5)*cos(beta))
        ) * cos(pi*m*x/l_4), [-l_4, 0]))

for m in range(T):
    for n in range(T):
        if m == n:
            Xm_11[m][n] = complex((1+delta(m))*(l_4/2)*(-1j*ksi(n,k,l_4)))

systemtwelve = np.hstack([Am_11, Bm_11, Cm_11, Dm_11, Fm_11, Lm_11, Km_11, Em_11, Gm_11, Hm_11, Rm_11, Xm_11])

'''
    Знайдемо загальну систему
'''

SYSTEM = np.vstack([systemone, systemtwo, systemthree, systemfour, systemfive, systemsix, systemseven, systemeight, systemnine, systemten, systemeleven, systemtwelve])

'''
    Побудуємо праву частину системи для подальшого розвʼязку
'''

G_1, G_2, G_3, G_4, G_5, G_6, G_7, G_8, G_9, G_10, G_11, G_12 = ([[0] for j in range(T)] for q in range(12))

for i in range(T):
    G_1[i][0] = complex(h_1*delta(i))

for i in range(T):
    G_2[i][0] = complex(1j*k*h_1*delta(i))

Y = np.vstack([G_1, G_2, G_3, G_4, G_5, G_6, G_7, G_8, G_9, G_10, G_11, G_12])

'''
    Знайдемо розвʼязки системи
'''

X = np.linalg.solve(SYSTEM, Y)

A_n, B_n, C_n, D_n, F_n, L_n, K_n, E_n, G_n, H_n, R_n, X_n = ([] for i  in range(12))

for i in range(0, T):
    A_n.append(X[i][0])

for i in range(T, 2*T):
    B_n.append(X[i][0])

for i in range(2*T, 3*T):
    C_n.append(X[i][0])

for i in range(3*T, 4*T):
    D_n.append(X[i][0])

for i in range(4*T, 5*T):
    F_n.append(X[i][0])

for i in range(5*T, 6*T):
    L_n.append(X[i][0])

for i in range(6*T, 7*T):
    K_n.append(X[i][0])

for i in range(7*T, 8*T):
    E_n.append(X[i][0])

for i in range(8*T, 9*T):
    G_n.append(X[i][0])

for i in range(9*T, 10*T):
    H_n.append(X[i][0])

for i in range(10*T, 11*T):
    R_n.append(X[i][0])

for i in range(11*T, 12*T):
    X_n.append(X[i][0])

'''
    Перейдемо до пошуку енергетичних характеристик
'''

W_n = 0

W = 0

for n in range(0,T):
    W_n = (eps(n)*h_3*(ksi(n, k, h_3).real)*((abs(B_n[n]))**2))/(eps(0)*h_1*np.real(ksi(0,k,h_1)))
    W += W_n

print('Коеф. проходження: ',W)

V_n = 0

V = 0

for n in range(0,T):
    V_n = (eps(n)*(ksi(n, k, h_1).real)*((abs(A_n[n]))**2))/(eps(0)*ksi(0,k,h_1))
    V += V_n

print('Коеф. відбиття: ',V)

print(W+V)

'''
    Перейдемо до пошуку похибок
        та побудови графіків
'''

'''
    Розпочнемо з x = - l_1
'''

def P_0(x):
    return np.exp(1j*k*(x + l_1))

def P_1(x,y):
    P_1_first = np.exp(1j*k*(x + l_1))
    for n in range(T):
        P_1_first += A_n[n]*np.cos(n*np.pi*y/h_1)*np.exp(-1j*ksi(n,k,h_1)*(x+l_1))
    return P_1_first

def P_2(x,y):
    P_2_first = 0
    for n in range(T):
        P_2_first += (F_n[n]*np.cos(n*np.pi*y/h_1)*np.exp(1j*ksi(n,k,h_1)*(x+l_1)) + K_n[n]*np.cos(np.pi*n*x/l_1)*np.exp(1j*ksi(n,k,l_1)*y)
        + L_n[n]*np.cos(n*np.pi*(-x*np.sin(alfa) + y*np.cos(alfa))/h_2)*np.exp(-1j*ksi(n,k,h_2)*(x*np.cos(alfa) + y*np.sin(alfa) - l_2))
        + E_n[n]*np.cos(n*np.pi*(x*np.cos(alfa) + y*np.sin(alfa))/l_2)*np.exp(1j*ksi(n,k,l_2)*(-x*np.sin(alfa) + y*np.cos(alfa))))
    return P_2_first

y_1 = np.linspace(0,h_1,100)
x_1 = abs(P_2(-l_1, y_1) - P_1(-l_1,y_1))/abs(P_0(-l_1))

plt.plot(y_1, x_1)
plt.xlabel('$y/\\lambda$',loc='right',fontsize=14,rotation=0)
plt.ylabel('$\\delta_{p}^{(1)}$',loc='top',fontsize=14,rotation=0)
plt.subplots_adjust(left=0.15)
plt.grid(True)
plt.show()

'''
    Знайдемо для х = l_2
'''

def P_0_1st(x,y):
    return np.exp(1j*k*(x*np.cos(alfa) - y*np.sin(alfa) + l_1))

def P_2_1st(x,y):
    P_2_1st_first = 0
    for n in range(T):
        P_2_1st_first += (F_n[n]*np.cos(np.pi*n*(x*np.sin(alfa) + y*np.cos(alfa))/h_1)*np.exp(1j*ksi(n,k,h_1)*(x*np.cos(alfa) - y*np.sin(alfa) + l_1))
        + K_n[n]*np.cos(np.pi*n*(x*np.cos(alfa) - y*np.sin(alfa))/l_1)*np.exp(1j*ksi(n,k,l_1)*(x*np.sin(alfa) + y*np.cos(alfa)))
        + L_n[n]*np.cos(np.pi*n*y/h_2)*np.exp(-1j*ksi(n,k,h_2)*(x - l_2)) + E_n[n]*np.cos(np.pi*n*x/l_2)*np.exp(1j*ksi(n,k,l_2)*y))
    return P_2_1st_first

def P_3(x,y):
    P_3_first = 0
    for n in range(T):
        P_3_first += (C_n[n]*np.cos(np.pi*n*y/h_2)*np.exp(1j*ksi(n,k,h_2)*(x - l_2))
        + D_n[n]*np.cos(np.pi*n*y/h_2)*np.exp(-1j*ksi(n,k,h_2)*(x - l_2 - l_3)))
    return P_3_first

y_2 = np.linspace(0, h_2, 100)
x_2 = abs(P_2_1st(l_2, y_2) - P_3(l_2, y_2))/abs(P_0_1st(l_2,y_2))

plt.plot(y_2, x_2)
plt.xlabel('$yʼ/\\lambda$',loc='right',fontsize=14,rotation=0)
plt.ylabel('$\\delta_{p}^{(2)}$',loc='top',fontsize=14,rotation=0)
plt.subplots_adjust(left=0.15)
plt.grid(True)
plt.show()

'''
    Знайдемо для l_2 + l_3
'''

def P_4_1st(x,y):
    P_4_1st_first = 0
    for n in range(T):
        P_4_1st_first += (G_n[n]*np.cos(np.pi*n*((x - l_2 - l_3)*np.sin(beta) - y*np.cos(beta))/h_3)*np.exp(1j*ksi(n,k,h_3)*((x - l_2 - l_3)*np.cos(beta) + y*np.sin(beta)))
        + X_n[n]*np.cos(np.pi*n*(-(x - l_2 - l_3)*np.cos(beta) - y*np.sin(beta))/l_4)*np.exp(-1j*ksi(n,k,l_4)*((x - l_2 - l_3)*np.sin(beta) - y*np.cos(beta) - h_3))
        + H_n[n]*np.cos(np.pi*n*y/h_2)*np.exp(1j*ksi(n,k,h_2)*(x - l_2 - l_3)) + R_n[n]*np.cos(np.pi*n*(x - l_2 - l_3)/l_5)*np.exp(-1j*ksi(n,k,l_5)*(y - h_2)))
    return P_4_1st_first

y_3 = np.linspace(0, h_2, 100)
x_3 = abs(P_4_1st(l_2+l_3, y_3) - P_3(l_2+l_3, y_3))/abs(P_0_1st(l_2+l_3, y_3))

plt.plot(y_3, x_3)
plt.xlabel('$yʼ/\\lambda$',loc='right',fontsize=14,rotation=0)
plt.ylabel('$\\delta_{p}^{(3)}$',loc='top',fontsize=14,rotation=0)
plt.subplots_adjust(left=0.15)
plt.grid(True)
plt.show()

'''
    Для 0
'''

def P_5(x,y):
    P_5_first = 0
    for n in range(T):
        P_5_first += B_n[n]*np.cos(np.pi*n*y/h_3)*np.exp(1j*ksi(n,k,h_3)*x)
    return P_5_first

def P_4_2st(x,y):
    P_4_2st_first = 0
    for n in range(T):
        P_4_2st_first += (G_n[n]*np.cos(np.pi*n*y/h_3)*np.exp(-1j*ksi(n,k,h_3)*x) + X_n[n]*np.cos(np.pi*n*x/l_4)*np.exp(-1j*ksi(n,k,l_4)*(y-h_3))
        + H_n[n]*np.cos(np.pi*n*(-x*np.sin(beta) - y*np.cos(beta))/h_2)*np.exp(1j*ksi(n,k,h_2)*(-x*np.cos(beta)+y*np.sin(beta)))
        + R_n[n]*np.cos(np.pi*n*(-x*np.cos(beta)+y*np.sin(beta))/l_5)*np.exp(-1j*ksi(n,k,l_5)*(-x*np.sin(beta) - y*np.cos(beta) - h_2)))
    return P_4_2st_first

def P_0_2st(x,y):
    return np.exp(1j*k*(np.cos(alfa)*(-x*np.cos(beta)+y*np.sin(beta)+l_2+l_3) - np.sin(alfa)*(-x*np.sin(beta)-y*np.cos(beta)) + l_1))

y_4 = np.linspace(0, h_3, 100)
x_4 = abs(P_5(0,y_4) - P_4_2st(0,y_4))/abs(P_0_2st(0, y_4))

plt.plot(y_4, x_4)
plt.xlabel('$yʼʼ/\\lambda$',loc='right',fontsize=14,rotation=0)
plt.ylabel('$\\delta_{p}^{(4)}$',loc='top',fontsize=14,rotation=0)
plt.subplots_adjust(left=0.15)
plt.grid(True)
plt.show()

'''
    Перейдемо до пошуку похибок за коливальною швидкістю
'''

'''
    Розпочнемо з x = -l_1
'''

def V_0(x):
    return 1j*k*np.exp(1j*k*(x + l_1))

def V_1(x,y):
    V_1_first = 1j*k*np.exp(1j*k*(x + l_1))
    for n in range(T):
        V_1_first += A_n[n]*np.cos(np.pi*n*y/h_1)*(-1j*ksi(n,k,h_1))*np.exp(-1j*ksi(n,k,h_1)*(x + l_1))
    return V_1_first

def V_2(x,y):
    V_2_first = 0
    for n in range(T):
        V_2_first += (F_n[n]*np.cos(np.pi*n*y/h_1)*(1j*ksi(n,k,h_1))*np.exp(1j*ksi(n,k,h_1)*(x + l_1)) + K_n[n]*(np.pi*n/l_1)*(-np.sin(np.pi*n*x/l_1))*np.exp(1j*ksi(n,k,l_1)*y)
        + L_n[n]*np.exp(-1j*ksi(n,k,h_2)*(x*np.cos(alfa) + y*np.sin(alfa) - l_2))*((np.sin(alfa)*np.pi*n/h_2)*np.sin(np.pi*n*(-x*np.sin(alfa) + y*np.cos(alfa))/h_2) + np.cos(np.pi*n*(-x*np.sin(alfa) + y*np.cos(alfa))/h_2)*(-1j*ksi(n,k,h_2)*np.cos(alfa)))
        + E_n[n]*np.exp(1j*ksi(n,k,l_2)*(-x*np.sin(alfa) + y*np.cos(alfa)))*((np.pi*n*np.cos(alfa)/l_2)*(-np.sin(np.pi*n*(x*np.cos(alfa) + y*np.sin(alfa))/l_2)) + np.cos(np.pi*n*(x*np.cos(alfa) + y*np.sin(alfa))/l_2)*(-1j*ksi(n,k,l_2)*np.sin(alfa))))
    return V_2_first

y_1v = np.linspace(0,h_1,100)
x_1v = abs(V_2(-l_1, y_1v) - V_1(-l_1,y_1v))/abs(V_0(-l_1))

plt.plot(y_1v, x_1v)
plt.xlabel('$y/\\lambda$',loc='right',fontsize=14,rotation=0)
plt.ylabel('$\\delta_{V}^{(1)}$',loc='top',fontsize=14,rotation=0)
plt.subplots_adjust(left=0.15)
plt.grid(True)
plt.show()

'''
    Знайдемо для x = l_2
'''

def V_0_1st(x,y):
    return (1j*k*np.cos(alfa))*np.exp(1j*k*(x*np.cos(alfa) - y*np.sin(alfa) + l_1))

def V_3(x,y):
    V_3_first = 0
    for n in range(T):
        V_3_first += (C_n[n]*np.cos(np.pi*n*y/h_2)*(1j*ksi(n,k,h_2))*np.exp(1j*ksi(n,k,h_2)*(x - l_2))
        + D_n[n]*np.cos(np.pi*n*y/h_2)*(-1j*ksi(n,k,h_2))*np.exp(-1j*ksi(n,k,h_2)*(x - l_2 - l_3)))
    return V_3_first

def V_2_1st(x,y):
    V_2_1st_first = 0
    for n in range(T):
        V_2_1st_first += (L_n[n]*np.cos(np.pi*n*y/h_2)*(-1j*ksi(n,k,h_2))*np.exp(-1j*ksi(n,k,h_2)*(x - l_2)) + E_n[n]*(np.pi*n/l_2)*(-np.sin(np.pi*n*x/l_2))*np.exp(1j*ksi(n,k,l_2)*y)
        + F_n[n]*np.exp(1j*ksi(n,k,h_1)*(x*np.cos(alfa) - y*np.sin(alfa) + l_1))*(
            (np.pi*n*np.sin(alfa)/h_1)*(-np.sin(np.pi*n*(x*np.sin(alfa) + y*np.cos(alfa))/h_1)) + np.cos(np.pi*n*(x*np.sin(alfa) + y*np.cos(alfa))/h_1)*(1j*ksi(n,k,h_1)*np.cos(alfa))
        ) + K_n[n]*np.exp(1j*ksi(n,k,l_1)*(x*np.sin(alfa) + y*np.cos(alfa)))*(
            (np.pi*n*np.cos(alfa)/l_1)*(-np.sin(np.pi*n*(x*np.cos(alfa) - y*np.sin(alfa))/l_1)) + np.cos(np.pi*n*(x*np.cos(alfa) - y*np.sin(alfa))/l_1)*(1j*ksi(n,k,l_1)*np.sin(alfa))
        ))
    return V_2_1st_first

y_2v = np.linspace(0, h_2, 100)
x_2v = abs(V_2_1st(l_2, y_2v) - V_3(l_2, y_2v))/abs(V_0_1st(l_2,y_2v))

plt.plot(y_2v, x_2v)
plt.xlabel('$yʼ/\\lambda$',loc='right',fontsize=14,rotation=0)
plt.ylabel('$\\delta_{V}^{(2)}$',loc='top',fontsize=14,rotation=0)
plt.subplots_adjust(left=0.15)
plt.grid(True)
plt.show()

'''
    Знайдемо для x = l_2 + l_3
'''

def V_4_1st(x,y):
    V_4_1st_first = 0
    for n in range(T):
        V_4_1st_first += (H_n[n]*np.cos(np.pi*n*y/h_2)*(1j*ksi(n,k,h_2))*np.exp(1j*ksi(n,k,h_2)*(x - l_2 - l_3)) + R_n[n]*(np.pi*n/l_5)*(-np.sin(np.pi*n*(x - l_2 - l_3)/l_5))*np.exp(-1j*ksi(n,k,l_5)*(y - h_2))
        + G_n[n]*np.exp(1j*ksi(n,k,h_3)*((x - l_2 - l_3)*np.cos(beta) + y*np.sin(beta)))*(
            (np.pi*n*np.sin(beta)/h_3)*(-np.sin(np.pi*n*((x - l_2 - l_3)*np.sin(beta) - y*np.cos(beta))/h_3)) + np.cos(np.pi*n*((x - l_2 - l_3)*np.sin(beta) - y*np.cos(beta))/h_3)*(1j*ksi(n,k,h_3)*np.cos(beta))
        ) + X_n[n]*np.exp(-1j*ksi(n,k,l_4)*((x - l_2 - l_3)*np.sin(beta) - y*np.cos(beta) - h_3))*(
            (np.pi*n*np.cos(beta)/l_4)*np.sin(np.pi*n*(-(x - l_2 - l_3)*np.cos(beta) - y*np.sin(beta))/l_4) + np.cos(np.pi*n*(-(x - l_2 - l_3)*np.cos(beta) - y*np.sin(beta))/l_4)*(-1j*ksi(n,k,l_4)*np.sin(beta))
        ))
    return V_4_1st_first

y_3v = np.linspace(0, h_2, 100)
x_3v = abs(V_4_1st(l_2+l_3, y_3v) - V_3(l_2+l_3, y_3v))/abs(V_0_1st(l_2+l_3, y_3v))

plt.plot(y_3v, x_3v)
plt.xlabel('$yʼ/\\lambda$',loc='right',fontsize=14,rotation=0)
plt.ylabel('$\\delta_{V}^{(3)}$',loc='top',fontsize=14,rotation=0)
plt.subplots_adjust(left=0.15)
plt.grid(True)
plt.show()

'''
    Знайдемо для x = 0
'''

def V_0_2st(x,y):
    return ((-1j*k*(np.cos(beta)*np.cos(alfa) - np.sin(beta)*np.sin(alfa)))*np.exp(1j*k*(
        (-x*np.cos(beta) + y*np.sin(beta) + l_2 + l_3)*np.cos(alfa) + (x*np.sin(beta) + y*np.cos(beta))*np.sin(alfa) + l_1
    )))

def V_5(x,y):
    V_5_first = 0
    for n in range(T):
        V_5_first += B_n[n]*np.cos(np.pi*n*y/h_3)*(1j*ksi(n,k,h_3))*np.exp(1j*ksi(n,k,h_3)*x)
    return V_5_first

def V_4_2st(x,y):
    V_4_2st_first= 0
    for n in range(T):
        V_4_2st_first += (G_n[n]*np.cos(np.pi*n*y/h_3)*(-1j*ksi(n,k,h_3))*np.exp(-1j*ksi(n,k,h_3)*x) + X_n[n]*(np.pi*n/l_4)*(-np.sin(np.pi*n*x/l_4))*np.exp(-1j*ksi(n,k,l_4)*(y - h_3))
        + H_n[n]*np.exp(1j*ksi(n,k,h_2)*(-x*np.cos(beta) + y*np.sin(beta)))*(
            (np.pi*n*np.sin(beta)/h_2)*(np.sin(np.pi*n*(-x*np.sin(beta) - y*np.cos(beta))/h_2)) + (np.cos(np.pi*n*(-x*np.sin(beta) - y*np.cos(beta))/h_2))*(-1j*ksi(n,k,h_2)*np.cos(beta))
        ) + R_n[n]*np.exp(1j*ksi(n,k,l_5)*(x*np.sin(beta) + y*np.cos(beta) + h_2))*(
            (np.pi*n*np.cos(beta)/l_5)*np.sin(np.pi*n*(-x*np.cos(beta) + y*np.sin(beta))/l_5) + np.cos(np.pi*n*(-x*np.cos(beta) + y*np.sin(beta))/l_5)*(1j*ksi(n,k,l_5)*np.sin(beta))
        ))
    return V_4_2st_first

y_4v = np.linspace(0, h_3, 100)
x_4v = abs(V_5(0,y_4v) - V_4_2st(0,y_4v))/abs(V_0_2st(0, y_4v))

plt.plot(y_4v, x_4v)
plt.xlabel('$yʼʼ/\\lambda$',loc='right',fontsize=14,rotation=0)
plt.ylabel('$\\delta_{V}^{(4)}$',loc='top',fontsize=14,rotation=0)
plt.subplots_adjust(left=0.15)
plt.grid(True)
plt.show()

'''
    Перейдемо до пошуку похибок на жорстких повернях
'''

'''
    Похибка по l_2
'''
def Vy_0_1st(x,y):
    return (-1j*k*np.sin(alfa))*np.exp(1j*k*(x*np.cos(alfa) - y*np.sin(alfa) + l_1))

def Vy_2_1st(x,y):
    V_2_1st_first = 0
    for n in range(T):
        V_2_1st_first += (L_n[n]*(np.pi*n/h_2)*(-np.sin(np.pi*n*y/h_2))*np.exp(-1j*ksi(n,k,h_2)*(x - l_2)) + E_n[n]*(np.cos(np.pi*n*x/l_2))*(1j*ksi(n,k,l_2))*np.exp(1j*ksi(n,k,l_2)*y)
        + F_n[n]*np.exp(1j*ksi(n,k,h_1)*(x*np.cos(alfa) - y*np.sin(alfa) + l_1))*(
            (np.pi*n*np.cos(alfa)/h_1)*(-np.sin(np.pi*n*(x*np.sin(alfa) + y*np.cos(alfa))/h_1)) + np.cos(np.pi*n*(x*np.sin(alfa) + y*np.cos(alfa))/h_1)*(-1j*ksi(n,k,h_1)*np.sin(alfa))
        ) + K_n[n]*np.exp(1j*ksi(n,k,l_1)*(x*np.sin(alfa) + y*np.cos(alfa)))*(
            (np.pi*n*np.sin(alfa)/l_1)*(np.sin(np.pi*n*(x*np.cos(alfa) - y*np.sin(alfa))/l_1)) + np.cos(np.pi*n*(x*np.cos(alfa) - y*np.sin(alfa))/l_1)*(1j*ksi(n,k,l_1)*np.cos(alfa))
        ))
    return V_2_1st_first

x_2vy = np.linspace(0, l_2, 100)
y_2vy = abs(Vy_2_1st(x_2vy, 0))/abs(V_0(-l_1))

plt.plot(x_2vy, y_2vy)
plt.xlabel('$xʼ/\\lambda$',loc='right',fontsize=14,rotation=0)
plt.ylabel('$\\delta_{V}^{(6)}$',loc='top',fontsize=14,rotation=0)
plt.subplots_adjust(left=0.15)
plt.grid(True)
plt.show()

'''
    Похибка по l_1
'''

def Vy_2(x,y):
    V_2_first = 0
    for n in range(T):
        V_2_first += (F_n[n]*(np.pi*n/h_1)*(-np.sin(np.pi*n*y/h_1))*np.exp(1j*ksi(n,k,h_1)*(x + l_1)) + K_n[n]*(np.cos(np.pi*n*x/l_1))*(1j*ksi(n,k,l_1))*np.exp(1j*ksi(n,k,l_1)*y)
        + L_n[n]*np.exp(-1j*ksi(n,k,h_2)*(x*np.cos(alfa) + y*np.sin(alfa) - l_2))*((np.cos(alfa)*np.pi*n/h_2)*(-np.sin(np.pi*n*(-x*np.sin(alfa) + y*np.cos(alfa))/h_2)) + np.cos(np.pi*n*(-x*np.sin(alfa) + y*np.cos(alfa))/h_2)*(-1j*ksi(n,k,h_2)*np.sin(alfa)))
        + E_n[n]*np.exp(1j*ksi(n,k,l_2)*(-x*np.sin(alfa) + y*np.cos(alfa)))*((np.pi*n*np.sin(alfa)/l_2)*(-np.sin(np.pi*n*(x*np.cos(alfa) + y*np.sin(alfa))/l_2)) + np.cos(np.pi*n*(x*np.cos(alfa) + y*np.sin(alfa))/l_2)*(1j*ksi(n,k,l_2)*np.cos(alfa))))
    return V_2_first

x_1vy = np.linspace(-l_1,0,100)
y_1vy = abs(Vy_2(x_1vy, 0))/abs(V_0(-l_1))

plt.plot(x_1vy, y_1vy)
plt.xlabel('$x/\\lambda$',loc='right',fontsize=14,rotation=0)
plt.ylabel('$\\delta_{V}^{(5)}$',loc='top',fontsize=14,rotation=0)
plt.subplots_adjust(left=0.15)
plt.grid(True)
plt.show()

'''
    Похибка по l_5
'''

def Vy_4_1st(x,y):
    V_4_1st_first = 0
    for n in range(T):
        V_4_1st_first += (H_n[n]*(np.pi*n/h_2)*(-np.sin(np.pi*n*y/h_2))*np.exp(1j*ksi(n,k,h_2)*(x - l_2 - l_3)) + R_n[n]*(np.cos(np.pi*n*(x - l_2 - l_3)/l_5))*(-1j*ksi(n,k,l_5))*np.exp(-1j*ksi(n,k,l_5)*(y - h_2))
        + G_n[n]*np.exp(1j*ksi(n,k,h_3)*((x - l_2 - l_3)*np.cos(beta) + y*np.sin(beta)))*(
            (np.pi*n*np.cos(beta)/h_3)*(np.sin(np.pi*n*((x - l_2 - l_3)*np.sin(beta) - y*np.cos(beta))/h_3)) + np.cos(np.pi*n*((x - l_2 - l_3)*np.sin(beta) - y*np.cos(beta))/h_3)*(1j*ksi(n,k,h_3)*np.sin(beta))
        ) + X_n[n]*np.exp(-1j*ksi(n,k,l_4)*((x - l_2 - l_3)*np.sin(beta) - y*np.cos(beta) - h_3))*(
            (np.pi*n*np.sin(beta)/l_4)*np.sin(np.pi*n*(-(x - l_2 - l_3)*np.cos(beta) - y*np.sin(beta))/l_4) + np.cos(np.pi*n*(-(x - l_2 - l_3)*np.cos(beta) - y*np.sin(beta))/l_4)*(1j*ksi(n,k,l_4)*np.cos(beta))
        ))
    return V_4_1st_first

x_3vy = np.linspace(l_2 + l_3, l_2 + l_3 + l_5, 100)
y_3vy = abs(Vy_4_1st(x_3vy, h_2))/abs(V_0(-l_1))

plt.plot(x_3vy, y_3vy)
plt.xlabel('$xʼ/\\lambda$',loc='right',fontsize=14,rotation=0)
plt.ylabel('$\\delta_{V}^{(7)}$',loc='top',fontsize=14,rotation=0)
plt.subplots_adjust(left=0.15)
plt.grid(True)
plt.show()

'''
    Похибка по l_4
'''

def Vy_4_2st(x,y):
    V_4_2st_first= 0
    for n in range(T):
        V_4_2st_first += (G_n[n]*(np.pi*n/h_3)*(-np.sin(np.pi*n*y/h_3))*np.exp(-1j*ksi(n,k,h_3)*x) + X_n[n]*(np.cos(np.pi*n*x/l_4))*(-1j*ksi(n,k,l_4))*np.exp(-1j*ksi(n,k,l_4)*(y - h_3))
        + H_n[n]*np.exp(1j*ksi(n,k,h_2)*(-x*np.cos(beta) + y*np.sin(beta)))*(
            (np.pi*n*np.cos(beta)/h_2)*(np.sin(np.pi*n*(-x*np.sin(beta) - y*np.cos(beta))/h_2)) + (np.cos(np.pi*n*(-x*np.sin(beta) - y*np.cos(beta))/h_2))*(1j*ksi(n,k,h_2)*np.sin(beta))
        ) + R_n[n]*np.exp(1j*ksi(n,k,l_5)*(x*np.sin(beta) + y*np.cos(beta) + h_2))*(
            (np.pi*n*np.sin(beta)/l_5)*(-np.sin(np.pi*n*(-x*np.cos(beta) + y*np.sin(beta))/l_5)) + np.cos(np.pi*n*(-x*np.cos(beta) + y*np.sin(beta))/l_5)*(1j*ksi(n,k,l_5)*np.cos(beta))
        ))
    return V_4_2st_first

x_4vy = np.linspace(-l_4, 0, 100)
y_4vy = abs(Vy_4_2st(x_4vy,h_3))/abs(V_0(-l_1))

plt.plot(x_4vy, y_4vy)
plt.xlabel('$xʼʼ/\\lambda$',loc='right',fontsize=14,rotation=0)
plt.ylabel('$\\delta_{V}^{(8)}$',loc='top',fontsize=14,rotation=0)
plt.subplots_adjust(left=0.15)
plt.grid(True)
plt.show()



'''
    Перейдемо  до побудови поля амплітуди тиску
'''

boundaries = [0]
x = 0
for i in range(20):
    x += 0.15
    boundaries.append(x)

cmap = plt.get_cmap('Blues', len(boundaries) - 1)
norm = matplotlib.colors.BoundaryNorm(boundaries, cmap.N)

xn = np.linspace(-1.5, -l_1, 500)
yn = np.linspace(0, h_1, 500)
X, Y = np.meshgrid(xn, yn)
fig, ax = plt.subplots(figsize = (5, 2))
ax.pcolormesh(X, Y, abs(P_1(X, Y)), cmap = cmap, norm = norm)

rotate_transform = trns.Affine2D().rotate_deg(45) + ax.transData
xn_4 = np.linspace(0, l_2, 500)
yn_4 = np.linspace(0, h_2, 500)
X_4, Y_4 = np.meshgrid(xn_4,yn_4)
ax.pcolormesh(X_4, Y_4, abs(P_2_1st(X_4, Y_4)), transform=rotate_transform, cmap = cmap, norm = norm)

xn_2 = np.linspace(-l_1, 0, 500)
yn_2 = np.linspace(0, h_1, 500)
X_1, Y_1 = np.meshgrid(xn_2, yn_2)
ax.pcolormesh(X_1, Y_1, abs(P_2(X_1, Y_1)), cmap = cmap, norm = norm)

xn_3_1 = np.linspace(l_2, l_2+l_3, 500)
yn_3_1 = np.linspace(0, h_2, 500)
X_2, Y_2 = np.meshgrid(xn_3_1, yn_3_1)
im = ax.pcolormesh(X_2, Y_2, abs(P_3(X_2, Y_2)), transform=rotate_transform, cmap = cmap, norm = norm)

xn_5 = np.linspace(l_2+l_3, l_2+l_3+l_5, 500)
yn_5 = np.linspace(0, h_2, 500)
X_5, Y_5 = np.meshgrid(xn_5, yn_5)
ax.pcolormesh(X_5, Y_5, abs(P_4_1st(X_5, Y_5)), transform=rotate_transform, cmap = cmap, norm = norm)

translate_to_third = trns.Affine2D().translate(l_2+l_3,0) + rotate_transform
second_rotate_transform = trns.Affine2D().rotate_deg(-45) + translate_to_third

'''xn_6 = np.linspace(0, 1, 100)
yn_6 = np.linspace(0, h_3, 100)
X_6, Y_6 = np.meshgrid(xn_6, yn_6)
ax.pcolormesh(X_6, Y_6, abs(P_5(X_6, Y_6)), transform=second_rotate_transform, cmap='Blues', vmin=0.5,vmax=1.5)'''

xn_7 = np.linspace(-l_4, 0, 500)
yn_7 = np.linspace(0, h_3, 500)
X_7, Y_7 = np.meshgrid(xn_7, yn_7)
ax.pcolormesh(X_7, Y_7, abs(P_4_2st(X_7, Y_7)), transform=second_rotate_transform, cmap = cmap, norm = norm)

xn_3_1 = np.linspace(l_2, l_2+l_3, 500)
yn_3_1 = np.linspace(0, h_2, 500)
X_2, Y_2 = np.meshgrid(xn_3_1, yn_3_1)
im = ax.pcolormesh(X_2, Y_2, abs(P_3(X_2, Y_2)), transform=rotate_transform, cmap = cmap, norm = norm)

xn = np.linspace(-1, -l_1, 500)
yn = np.linspace(0, h_1, 500)
X, Y = np.meshgrid(xn, yn)
ax.pcolormesh(X, Y, abs(P_1(X, Y)), cmap = cmap, norm = norm)

xn_6 = np.linspace(0, 1.1, 500)
yn_6 = np.linspace(0, h_3, 500)
X_6, Y_6 = np.meshgrid(xn_6, yn_6)
ax.pcolormesh(X_6, Y_6, abs(P_5(X_6, Y_6)), transform=second_rotate_transform, cmap = cmap, norm = norm)

cbar = fig.colorbar(im)
plt.xlabel('$x/\\lambda$',loc='right',fontsize=14,rotation=0)
plt.ylabel('$y/\\lambda$',loc='top',fontsize=14,rotation=0)
plt.show()


'''
    Перейдемо до пошуку поля вектора інтенсивності
'''
def Vy_1(x,y):
    Vy_1_first = 0
    for n in range(T):
        Vy_1_first += A_n[n]*(np.pi*n/h_1)*(-np.sin(np.pi*n*y/h_1))*np.exp(-1j*ksi(n,k,h_1)*(x + l_1))
    return Vy_1_first

def Vy_3(x,y):
    Vy_3_first = 0
    for n in range(T):
        Vy_3_first += (C_n[n]*(np.pi*n/h_2)*(-np.sin(np.pi*n*y/h_2))*np.exp(1j*ksi(n,k,h_2)*(x - l_2))
        + D_n[n]*(np.pi*n/h_2)*(-np.sin(np.pi*n*y/h_2))*np.exp(-1j*ksi(n,k,h_2)*(x - l_2 - l_3)))
    return Vy_3_first

def Vy_5(x,y):
    Vy_5_first = 0
    for n in range(T):
        Vy_5_first += B_n[n]*(np.pi*n*y/h_3)*(-np.sin(np.pi*n*y/h_3))*np.exp(1j*ksi(n,k,h_3)*x)
    return Vy_5_first

fig2, ax2 = plt.subplots()
rotate_transform = trns.Affine2D().rotate_deg(45) + ax2.transData
translate_to_third = trns.Affine2D().translate(l_2+l_3,0) + rotate_transform
second_rotate_transform = trns.Affine2D().rotate_deg(-45) + translate_to_third

xt_2 = np.linspace(-l_1,0,6)
yt_2 = np.linspace(0,h_1,10)
X_2, Y_2 = np.meshgrid(xt_2, yt_2)
U_1 = ((P_2(X_2, Y_2) * ((1/1j*V_2(X_2, Y_2)).conjugate())).real)/abs((P_0(X_2) * (1/1j*V_0(X_2)).conjugate()).real)
W_1 = ((P_2(X_2, Y_2) * ((1/1j*Vy_2(X_2, Y_2)).conjugate())).real)/abs((P_0(X_2) * (1/1j*V_0(X_2)).conjugate()).real)
U_1 = np.ma.masked_where(Y_2 > abs((h_1/l_1)*X_2), U_1)
W_1 = np.ma.masked_where(Y_2 > abs((h_1/l_1)*X_2), W_1)
ax2.quiver(X_2, Y_2, U_1, W_1, scale = 40, width = 0.003) #масштаб 50

xt_3 = np.linspace(0,l_2,6)
yt_3 = np.linspace(0,h_2,10)
X_3, Y_3 = np.meshgrid(xt_3, yt_3)
U_2_1 = ((P_2_1st(X_3, Y_3) * ((1/1j*V_2_1st(X_3, Y_3)).conjugate())).real)/abs((P_0_1st(X_3,Y_3) * (1/1j*V_0_1st(X_3,Y_3)).conjugate()).real)
W_2_1 = ((P_2_1st(X_3, Y_3) * ((1/1j*Vy_2_1st(X_3, Y_3)).conjugate())).real)/abs((P_0_1st(X_3,Y_3) * (1/1j*V_0_1st(X_3,Y_3)).conjugate()).real)
U_2_1 = np.ma.masked_where(X_3 < abs((l_2/h_2)*Y_3), U_2_1)
W_2_1 = np.ma.masked_where(X_3 < abs((l_2/h_2)*Y_3), W_2_1)
U_2 = np.zeros_like(U_2_1)
W_2 = np.zeros_like(W_2_1)
for i in range(len(U_2_1)):
    for j in range(len(U_2_1[1])):
        U_2[i,j] = U_2_1[i,j]*np.cos(alfa) - W_2_1[i,j]*np.sin(alfa)
        W_2[i,j] = U_2_1[i,j]*np.sin(alfa) + W_2_1[i,j]*np.cos(alfa)
ax2.quiver(X_3, Y_3, U_2, W_2, transform=rotate_transform, scale = 40, width = 0.003) #масштаб 70

xt_1 = np.linspace(-1, -l_1, 10)
yt_1 = np.linspace(0, h_1, 10)
X_4, Y_4 = np.meshgrid(xt_1, yt_1)
U_3 = ((P_1(X_4, Y_4) * ((1/1j*V_1(X_4, Y_4)).conjugate())).real)/abs((P_0(X_4) * (1/1j*V_0(X_4)).conjugate()).real)
W_3 = ((P_1(X_4, Y_4) * ((1/1j*Vy_1(X_4, Y_4)).conjugate())).real)/abs((P_0(X_4) * (1/1j*V_0(X_4)).conjugate()).real)
ax2.quiver(X_4, Y_4, U_3, W_3, scale = 40, width = 0.003) #масштаб 40

xt_4_1 = np.linspace(l_2, l_2+l_3, 10)
yt_4_1 = np.linspace(0, h_2, 10)
X_5, Y_5 = np.meshgrid(xt_4_1, yt_4_1)
U_4_1 = ((P_3(X_5, Y_5) * ((1/1j*V_3(X_5, Y_5)).conjugate())).real)/abs((P_0_1st(X_5,Y_5) * (1/1j*V_0_1st(X_5,Y_5)).conjugate()).real)
W_4_1 = ((P_3(X_5, Y_5) * ((1/1j*Vy_3(X_5, Y_5)).conjugate())).real)/abs((P_0_1st(X_5,Y_5) * (1/1j*V_0_1st(X_5,Y_5)).conjugate()).real)
U_4 = np.zeros_like(U_4_1)
W_4 = np.zeros_like(W_4_1)
for i in range(len(U_4_1)):
    for j in range(len(U_4_1[1])):
        U_4[i,j] = U_4_1[i,j]*np.cos(alfa) - W_4_1[i,j]*np.sin(alfa)
        W_4[i,j] = U_4_1[i,j]*np.sin(alfa) + W_4_1[i,j]*np.cos(alfa)
ax2.quiver(X_5, Y_5, U_4, W_4,transform=rotate_transform,  scale = 40, width = 0.003) #масштаб 70, transform=rotate_transform,


xt_5 = np.linspace(l_2+l_3, l_2+l_3+l_5, 6)
yt_5 = np.linspace(0, h_2, 10)
X_6, Y_6 = np.meshgrid(xt_5, yt_5)
U_5_1 = ((P_4_1st(X_6, Y_6) * ((1/1j*V_4_1st(X_6, Y_6)).conjugate())).real)/abs((P_0_1st(X_6,Y_6) * (1/1j*V_0_1st(X_6,Y_6)).conjugate()).real)
W_5_1 = ((P_4_1st(X_6, Y_6) * ((1/1j*Vy_4_1st(X_6, Y_6)).conjugate())).real)/abs((P_0_1st(X_6,Y_6) * (1/1j*V_0_1st(X_6,Y_6)).conjugate()).real)
U_5_1 = np.ma.masked_where((X_6 - l_2 - l_3) > (l_5/h_2)*Y_6, U_5_1)
W_5_1 = np.ma.masked_where((X_6 - l_2 - l_3) > (l_5/h_2)*Y_6, W_5_1)
U_5 = np.zeros_like(U_5_1)
W_5 = np.zeros_like(W_5_1)
for i in range(len(U_5_1)):
    for j in range(len(U_5_1[1])):
        U_5[i,j] = U_5_1[i,j]*np.cos(alfa) - W_5_1[i,j]*np.sin(alfa)
        W_5[i,j] = U_5_1[i,j]*np.sin(alfa) + W_5_1[i,j]*np.cos(alfa)
ax2.quiver(X_6, Y_6, U_5, W_5, transform=rotate_transform, scale = 40, width = 0.003) #масштаб 80

xt_6 = np.linspace(-l_4, 0, 6)
yt_6 = np.linspace(0, h_3, 10)
X_7, Y_7 = np.meshgrid(xt_6, yt_6)
U_6 = ((P_4_2st(X_7, Y_7) * ((1/1j*V_4_2st(X_7, Y_7)).conjugate())).real)/abs((P_0_2st(X_7,Y_7) * ((1/1j*V_0_2st(X_7,Y_7)).conjugate())).real)
W_6 = ((P_4_2st(X_7, Y_7) * ((1/1j*Vy_4_2st(X_7, Y_7)).conjugate())).real)/abs((P_0_2st(X_7,Y_7) * ((1/1j*V_0_2st(X_7,Y_7)).conjugate())).real)
U_6 = np.ma.masked_where(Y_7 < abs((h_3/l_4)*X_7), U_6)
W_6 = np.ma.masked_where(Y_7 < abs((h_3/l_4)*X_7), W_6)
ax2.quiver(X_7, Y_7, U_6, W_6, transform=second_rotate_transform, scale = 40, width = 0.003) #масштаб 200

xt_7 = np.linspace(0, 0.5, 10)
yt_7 = np.linspace(0, h_3, 10)
X_8, Y_8 = np.meshgrid(xt_7, yt_7)
U_7 = ((P_5(X_8, Y_8) * ((1/1j*V_5(X_8, Y_8)).conjugate())).real)/abs((P_0_2st(X_8,Y_8) * ((1/1j*V_0_2st(X_8,Y_8)).conjugate())).real)
W_7 = ((P_5(X_8, Y_8) * ((1/1j*Vy_5(X_8, Y_8)).conjugate())).real)/abs((P_0_2st(X_8,Y_8) * ((1/1j*V_0_2st(X_8,Y_8)).conjugate())).real)
ax2.quiver(X_8, Y_8, U_7, W_7, transform=second_rotate_transform, scale = 40, width = 0.003) #масштаб 200

plt.xlabel('$x/\\lambda$',loc='right',fontsize=14,rotation=0)
plt.ylabel('$y/\\lambda$',loc='top',fontsize=14,rotation=0)
plt.show()

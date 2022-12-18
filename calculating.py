#Designing a Laval nozzle for a given value of thrust in the design mode

import math#standart libriry
import eel #pip install eel

from array import*
import numpy as np
import matplotlib.pyplot as plt
#eel.init("web")

def DigitsFixed(num, digits):
    return f"{num:.{digits}f}"


def grad(val):
    return val*180/np.pi


def rad(val):
    return val*np.pi/180


def mh_forM(a, b, ε, k):
    fa = nu_forroot(a,k)
    fb = nu_forroot(b,k)
    fx = 0.0
    x = 0.0
    if fa*fb>0:

        print("Error mh \n", fa, fb)
        return 0.0
    else:
        for i in range(1, 100):
            x = a-fa*((b-a)/(fb-fa))
            fx = nu_forroot(x,k)
            if fx*fa > 0:
                a = x
                fa = fx
            else:
                b = x
                fb = fx
            if b - a < ε:
                return (a+b)/2
            if abs(fx) < ε:
                return x


def mh_forl(r, rcr,a, b, ε, k):
    q1 = np.power(rcr/r, 2)
    fa = q(a,k, q1)
    fb = q(b,k, q1)
    fx = 0.0
    x = 0.0
    if fa*fb>0:

        print("Error mh \n", fa, fb)
        return 0.0
    else:
        for i in range(1, 500):
            x = a-fa*((b-a)/(fb-fa))
            fx = q(x,k, q1)
            if fx*fa > 0:
                a = x
                fa = fx
            else:
                b = x
                fb = fx
            if b - a < ε:
                return (a+b)/2
            if abs(fx) < ε:
                return x

def q(f, k, q1):
    return f*np.power((k+1)/2, 1/(k-1))*np.power(1-np.power(f, 2)*(k-1)/(k+1), 1/(k-1)) - q1


def m(k, R):
    return pow(2/(k+1), (k+1)/2*(k-1))*pow(k/R, 1/2)  #Расходный комплекс


def nu(M, k):
    return (np.power((k+1)/(k-1), 1/2)*np.arctan(np.power((k-1)/(k+1)*(M*M-1), 1/2))-np.arctan(np.power(M*M-1, 1/2)))


def nu_forroot(M,k):
    return grad(1/3*(np.power((k+1)/(k-1), 1/2)*np.arctan(np.power((k-1)/(k+1)*(M*M-1), 1/2))-np.arctan(np.power(M*M-1, 1/2)))) - 35

def solve(M, k, R, Tk, P, hi, muc, pk, R2, aa):
    p = []
    p.append(pow((k+1)*(M**2)/(2+(M**2)*(k-1)), 1/2)) #λ определение 0
    p.append(p[0]*pow((k+1)/2, 1/(k-1))*pow(1-pow(p[0], 2)*(k-1)/(k+1), 1/(k-1))) #q(λ) определение 1
    p.append(1/pow(p[1], 1/2)) #Ra определение 2
    p.append(pow(2*k/(k+1)*R*Tk, 1/2)) #acrit 3
    p.append(p[0]*p[3]) #wa 4
    p.append(P/p[4]) #G 5
    p.append(p[5]*pow(hi*Tk, 1/2)/muc/pk/m(k,R)) #Fcr 6
    p.append(pow(p[6]/math.pi, 1/2)) #Rcr 7
    p.append(1/3*nu(M,k)) #am 8
    p.append(1 + R2*(1-math.cos(p[8]))) #rm 9
    p.append(2*math.tan(p[8])*math.tan(aa)*(p[2]-p[9])/(math.tan(p[8])-math.tan(aa))) #p 10
    p.append((p[9]*math.tan(p[8])-p[2]*math.tan(aa))/(math.tan(p[8])-math.tan(aa))) #c1 11
    p.append((p[9]-p[11])**2/p[10]) #c2 12
    p.append((p[2]-p[11])**2/p[10]) #xa 13
    p.append(R2*math.sin(p[8])) #ksi1 14
    p.append(abs(p[13]) + abs(p[14])) #L 15
    p.append(M) #16
    return p

def shortsolve(s,Ra):
    s.append(np.arctan(s[10]/2/(Ra-s[11]))) #aai 17
    s.append((Ra-1)/((s[2]-1.0))) #mi 18
    s[13] = np.power(Ra-s[11], 2)/s[10] - s[12]
    s[15] = s[13] + s[14]
    return s


def solvewastes(w, k, R, mu, Tk, pk, Tw, Ra):
    w.append(np.power((2*k/(k-1))*R*Tk, 1/2)) #wmaxi 19
    print('wmax = ', w[19])
    ro = pk/R/Tk
    print("ро = ", ro)
    w.append(w[19]*w[15]*ro*w[7]/mu) #Rewi 20
    print('Re = ', w[20]/np.power(10, 8))
    print("ksim0 = ", 0.008*(2.62*k*k/np.power(Tw, 1/3) - 1)*np.power(Ra-1, 1/2))
    ksim = 0.008*(2.62*k*k/np.power(Tw, 1/3) - 1)*np.power(Ra-1, 1/2)*np.power(w[18], 1/10)*(0.3+0.035*np.exp(3*np.power(w[18],3)))*(np.power(np.power(10, 8)/w[20], 2/10)) #ksimp
    if  0.000 < w[17] <= 0.3:
        ksip = 0.05*w[17]*(0.1+w[17])
    elif 0.3 < w[17] <= 0.7:
        ksip = 0.00375-0.15*w[17] + 0.075*w[17]*w[17]
    else:
        print("С углом чёт", w[17])
        ksip = 1

    w.append(ksim+ksip) #ksisum 21
    print(ksim, "+", ksip, " = ", w[21])

    return w
    

P = 10.0*pow(10,4) #тяга
pk = 40.0*pow(10,5) #давление в камере сгорания
Tk = 3100.0 #температура в камере сгорания
hi = 0.98 #коэффициент хи
pa = 0.8*pow(10,5) #давление на срезе сопла
R = 310.0 #Газовая постоянная
k = 1.20 #Коэффициент аддиабаты
mu = 6.0/np.power(10, 5) #коэффициент мю
muc = 1.0
Tw = 0.6 #
R1 = 0.6 # r1/rk
aa = 4.0*math.pi/180 #
tethavh = math.pi/4
R2 = 0.6
pil = pa/pk #
l = pow((1.0 - pow(pil, (k-1)/k))*((k+1)/(k-1)), 1/2)
M = pow(2/(k+1)*(l**2)/(1-(l**2)*(k-1)/(k+1)), 1/2)
param = solve(M, k, R, Tk, P, hi, muc, pk, R2, aa)
Mmax = mh_forM(M, 10*M, 0.1, k)
n = 6
family = [0] * n

for i in range(n):
    if i == 0:
        family[i] = (param)
    else:
        M1 = M+(Mmax-M)*(i+1)/n
        family[i] = solve(M1, k, R, Tk, P, hi, muc, pk, R2, aa)

for row in family:
    print(' '.join([str(DigitsFixed(elem,3)) for elem in row]))
print("\n")


t = np.arange(np.pi-param[8], 3/2*np.pi-tethavh, 0.01)
x = np.arange(0, param[13], 0.01)
xx = np.arange(-param[14]-R2,param[13], 0.01)
xxx = np.arange(param[13], param[13]+2, 0.01)
m = np.arange(M, 1.25*Mmax, 0.01)
sp = plt.subplot(221)
plt.plot(x, param[11] + pow(param[10]*(x+param[12]), 1/2), -param[14]+R2*np.sin(t), param[9]+R2*np.cos(param[8])+R2*np.cos(t), xx, 0*xx, '--', color = "black")
plt.plot(xxx, param[2] - param[13]*np.tan(aa) + (xxx+param[12])*np.tan(aa), '--', color = 'black')
plt.axis('equal')


sp = plt.subplot(222)
plt.title('αm = f(M)', fontsize=20, fontname='Times New Roman')
plt.xlabel('M', color='gray')
plt.ylabel('αm, град.',color='gray')
plt.plot(m, grad(1/3*nu(m, k)), '--', color = 'black')
plt.plot([M, Mmax], [grad(param[8]), 35], 'ro')
plt.text(Mmax, 35, f'Mmax = {DigitsFixed(Mmax, 3)}')
plt.text(M, grad(param[8]), f'M = {DigitsFixed(M, 3)}')

sp = plt.subplot(223)
xx = np.arange(-family[n-1][14]-R2,family[n-1][13], 0.01)
plt.plot(xx, 0*xx, '--', color="black")
for i in range(n):
    if i == 0:
        plt.plot(x, family[i][11] + pow(family[i][10]*(x+family[i][12]), 1/2), '--', -family[i][14]+R2*np.sin(t), family[i][9]+R2*np.cos(family[i][8])+R2*np.cos(t), '--', linewidth= 1,  color ="red")
    else:
        t = np.arange(np.pi-family[i][8], 3/2*np.pi-tethavh, 0.01)
        x = np.arange(0, family[i][13], 0.01)
        plt.plot(x, family[i][11] + pow(family[i][10]*(x+family[i][12]), 1/2), '--', -family[i][14]+R2*np.sin(t), family[i][9]+R2*np.cos(family[i][8])+R2*np.cos(t), '--', linewidth= 1,  color ="black")
#plt.axis('equal')
for i in range(n):
    if i == 0:
        family[i] = param
        family[i].append(aa)
        family[i].append(1.00)
    else:
        family[i] = shortsolve(family[i], param[2])
for row in family:
    print(' '.join([str(DigitsFixed(elem,2)) for elem in row]))
print("\n")
for i in range(n):
    if i == 0:
        x = np.arange(0, family[i][13], 0.01)
        plt.plot(x, family[i][11] + pow(family[i][10]*(x+family[i][12]), 1/2), -family[i][14]+R2*np.sin(t), family[i][9]+R2*np.cos(family[i][8])+R2*np.cos(t), color ="red")
    else:
        t = np.arange(np.pi-family[i][8], 3/2*np.pi-tethavh, 0.01)
        x = np.arange(0, family[i][13], 0.01)
        plt.plot(x, family[i][11] + pow(family[i][10]*(x+family[i][12]), 1/2), -family[i][14]+R2*np.sin(t), family[i][9]+R2*np.cos(family[i][8])+R2*np.cos(t), color ="red")
wastes = [0]*n
for i in range(n):
    family[i] = solvewastes(family[i], k, R, mu, Tk, pk, Tw, param[2])
    wastes[i] = family[i][21]

param = family[wastes.index(min(wastes))]

disc = 100
Max = [0.0]*disc
tl = [0.0]*disc
pl = [0.0]*disc
el = [0.0]*disc
la = [0.0]*disc
for i in range(disc):
    r = (param[11] + np.power(param[10]*(param[13]*i/disc+param[12]), 1/2))*param[7]
    la[i] = mh_forl(r,param[7],0.99, 2.43, 0.5, k)
    f = la[i]
    tl[i] = 1.0-pow(f,2)*(k-1.0)/(k+1.0)
    pl[i] = pow(1-pow(f, 2)*(k-1)/(k+1), k/(k-1))
    el[i] = pow(1-(k-1)/(k+1)*pow(f, 2), 1/(k-1))
    Max[i] = pow(2/(k+1)*(f**2)/(1-(f**2)*(k-1)/(k+1)), 1/2)
k1 = param[13]*param[7]/(la[disc-1] - la[0])
print( param[13],"*", param[7],"//",la[disc-1], "=", k1)
for i in range(disc):
    la[i] = k1*la[i]
    tl[i] = k1*tl[i]
    pl[i] = k1*pl[i]
    el[i] = k1*el[i]
    Max[i] = Max[i]*k1/3
sp = plt.subplot(224)
t = np.arange(np.pi-param[8], 3/2*np.pi-tethavh, 0.01)
x = np.arange(0, param[13], 0.01)
x1 = np.arange(0, param[13], param[13]/100)
xx = np.arange(-param[14]*param[7]-R2*param[7],param[13]*param[7], 0.01)
plt.plot(x*param[7], (param[11] + pow(param[10]*(x+param[12]), 1/2))*param[7], (-param[14]+R2*np.sin(t))*param[7], (param[9]+R2*np.cos(param[8])+R2*np.cos(t))*param[7], xx, 0*xx, '--', color = "black")
plt.plot(la-la[0], tl, '--', la-la[0], el, '--', la-la[0], pl, '--', la-la[0], Max, "--", linewidth = 1)
plt.text(la[6]-la[0],pl[20],  "π(λ)")
plt.text(la[2]-la[0], el[2], "ε(λ)")
plt.text( 0, tl[0], "τ(λ)")
plt.text(la[disc-1]-la[0], Max[disc-1], "M(λ)")
#plt.plot(x1*param[7] - param[14]*param[7], pl, '--', x1*param[7] - param[14]*param[7], el, '--', linewidth= 1)
plt.axis('equal')
plt.show()

#!/usr/bin/env python
# coding: utf-8

# In[133]:


# ME 320 Lab 
# Independent Project Program: Lumped Capacitance
# Created by Amanda Yaklin 
# December 6, 2024

import numpy as np
import matplotlib.pyplot as plt

# test airspeeds
v_air=np.array([1,3,5]) # m/s

# designated room temperature
Tinf = 22 #degrees C

# properties of air (300K)
k_air=26.3e-3 #W/mK
rho_air=1.1614 #kg/m3
cp_air=1.007 #kJ/kgK
mu_air=184.63e-7 #Ns/m2
alpha_air=22.5e-6 #m2/s
Pr_air=0.707 #unitless

# properties of brass (300K)
k_brass=110 #W/mK
rho_brass=8530 #kg/m3
cp_brass=380 #J/kgK
alpha_brass=33.9e-6 #m2/s

# properties of aluminum 2024-T6
k_al=177 #W/mK
rho_al=2770 #kg/m3
cp_al=875 #J/kgK
alpha_al=73e-7 #m2/s

# cylinder
cyvol=14476.68671e-9 #m3
cyD=25e-3 #m, diameter
cyL=29e-3 #m, height
cyAs=np.pi*cyD*cyL+2*np.pi*cyD**2/4


# cube
cbvol=15422e-09 #m3
cbL=25e-3 #m, side
cbAs=6*cbL**2

# hexagonal cylinder (aluminum)
hcvol=14359.84675e-09 #m3
hcH=47e-3 #m, height 
hcP_P=22e-3 #m point-to-point max of hexagon
hcAs=3406e-6 #m2

# hexagonal cylinder (brass)
hcvol2=1947e-09
hcP_P2=11.5e-3 
hcAs2=952.631e-6

# lumped capacitance equation
def theta(rho,volume,cp,h,As,t,Ti):
    return np.exp(-h*As*t/(rho*volume*cp))*(Ti-Tinf)+Tinf

# Reynolds number
def Re(l_char, velocity=v_air, mu=mu_air, rho=rho_air): # mu and rho are for the fluid
    out=rho*velocity*l_char/mu
    return out

# Find h
# Nu_D = CRe^m * Pr^(1/3)
Nu_D_cyl=0.683*Pr_air*Re(cyD)**0.466
Nu_D_cube=0.158*Pr_air*Re(cbL)**0.66
Nu_D_hex=0.164*Pr_air*Re(hcP_P)**0.638 # check model for point-to-point dimension of hex
Nu_D_hex2=0.164*Pr_air*Re(hcP_P2)**0.638 # check model for point-to-point dimension of hex


h_cyl=Nu_D_cyl*k_air/cyD
h_cb=Nu_D_cube*k_air/cbL
h_hex=Nu_D_hex*k_air/hcP_P
h_hex2=Nu_D_hex2*k_air/hcP_P2

# for 3m/s: 72, 70.4, 55, 61.3, 76.2, 73.8

biot_alcyl = h_cyl*cyD/k_al
biot_alcb = h_cb*cbL/k_al
biot_alhc = h_hex*hcP_P/k_al

biot_brcyl = h_cyl*cyD/k_brass
biot_brcb = h_cb*cbL/k_brass
biot_brhc=h_hex2*hcP_P2/k_brass

t=np.linspace(0,800)
fig=plt.figure()
ax=fig.add_subplot()
# m/s
ax.plot(t,theta(rho_al,cyvol,cp_al,h_cyl[0],cyAs,t,Ti=76.2),
        label="Al cylinder",linestyle='dashed')
ax.plot(t,theta(rho_al,cbvol,cp_al,h_cb[0],cbAs,t,Ti=76.2),
        label="Al cube",linestyle='dotted')
ax.plot(t,theta(rho_al,hcvol,cp_al,h_hex[0],hcAs,t,Ti=63.6),
        label="Al hex cylinder",linestyle='dashdot')

ax.plot(t,theta(rho_brass,cyvol,cp_brass,h_cyl[0],cyAs,t,Ti=65.3),
        label="Brass cylinder",linestyle='dashed', linewidth=2.2)
ax.plot(t,theta(rho_brass,cbvol,cp_brass,h_cb[0],cbAs,t,Ti=76.2),
        label="Brass cube",linestyle='dotted', linewidth=2.2)
ax.plot(t,theta(rho_brass,hcvol2,cp_brass,h_hex2[0],hcAs,t,Ti=67.8),
        label="Brass hex cylinder",linestyle='dashdot',linewidth=2.2)

ax.set_title("Theoretical Forced Convection Cooling 1m/s Flow Rate")
ax.set_xlabel("Time (s)")
ax.set_ylabel("Temperature (C)")
ax.legend()
ax.set_ylim(35,80)
plt.show()

print("RE")
print(Re(cyD),"cyl")
print(Re(cbL),"cube")
print(Re(hcP_P),"cyl al")
print(Re(hcP_P2),"cyl,br")

print("NU")
print(Nu_D_cyl,"cyl")
print(Nu_D_cube,"cube")
print(Nu_D_hex,"cyl al")
print(Re(hcP_P2),"cyl,br")

print("H")
print(h_cyl,"cyl")
print(h_cb,"cube")
print(h_hex,"cyl al")
print(h_hex2,"cyl,br")

print("BIOT")
print(biot_alcyl,"al cyl")
print(biot_alcb,"al cube")
print(biot_alhc,"al hex")
print(biot_brcyl,"br cyl")
print(biot_brcb,"br cube")
print(biot_brhc,"br hex")


# In[126]:


t=np.linspace(0,800)
fig=plt.figure()
ax=fig.add_subplot()
# for 3m/s: 72, 70.4, 55, 61.3, 76.2, 73.8
ax.plot(t,theta(rho_al,cyvol,cp_al,h_cyl[1],cyAs,t,Ti=72),
        label="Al cylinder",linestyle='dashed')
ax.plot(t,theta(rho_al,cbvol,cp_al,h_cb[1],cbAs,t,Ti=70.4),
        label="Al cube",linestyle='dotted')
ax.plot(t,theta(rho_al,hcvol,cp_al,h_hex[1],hcAs,t,Ti=55),
        label="Al hex cylinder",linestyle='dashdot')

ax.plot(t,theta(rho_brass,cyvol,cp_brass,h_cyl[1],cyAs,t,Ti=61.3),
        label="Brass cylinder",linestyle='dashed', linewidth=2.2)
ax.plot(t,theta(rho_brass,cbvol,cp_brass,h_cb[1],cbAs,t,Ti=76.2),
        label="Brass cube",linestyle='dotted', linewidth=2.2)
ax.plot(t,theta(rho_brass,hcvol2,cp_brass,h_hex2[1],hcAs,t,Ti=73.8),
        label="Brass hex cylinder",linestyle='dashdot',linewidth=2.2)

ax.set_title("Theoretical Forced Convection Cooling 3m/s Flow Rate")
ax.set_xlabel("Time (s)")
ax.set_ylabel("Temperature (C)")
ax.legend()
ax.set_ylim(35,80)
plt.show()


# In[128]:


t=np.linspace(0,800)
fig=plt.figure()
ax=fig.add_subplot()
# for 3m/s: 72, 70.4, 55, 61.3, 76.2, 73.8
ax.plot(t,theta(rho_al,cyvol,cp_al,h_cyl[2],cyAs,t,Ti=64),
        label="Al cylinder",linestyle='dashed')
ax.plot(t,theta(rho_al,cbvol,cp_al,h_cb[2],cbAs,t,Ti=69),
        label="Al cube",linestyle='dotted')
ax.plot(t,theta(rho_al,hcvol,cp_al,h_hex[2],hcAs,t,Ti=64),
        label="Al hex cylinder",linestyle='dashdot')

ax.plot(t,theta(rho_brass,cyvol,cp_brass,h_cyl[2],cyAs,t,Ti=74.2),
        label="Brass cylinder",linestyle='dashed', linewidth=2.2)
ax.plot(t,theta(rho_brass,cbvol,cp_brass,h_cb[2],cbAs,t,Ti=76.2),
        label="Brass cube",linestyle='dotted', linewidth=2.2)
ax.plot(t,theta(rho_brass,hcvol2,cp_brass,h_hex2[2],hcAs,t,Ti=65),
        label="Brass hex cylinder",linestyle='dashdot',linewidth=2.2)

ax.set_title("Theoretical Forced Convection Cooling 5m/s Flow Rate")
ax.set_xlabel("Time (s)")
ax.set_ylabel("Temperature (C)")
ax.legend()
ax.set_ylim(35,80)
plt.show()


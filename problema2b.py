#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 02:04:18 2019

@author: paz
"""


import numpy as np
import math as mt
import matplotlib.pyplot as plt

plt.ion()



dim=""

#term=np.loadtxt('termalizacion2b')
datos=np.loadtxt('2b')


###################################################################
##
##                Correlación en H y M
##
###################################################################	
"""
#a=len(term[:,0])
#for  l in range(0,8):
l=2
m=[]
h=[]
for i in range(1,len(term[0,:]),2):
	m.append(term[l,i])
	h.append(term[l,i+1])
#m=np.arrray(m_i)
#h=np.arrray(h_i)

k_max=len(m)
nat=np.zeros(k_max)
for i in range(0,len(nat)):
	nat[i]=i+1



mm_k=np.zeros(k_max//20)
rho=np.zeros(k_max//20)

hh_k=np.zeros(k_max//20)
rho_h=np.zeros(k_max//20)

k=len(mm_k)

for i in range(0,k):
	for j in range(0,k-i):
		mm_k[i]=mm_k[i]+m[j]*m[j+i]/(k-i)    #      Calcula <M_i*M_{i+k}>
		hh_k[i]=hh_k[i]+h[j]*h[j+i]/(k-i)	#      Calcula <H_i*H_{i+k}>	
	
m_med=np.mean(m)
h_med=np.mean(h)
m_cuad_med=np.mean(np.power(m,2))
h_cuad_med=np.mean(np.power(h,2))

for i in range(0,len(mm_k)):
	rho[i]=(mm_k[i]-m_med**2)/(m_cuad_med-m_med**2)      #  (<M_i*M_{i+k}>-<M_i>**2)/(<M_i**2>-<M_i>**2)
	rho_h[i]=(hh_k[i]-h_med**2)/(h_cuad_med-h_med**2)



plt.figure()
line1=plt.plot(rho,label='C(M) J$^*={}$'.format('%.2f' %term[l,0]))
line2=plt.plot(rho_h,label='C(H) J$^*={}$'.format('%.2f' %term[l,0]))
plt.xlabel('k')
plt.ylabel('$\\rho(k)$')
#plt.xlim([-100,3000])
plt.ylim([-0.5,1.2])
plt.legend()
plt.show()

"""
###################################################################
##
##                B vs M
##
###################################################################


T=datos[:,0]
J=datos[:,1]
H=datos[:,2]
M=datos[:,3]
H_cuad=datos[:,4]
M_cuad=datos[:,5]
xi=M_cuad-M*M
cv=H_cuad-H*H	

indice_optimo=np.where(xi == max(xi))[0][0]
Tc=T[indice_optimo]

indice_optimo2=np.where(cv == max(cv))[0][0]
Tc2=T[indice_optimo2]

"""
plt.figure(2)
line1 = plt.scatter(J,M,label='$<M>$')
#line2 = plt.plot(J,np.tanh(J),label='$\\tanh(B*)$',color='r')
plt.xlabel('$J^*$')
plt.ylabel('$M$')
#plt.xlim([-10,10])
plt.legend()
plt.show()
"""


plt.figure(3)
line1 = plt.scatter(T,xi,label='$<M^2>-<M>^2$')
line2 = plt.plot([Tc2,Tc2],np.linspace(min(xi),max(xi),2),color='r',label='$T_c={}$  '.format('%.4f' %Tc2)+dim)
plt.xlabel('T')
plt.ylabel('$\\chi$')
#plt.xlim([-10,10])
#plt.xscale('log')
#plt.yscale('log')
plt.legend()
plt.show()


plt.figure(4)
line1 = plt.scatter(T,M,label='Magnetización '+dim)
plt.xlabel('T')
plt.ylabel('$M$')
#plt.xlim([-10,10])
plt.legend()
plt.show()






plt.figure(5)
line1 = plt.scatter(T,cv,label='$<H^2>-<H>^2$')
line2 = plt.plot([Tc2,Tc2],np.linspace(min(cv),max(cv),2),color='r',label='$T_c={}$  '.format('%.4f' %Tc2)+dim)
plt.xlabel('T')
plt.ylabel('Cv')
#plt.xlim([-10,10])
plt.legend()
plt.show()


plt.figure(6)
line1 = plt.scatter(T,(H),label='Energia '+dim)
plt.xlabel('T')
plt.ylabel('$H$')
#plt.xlim([-10,10])
plt.legend()
plt.show()

"""
plt.figure(7)
line1 = plt.scatter(np.log(abs(T)),np.log(cv),label='$<H^2>-<H>^2$')
plt.xlabel('T')
plt.ylabel('Cv')
#plt.xlim([-10,10])
plt.legend()
plt.show()


plt.figure(8)
line1 = plt.scatter(np.log(T),np.log(xi),label='$<M^2>-<M>^2$')
plt.xlabel('T')
plt.ylabel('$\\chi$')
#plt.xlim([-10,10])
#plt.xscale('log')
#plt.yscale('log')
plt.legend()
plt.show()
"""

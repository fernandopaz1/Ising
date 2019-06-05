#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 02:04:18 2019

@author: paz
"""


import numpy as np
import math as mt
import matplotlib.pyplot as plt





datos=np.loadtxt('1a_dim10')
correl=np.loadtxt('correlacion2a')

T=datos[:,0]
B=datos[:,1]
M=datos[:,3]
M_cuad=datos[:,5]


###################################################################
##
##                Correlación en H y M
##
###################################################################	





k_max=len(correl[0,:])
nat=np.zeros(len(correl[0,1:k_max]))
for i in range(0,len(nat)):
	nat[i]=i+1

"""
plt.figure(2)
for i in range(0,len(correl[:,0])):
	plt.plot(nat,correl[i,1:k_max],label='$<M>$ con $B^*={}$'.format('%.2f' %correl[i,0]))
	line2=plt.plot(nat,np.ones(len(nat))*np.tanh(correl[i,0]))#,label='$\\tanh(B*)$')
plt.legend(loc=1)
plt.xlabel('k pasos')
plt.ylabel('<M>')
plt.show()

"""
"""
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
m_cuad_med=np.mean(m**2)
h_cuad_med=np.mean(h**2)
print(m_med**2, m_cuad_med)
print(h_med**2, h_cuad_med)

nat=np.zeros(len(mm_k))
for i in range(0,len(mm_k)):
	nat[i]=i+1
	rho[i]=(mm_k[i]-m_med**2)/(m_cuad_med-m_med**2)      #  (<M_i*M_{i+k}>-<M_i>**2)/(<M_i**2>-<M_i>**2)
	rho_h[i]=(hh_k[i]-h_med**2)/(h_cuad_med-h_med**2)



plt.figure(1)
line1=plt.plot(nat,rho,label='C(M)')
line2=plt.plot(nat,rho_h,label='C(H)')
plt.xlabel('k')
plt.ylabel('$\\rho(k)$')
plt.legend()
plt.show()

"""
###################################################################
##
##                B vs M
##
###################################################################	

plt.figure(2)
line1 = plt.scatter(B,M,label='$<M>$')
line2 = plt.plot(B,np.tanh(B),label='$\\tanh(B*)$',color='r')
plt.xlabel('$B^*$')
plt.ylabel('$M$')
#plt.xlim([-10,10])
plt.legend()
plt.show()



plt.figure(3)
line1 = plt.scatter(T,np.abs(M_cuad-M*M),label='$<M^2>-<M>^2$')
plt.xlabel('T')
plt.ylabel('$\\chi$')
#plt.xlim([-10,10])
plt.legend()
plt.show()


plt.figure(4)
line1 = plt.scatter(T,np.abs(M),label='Magnetización')
plt.xlabel('T')
plt.ylabel('$M$')
#plt.xlim([-10,10])
plt.legend()
plt.show()


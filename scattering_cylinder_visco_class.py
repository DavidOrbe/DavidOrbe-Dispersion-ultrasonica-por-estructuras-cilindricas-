# -*- coding: utf-8 -*-
"""
Created on Fri Feb 18 20:51:38 2022

@author: david
"""
'''
En este codigo se presentan las clases que resuelven diferentes propiedades
en el problema de dispersion de estructuras cilindricas
'''

import scipy.special as sp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

P= lambda x,y: np.power(x,y) #Potencia de un valor x a un valor y

class medio_dispersor():
    
    '''
    Esta clase a traves de el material y medio seleccionado nos regresa sus
    propiedades asi como los numeros de onda del medio 
    '''
    
    def __init__(self,material,medio):
        self.material=material 
        self.medio=medio 
        
    def eleccion_material(self):
        archivo = './DatosPolimeros.xlsx'
        dispersor=pd.read_excel(archivo,sheet_name=str(self.material))
        return dispersor
    
    def eleccion_medio_numero_onda(self):
        
        '''ks numero de onda de corte (shear) y kc numero de onda longitudinal
        (compresión) en el fluido '''
        
        archivo = './DatosPolimeros.xlsx' 
        medio= pd.read_excel(archivo, sheet_name=str(self.medio))
        self.c,self.rho_l,self.mu,self.mu_b=medio.values.T
        w=(self.ka*self.c/self.a)
        kc=(w/self.c)*(1+((1j*w*self.mu)/(2*self.rho_l*P(self.c,2)))*(4/3+(self.mu_b/self.mu)))
        ks=(1+1j)*np.sqrt((w*self.rho_l)/(2*self.mu))
        self.w=w
        self.kc=kc
        self.ks=ks
        return kc, ks

class elastico(medio_dispersor):
    
    '''
    Esta clase nos devuelve las constantes de lame y los numeros de onda para
    un dispersor elastico
    '''
    
    def __init__(self,material,medio,ka,a):
        
        medio_dispersor.__init__(self, material, medio)
        
        self.ka=ka
        self.a=a
        
    def lame_elastico(self):    
       
        self.rho_s,nu,E=self.eleccion_material().values.T  
        M_corte=E/(2*(1+nu))
        Lambda=(nu*E)/((1+nu)*(1-2*nu))
        Mu=M_corte
        return Lambda, Mu
    
    def numeros_onda_elastico(self):
        ''' Kc es el numero de onda de compresion y Ks el de corte en el
        solido elastico'''
        
        Lambda,Mu=self.lame_elastico()
        self.eleccion_medio_numero_onda()
      
        Kc=self.w/np.sqrt((Lambda+2*Mu)/self.rho_s)
        Ks=self.w/np.sqrt(Mu/self.rho_s)
        self.Kc=Kc
        self.Ks=Ks 
        
        return Kc, Ks
    
class viscoelastico(medio_dispersor):
    
    '''
    Esta clase nos devuelve las constantes de lame y los numeros de onda para
    un dispersor viscoelástico, asimismo permite graficar el comportamiento a
    frecuencia de los modulos de almacenamiento y de perdidas, adicionalmente 
    se puede graficar el factor de perdida
    '''
    
    def __init__(self,material,medio,ka,a):
        
        medio_dispersor.__init__(self, material, medio)
        
        self.ka=ka
        self.a=a
        
    def lame_viscoelasticos(self):
            
        G0,G_inf,alpha,beta,tau,self.rho_s,nu=self.eleccion_material().values.T
        self.eleccion_medio_numero_onda()
        w=self.w
        k=np.arctan((P(w,alpha)*P(tau,alpha)*np.sin(alpha*np.pi/2))/(1+P(w,alpha)\
          *P(tau,alpha)*np.cos(alpha*np.pi/2)))  
        
        self.k=k   #puede ser una opcion para llamar los atributos sin tanto self                                                
                      
        M_almacenamiento=G_inf+((G0-G_inf)*np.cos(beta*k))/(P(1+2*P(w,alpha)\
        *P(tau,alpha)*np.cos(alpha*np.pi/2)+P(w,2*alpha)*P(tau,2*alpha),beta/2))
        self.M_almacenamiento=M_almacenamiento
                
        M_perdidas=((G_inf-G0)*np.sin(beta*k))/(P(1+2*P(w,alpha)*P(tau,alpha)\
                    *np.cos(alpha*np.pi/2)+P(w,2*alpha)*P(tau,2*alpha),beta/2))
        self.M_perdidas=M_perdidas
                    
        M_corte=M_almacenamiento+1j*M_perdidas
        self.M_corte=M_corte
        
        T_perdidas=M_perdidas/M_almacenamiento
        self.T_perdidas=T_perdidas
        
        Lambda=(2*nu*M_corte)/(1-2*nu)
        self.Lambda=Lambda
        
        Mu=M_corte
        self.Mu=Mu
        
        return Lambda, Mu
        
    def numeros_onda_viscoelastico(self):
        ''' Kc es el numero de onda de corte y Ks el de compresion en el
        solido viscoelastico'''
        
        Lambda,Mu=self.lame_viscoelasticos()
        # G0,G_inf,alpha,beta,tau,rho_s,nu=self.eleccion_material().values.T
        Kc=self.w/np.sqrt((Lambda+2*Mu)/self.rho_s)
        Ks=self.w/np.sqrt(Mu/self.rho_s)
        self.Kc=Kc
        self.Ks=Ks 
        return Kc, Ks
    
    def plot_modulus(self):
        plt.figure()
        plt.loglog(self.ka,self.M_almacenamiento, linestyle='dashed', label='G\'')    
        plt.loglog(self.ka,self.M_perdidas, linestyle='dotted', label='G\'\'')
        # plt.loglog(self.ka,self.T_perdidas)
      
        plt.vlines(x=1e-1, ymin=1, ymax=1e10,linestyle='dashed', color='black')
        plt.vlines(x=10, ymin=1, ymax=1e10,linestyle='dashed', color='black')
        plt.xlim(1e-11,1e13)
        plt.ylim(1e4,1e10)
        # plt.title(str(self.material) 
        plt.xlabel('Frecuencia adimensional, ka')
        plt.ylabel('Módulo de Corte [Pa]')
        plt.legend(loc='lower right')
        plt.annotate('ka=10', xy=(50,1e5))
        plt.annotate('ka=0.1', xy=(0.00005,1e5))
        
    def plot_factor_perdida(self):
        plt.figure()
        plt.loglog(self.ka,self.T_perdidas, linestyle='dashed', label=r'$tan(\delta)$')
        plt.vlines(x=1e-1, ymin=0, ymax=2,linestyle='dashed', color='black')
        plt.vlines(x=10, ymin=0, ymax=2,linestyle='dashed', color='black')
        plt.xlim(1e-11,1e13)
        plt.ylim(1e-4,1.5)
        # plt.title(str(self.material) 
        plt.xlabel('Frecuencia adimensional, ka')
        plt.ylabel('Tangente de pérdidas')
        plt.legend(loc='lower right')
        plt.annotate('ka=10', xy=(50,1.5e-4))
        plt.annotate('ka=0.1', xy=(0.00005,1.5e-4))

class cilindro_visco_elastico(elastico,viscoelastico):
    
    '''
    Esta clase resuleve el problema de dispersion para un cilindro (visco)elastico
    y grafica sus patrones de dispersion
    '''
    
    
    def __init__(self,material,medio,ka,a):
        
        medio_dispersor.__init__(self, material, medio)
        elastico.__init__(self,material,medio,ka,a)
        viscoelastico.__init__(self,material,medio,ka,a)
        self.n=np.array([np.arange(1,40)])
        self.theta=np.array([np.linspace(0,2*np.pi,200)])
        
        
    def sistema_ecuaciones(self):
        '''
        Definimos los coeficientes a encontrar
        '''
        
        if 'POLIMERO' in self.material:
            Kc,Ks=self.numeros_onda_viscoelastico()
            Lambda, Mu=self.lame_viscoelasticos()
        else:
            Kc,Ks=self.numeros_onda_elastico()
            Lambda, Mu =self.lame_elastico()
        
        kc,ks=self.eleccion_medio_numero_onda()
        a=self.a
        w=self.w
        phi=1
        n=self.n
        rho_l=self.rho_l 
        mu=self.mu
        
        
        xkc=kc*a
        xks=ks*a
        xKc=Kc*a
        xKs=Ks*a
        self.xkc=xkc
        self.xks=xks
        self.xKc=xKc
        self.xKs=xKs
        # #termino cero 
        A01=-kc*sp.h1vp(0,xkc,1)
        C01=1j*w*Kc*sp.jvp(0,xKc,1)
        I01=phi*kc*sp.jvp(0,xkc,1)
        A03=(1j*w*rho_l-2*mu*P(kc,2))*sp.hankel1(0,xkc)-2*mu*P(kc,2)*sp.h1vp(0,xkc,2)
        C03=Lambda*P(Kc,2)*sp.jv(0,xKc)-2*Mu*P(Kc,2)*sp.jvp(0,xKc,2)
        I03=-phi*((1j*w*rho_l-2*mu*P(kc,2))*sp.jv(0,xkc)-2*mu*P(kc,2)*sp.jvp(0,xkc,2))
        self.A01=A01
        self.C01=C01
        self.I01=I01
        self.A03=A03
        self.C03=C03
        self.I03=I03
        
        # terminos mayores a cero
        A1=-kc*sp.h1vp(n,xkc,1)
        B1=(n/a)*sp.hankel1(n,xks)
        C1=1j*w*Kc*sp.jvp(n,xKc,1)
        D1=((1j*w*n)/a)*sp.jv(n,xKs)
        I1=2*phi*P(1j,n)*kc*sp.jvp(n,xkc,1)
        
        self.A1=A1
        self.B1=B1
        self.C1=C1
        self.D1=D1
        self.I1=I1
        
        A2=(n/a)*sp.hankel1(n,xkc)
        B2=-ks*sp.h1vp(n,xks,1)
        C2=((-1j*w*n)/a)*sp.jv(n,xKc)
        D2=-1j*w*Ks*sp.jvp(n,xKs,1)
        I2=((-2*phi*n)/a)*P(1j,n)*sp.jv(n,xkc)    
        
        self.A2=A2
        self.B2=B2
        self.C2=C2
        self.D2=D2
        self.I2=I2
            
        A3=(1j*w*rho_l-2*mu*P(kc,2))*sp.hankel1(n,xkc)-2*mu*P(kc,2)*sp.h1vp(n,xkc,2)
        B3=((-2*mu*n)/P(a,2))*(sp.hankel1(n,xks)-xks*sp.h1vp(n,xks,1))    
        C3=-2*Mu*P(Kc,2)*sp.jvp(n,xKc,2)+Lambda*P(Kc,2)*sp.jv(n,xKc)
        D3=((2*Mu*n)/P(a,2))*(sp.jv(n,xKs)-xKs*sp.jvp(n,xKs,1))
        I3=(-2*phi*P(1j,n))*((1j*w*rho_l-2*mu*P(kc,2))*sp.jv(n,xkc)-2*mu*P(kc,2)*sp.jvp(n,xkc,2))
        
        self.A3=A3
        self.B3=B3
        self.C3=C3
        self.D3=D3
        self.I3=I3
            
        A4=(-2*mu*n)*(sp.hankel1(n,xkc)-xkc*sp.h1vp(n,xkc,1))
        B4=-mu*(P(n,2)*sp.hankel1(n,xks)+P(xks,2)*sp.h1vp(n,xks,2)-xks*sp.h1vp(n,xks,1))
        C4=-2*Mu*n*(sp.jv(n,xKc)-xKc*sp.jvp(n,xKc,1))
        D4=Mu*(P(n,2)*sp.jv(n,xKs)+P(xKs,2)*sp.jvp(n,xKs,2)-xKs*sp.jvp(n,xKs,1))
        I4=4*phi*mu*n*P(1j,n)*(sp.jv(n,xkc)-xkc*sp.jvp(n,xkc,1))
        
        self.A4=A4
        self.B4=B4
        self.C4=C4
        self.D4=D4
        self.I4=I4
             
        
        #construimos las matrices para resolver el sistema de ecuaciones
        
        
        term_cero=np.array([[A01,C01],[A03,C03]]).reshape(2,2)
        ind_cero=np.array([I01,I03]).reshape(2,1)
        coef_cero=np.linalg.solve(term_cero,ind_cero)
        
        coef_sist_ecu=np.array([[A1,B1,C1,D1],[A2,B2,C2,D2],[A3,B3,C3,D3],[A4,B4,C4,D4]]).reshape(4,4,len(n.T)).transpose(2,0,1)
        self.coef_sist_ecu=coef_sist_ecu
        terminos_ind=np.array([I1,I2,I3,I4]).transpose(2,0,1)
        self.terminos_ind=terminos_ind
        sol_coef=np.linalg.solve(coef_sist_ecu,terminos_ind)
        self.sol_coef=sol_coef
        
        coeficientes_An=np.array([np.append(coef_cero[0,0],sol_coef[:,0])])
        self.coeficientes_An=coeficientes_An
        
        ########error####
        self.cambio_coef=cambio_coef(coeficientes_An)
        
        
    def presion_dispersada(self):
        plot_presion_dispersada(self)
        

class cilindro_compuesto(elastico,viscoelastico):
    
    '''
    Esta clase resuleve el problema de dispersion para un cilindro compuesto
    y grafica sus patrones de dispersion
    '''
    
    def __init__(self,material,material2,medio,ka,a,b_a):
        
        # medio_dispersor.__init__(self, material, medio)
        # elastico.__init__(self,material,medio,ka,a)
        # viscoelastico.__init__(self,material,medio,ka,a)
        self.material=material
        self.material2=material2
        self.medio=medio
        self.ka=ka
        self.a=a
        self.b_a=b_a
        self.b=self.b_a*self.a
        self.n=np.array([np.arange(1,40)])      #maybe we can put in the function data
        self.theta=np.array([np.linspace(0,2*np.pi,200)])    #Maybe we can put in the function data np.array(np.linspace(0,2*mt.pi,10)) np.array([np.pi/2,np.pi])      
        
    def cubierta(self):
        '''
        
        '''
        medio_dispersor.__init__(self, self.material, self.medio)
        elastico.__init__(self,self.material,self.medio,self.ka,self.a)
        viscoelastico.__init__(self,self.material,self.medio,self.ka,self.a)
        
        if 'POLIMERO' in self.material:
            Kc1,Ks1=self.numeros_onda_viscoelastico()
            Lambda1,Mu1 =self.lame_viscoelasticos()
        else:
            Kc1,Ks1=self.numeros_onda_elastico()
            Lambda1, Mu1 =self.lame_elastico()
            
        return Kc1,Ks1,Lambda1,Mu1
            
    def nucleo(self):
        '''
        Definimos los coeficientes a encontrar
        '''
        medio_dispersor.__init__(self, self.material2, self.medio)
        elastico.__init__(self,self.material2,self.medio,self.ka,self.a)
        viscoelastico.__init__(self,self.material2,self.medio,self.ka,self.a)
        
        if 'POLIMERO' in self.material2:
            Kc2,Ks2=self.numeros_onda_viscoelastico()
            Lambda2,Mu2 =self.lame_viscoelasticos()
        else:
            Kc2,Ks2=self.numeros_onda_elastico()
            Lambda2, Mu2 =self.lame_elastico()
        
        return Kc2,Ks2,Lambda2,Mu2
    
    def sistema_ecuaciones(self):
        
        
        kc,ks=self.eleccion_medio_numero_onda()
        Kc1,Ks1,Lambda1,Mu1=self.cubierta()
        Kc2,Ks2,Lambda2,Mu2=self.nucleo()
        a=self.a
        b=self.b
        w=self.w
        phi=1
        n=self.n
        rho_l=self.rho_l 
        mu=self.mu
        
        xkc=kc*a
        xks=ks*a
        xKc1A=Kc1*a
        xKc1B=Kc1*b
        xKs1A=Ks1*a
        xKs1B=Ks1*b
        
        xKc2A=Kc2*a
        xKc2B=Kc2*b
        xKs2A=Ks2*a
        xKs2B=Ks2*b
        
        self.xkc=xkc
        self.xks=xks
        self.xKc1A=xKc1A
        self.xKc1B=xKc1B
        self.xKs1A=xKs1A
        self.xKs1B=xKs1B
        
        self.xKc2A=xKc2A
        self.xKc2B=xKc2B
        self.xKs2A=xKs2A
        self.xKs2B=xKs2B
        
        #Termino cero
            
        A01=-xkc*sp.h1vp(0,xkc,1) 
        C01=1j*w*xKc1A*sp.h1vp(0,xKc1A,1)
        D01=1j*w*xKc1A*sp.h2vp(0,xKc1A,1)
        G01=np.zeros(1)
        I01=xkc*phi*sp.jvp(0,xkc,1)
        
        A03=np.zeros(1)
        C03=xKc1B*sp.h1vp(0,xKc1B,1)
        D03=xKc1B*sp.h2vp(0,xKc1B,1)
        G03=-xKc2B*sp.h1vp(0,xKc2B,1)
        I03=np.zeros(1)
        
        A05=(1j*w*rho_l*P(a,2)-2*mu*P(xkc,2))*sp.hankel1(0,xkc)-2*mu*P(xkc,2)\
            *sp.h2vp(0,xkc,2)
        C05=P(xKc1A,2)*(Lambda1*sp.hankel1(0,xKc1A)-2*Mu1*sp.h1vp(0,xKc1A,2))
        D05=P(xKc1A,2)*(Lambda1*sp.hankel2(0,xKc1A)-2*Mu1*sp.h2vp(0,xKc1A,2))  
        G05=np.zeros(1)
        I05=-((1j*w*rho_l*P(a,2)-2*mu*P(xkc,2))*sp.jv(0,xkc)-2*mu*P(xkc,2)\
              *sp.jvp(0,xkc,2))*phi
        
        A07=np.zeros(1)
        C07=P(xKc1B,2)*(2*Mu1*sp.h1vp(0,xKc1B,2)-Lambda1*sp.hankel1(0,xKc1B))
        D07=P(xKc1B,2)*(2*Mu1*sp.h2vp(0,xKc1B,2)-Lambda1*sp.hankel2(0,xKc1B))
        G07=-P(xKc2B,2)*(2*Mu2*sp.jvp(0,xKc2B,2)-Lambda2*sp.jv(0,xKc2B))
        I07=np.zeros(1)
        
        self.A01=A01
        self.C01=C01
        self.D01=D01
        self.G01=G01
        self.I01=I01
        
        self.A03=A03
        self.C03=C03
        self.D03=D03
        self.G03=G03
        self.I03=I03
        
        self.A05=A05
        self.C05=C05
        self.D05=D05
        self.G05=G05
        self.I05=I05
        
        self.A07=A07
        self.C07=C07
        self.D07=D07
        self.G07=G07
        self.I07=I07
       

    # def coefficients_a    
        
        A1=-xkc*sp.h1vp(n,xkc)
        B1=n*sp.hankel1(n,xks)
        C1=1j*w*xKc1A*sp.h1vp(n,xKc1A)
        D1=(1j*w*xKc1A)*(sp.h2vp(n,xKc1A))
        E1=1j*w*n*sp.hankel1(n,xKs1A)
        F1=1j*w*n*sp.hankel2(n,xKs1A)
        G1=np.array([np.zeros(len(self.n.T))])
        K1=np.array([np.zeros(len(self.n.T))])
        I1=2*xkc*phi*P(1j,n)*sp.jvp(n,xkc,1)
        
        A2=n*sp.hankel1(n,xkc)
        B2=-xks*sp.h1vp(n,xks)
        C2=-1j*w*n*sp.hankel1(n,xKc1A)
        D2=-1j*w*n*sp.hankel2(n,xKc1A)
        E2=-1j*w*xKs1A*sp.h1vp(n,xKs1A)
        F2=-1j*w*xKs1A*sp.h2vp(n,xKs1A)
        G2=np.array([np.zeros(len(self.n.T))])
        K2=np.array([np.zeros(len(self.n.T))])
        I2=-2*phi*n*P(1j,n)*sp.jv(n,xkc)
        
        A3=np.array([np.zeros(len(self.n.T))])
        B3=np.array([np.zeros(len(self.n.T))])
        C3=xKc1B*sp.h1vp(n,xKc1B)
        D3=xKc1B*sp.h2vp(n,xKc1B)
        E3=n*sp.hankel1(n,xKs1B)
        F3=n*sp.hankel2(n,xKs1B)
        G3=-xKc2B*sp.jvp(n,xKc2B,1)
        K3=-n*sp.jv(n,xKs2B)
        I3=np.array([np.zeros(len(self.n.T))])
        
        
        A4=np.array([np.zeros(len(self.n.T))])
        B4=np.array([np.zeros(len(self.n.T))])
        C4=-n*sp.hankel1(n,xKc1B)
        D4=-n*sp.hankel2(n,xKc1B)
        E4=-xKs1B*sp.h1vp(n,xKs1B,1)
        F4=-xKs1B*sp.h2vp(n,xKs1B,1)
        G4=n*sp.jv(n,xKc2B)
        K4=xKs2B*sp.jvp(n,xKs2B,1)
        I4=np.array([np.zeros(len(self.n.T))])
                 
        A5=(1j*w*rho_l*P(a,2)-2*mu*P(xkc,2))*sp.hankel1(n,xkc)-2*mu*P(xkc,2)\
            *sp.h2vp(n,xkc,2)
        B5=(-2*mu*n)*(sp.hankel1(n,xks)-xks*sp.h1vp(n,xks,1))
        C5=-P(xKc1A,2)*(2*Mu1*sp.h1vp(n,xKc1A,2)-Lambda1*sp.hankel1(n,xKc1A))
        D5=-P(xKc1A,2)*(2*Mu1*sp.h2vp(n,xKc1A,2)-Lambda1*sp.hankel2(n,xKc1A))
        E5=2*Mu1*n*(sp.hankel1(n,xKs1A)-xKs1A*sp.h1vp(n,xKs1A,1))
        F5=2*Mu1*n*(sp.hankel2(n,xKs1A)-xKs1A*sp.h2vp(n,xKs1A,1))
        G5=np.array([np.zeros(len(self.n.T))])
        K5=np.array([np.zeros(len(self.n.T))])
        I5=-((1j*w*rho_l*P(a,2)-2*mu*P(xkc,2))*sp.jv(n,xkc)-2*mu*P(xkc,2)\
             *sp.jvp(n,xkc,2))*2*phi*P(1j,n)
        
        A6=-2*mu*n*(sp.hankel1(n,xkc)-xkc*sp.h1vp(n,xkc))
        B6=-mu*(P(n,2)*sp.hankel1(n,xks)+P(xks,2)*sp.h1vp(n,xks,2)\
                -xks*sp.h1vp(n,xks))
        C6=-2*Mu1*n*(sp.hankel1(n,xKc1A)-xKc1A*sp.h1vp(n,xKc1A,1))
        D6=-2*Mu1*n*(sp.hankel2(n,xKs1A)-xKs1A*sp.h2vp(n,xKs1A,1))
        E6=Mu1*(P(n,2)*sp.hankel1(n,xKs1A)+P(xKs1A,2)*sp.h1vp(n,xKs1A,2)\
                -xKs1A*sp.h1vp(n,xKs1A))
        F6=Mu1*(P(n,2)*sp.hankel2(n,xKs1A)+P(xKs1A,2)*sp.h2vp(n,xKs1A,2)\
                -xKs1A*sp.h2vp(n,xKs1A))
        G6=np.array([np.zeros(len(self.n.T))])
        K6=np.array([np.zeros(len(self.n.T))])
        I6=4*mu*n*phi*P(1j,n)*(sp.jv(n,xkc)-xkc*sp.jvp(n,xkc))
        
        A7=np.array([np.zeros(len(self.n.T))])
        B7=np.array([np.zeros(len(self.n.T))])
        C7=P(xKc1B,2)*(-Lambda1*sp.hankel1(n,xKc1B)+2*Mu1*sp.h1vp(n,xKc1B,2))
        D7=P(xKc1B,2)*(-Lambda1*sp.hankel2(n,xKc1B)+2*Mu1*sp.h2vp(n,xKc1B,2))    
        E7=-2*Mu1*n*(sp.hankel1(n,xKs1B)-xKs1B*sp.h1vp(n,xKs1B))
        F7=-2*Mu1*n*(sp.hankel2(n,xKs1B)-xKs1B*sp.h2vp(n,xKs1B))
        G7=P(xKc2B,2)*(Lambda2*sp.jv(n,xKc2B)-2*Mu2*sp.jvp(n,xKc2B,2))
        K7=2*Mu2*n*(sp.jv(n,xKs2B)-xKs2B*sp.jvp(n,xKs2B))
        I7=np.array([np.zeros(len(self.n.T))])           
        
        A8=np.array([np.zeros(len(self.n.T))])
        B8=np.array([np.zeros(len(self.n.T))])
        C8=2*Mu1*(sp.hankel1(n,xKc1B)-xKc1B*sp.h1vp(n,xKc1B))
        D8=2*Mu1*(sp.hankel2(n,xKc1B)-xKc1B*sp.h2vp(n,xKc1B))
        E8=-Mu1*(P(n,2)*sp.hankel1(n,xKs1B)+P(xKs1B,2)*sp.h1vp(n,xKs1B,2)\
                -xKs1B*sp.h1vp(n,xKs1B))
        F8=-Mu1*(P(n,2)*sp.hankel2(n,xKs1B)+P(xKs1B,2)*sp.h2vp(n,xKs1B,2)\
                -xKs1B*sp.h2vp(n,xKs1B))
        G8=-2*Mu2*n*(sp.jv(n,xKc2B)-xKc2B*sp.jvp(n,xKc2B))
        K8=Mu2*(P(n,2)*sp.jv(n,xKs2B)+P(xKs2B,2)*sp.jvp(n,xKs2B,2)\
                -xKs2B*sp.jvp(n,xKs2B))
        I8=np.array([np.zeros(len(self.n.T))])
        
        
        #construimos las matrices para resolver el sistema de ecuaciones
        
        
        term_cero=np.array([[A01,C01,D01,G01],[A03,C03,D03,G03],[A05,C05,D05,G05],[A07,C07,D07,G07]]).reshape(4,4)
        ind_cero=np.array([I01,I03,I05,I07]).reshape(4,1)
        coef_cero=np.linalg.solve(term_cero,ind_cero)
        self.term_cero=term_cero
        
        coef_sist_ecu=np.array([[A1,B1,C1,D1,E1,F1,G1,K1],[A2,B2,C2,D2,E2,F2,G2,K2],[A3,B3,C3,D3,E3,F3,G3,K3],[A4,B4,C4,D4,E4,F4,G4,K4],[A5,B5,C5,D5,E5,F5,G5,K5],[A6,B6,C6,D6,E6,F6,G6,K6],[A7,B7,C7,D7,E7,F7,G7,K7],[A8,B8,C8,D8,E8,F8,G8,K8]]).reshape(8,8,len(n.T)).transpose(2,0,1)
        self.coef_sist_ecu=coef_sist_ecu
        terminos_ind=np.array([I1,I2,I3,I4,I5,I6,I7,I8]).transpose(2,0,1)
        self.terminos_ind=terminos_ind
        sol_coef=np.linalg.solve(coef_sist_ecu,terminos_ind)
        self.sol_coef=sol_coef
        
        coeficientes_An=np.array([np.append(coef_cero[0,0],sol_coef[:,0])])
        self.coeficientes_An=coeficientes_An
        
        ########error####
        self.cambio_coef=cambio_coef(coeficientes_An)
        
        
    def presion_dispersada(self):
        plot_presion_dispersada(self)

 
def cambio_coef(coeficientes_An):
    """
    Funcion que recibe el arreglo de coeficientes An y calcula la diferencia en el cambio de coeficientes,
    es decir, se obtiene la diferencia entre el orden n y el n-1, siendo que 
    si su diferencia es menor a 1e-30 la contribucion de los siguientes ordenes n 
    sera muy pequeña
    """

    cambio=1
    orden=0
    while cambio>1e-30:
       cambio=abs(coeficientes_An[0,orden+1]-coeficientes_An[0,orden])
       orden+=1
    return orden

def plot_presion_dispersada(self):
    
    '''
    Funcion para graficar la presion dispersada 
    '''
    
    self.sistema_ecuaciones()
    n0=np.append(0,self.n)
    self.n0=n0
    mu_b=self.mu_b
    c=self.c
    theta=self.theta
       
    P_scat=(P(2/(np.pi*self.xkc),1/2)*(-1j*self.w*self.rho_l+P(self.kc,2)\
            *(mu_b+(4/3)*self.mu))*self.coeficientes_An\
            *np.exp(1j*(self.xkc-(n0*np.pi)/2-np.pi/4))*np.cos(n0*self.theta.T)).sum(axis=1)
    self.P_scat=P_scat                               
    sigma=abs(P_scat)/(self.rho_l*P(c,2))
    self.sigma=sigma
    
    self.maximo_sigma=max(self.sigma)
    self.normalizado=self.sigma/self.maximo_sigma

    plt.figure()
    plt.polar(theta.reshape(-1),self.normalizado,'r')
    
    
       
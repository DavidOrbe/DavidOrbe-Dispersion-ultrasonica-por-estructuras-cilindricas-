# -*- coding: utf-8 -*-
"""
Created on Fri Feb 18 22:27:13 2022

@author: david
"""

'''
El presente programa resuelve el problema de dispersion de una onda plana incidente
en un dispersor cilindrico circular de longitud infinita inmerso en un fluido viscoso,
dicho dispersor tiene como propiedad de ser (visco)elastico

Se desarrolla un scrip llamado 'scattering_cylinder_visco_class.py'  donde estaran todas
las clases y metodos para resolver el problema por el metodo de condiciones a la frontera
y generar los patrones de dispersion.

Finalmente, en este scrip llamado 'scattering_visco_cylinder.py' llamamos al resto de metodos
y especificamos el rango de frecuencias y radios  para obtener los patrones de dispersion

'''

import numpy as np
import scattering_cylinder_visco_class

'''
A Continuacion se presenta un menu para elegir el tipo de dispersor, el medio 
e ingresar los valores como radios, frecuencia etc.

'''
 
def convertidor_eleccion(texto):
    '''
    Funcion que recibe un caracter de numero por las opciones y les da su
    correcpondiente nombre
    '''
    
    if texto=='1':
        texto='ACERO'
    elif texto=='2':
        texto='LATON'
    elif texto=='3':
        texto='POLIMERO18'
    elif texto=='4':
        texto='POLIMERO19'
    elif texto=='5':
        texto='POLIMERO20'
    else:
        print('No esta en las opciones')
    return texto

opcion=''

while opcion!='n': 

    print('Elija el medio en el que estara inmerso el dispersor:',
          '\n1. GLICERINA \n2. AGUA \nIngrese n para salir')

    medio=input('Opcion:')
    if medio=='1':
        medio='GLICERINA'
    elif medio=='2':
        medio='AGUA'
    else:
        opcion='n'
        break
    print('\n')
    print('Elija el tipo de dispersor cilindrico:',
          '\n1. (visco)elastico \n2. Compuesto \nIngrese n para salir')

    dispersor=input('Opcion:')
    print('\n')

    if dispersor=='1':
        a=float(input('Ingrese el radio del cilindro:'))
        print('\n')
        print('Elija el tipo de material para el cilindro:',
              '\n1. ACERO \n2. LATON \n3. POLIMERO18 \n4. POLIMERO19',
              '\n5. POLIMERO20 \nIngrese n para salir')

        material_simple=convertidor_eleccion(input('Opcion:'))
        print('\n')
        if 'POLIMERO' in material_simple:
            print('Eligio un polímero como material', 
                  '\n¿Quiere ver su comportamiento a frecuencia? s/n')
            eleccion=input('Opcion:')
            print('\n')
            if eleccion=='s':
                ka=np.logspace(-11,13,100)
                cilindro_visco_elastico=scattering_cylinder_visco_class.cilindro_visco_elastico(material_simple,medio, ka, a)
                cilindro_visco_elastico.lame_viscoelasticos()
                cilindro_visco_elastico.plot_modulus()
                cilindro_visco_elastico.plot_factor_perdida()            
        elif material_simple!='ACERO' and material_simple!='LATON':
            opcion='n'
            break
            
        print('Ingrese el las frecuencias de la onda de ultrasonido:',
              '\nseparadas por comas. Ej: 1,2,3')
        ka=input('Frecuencias: ').split(',')
        print('\n')
        for i in ka:
            cilindro_visco_elastico=scattering_cylinder_visco_class.cilindro_visco_elastico(material_simple,medio, float(i), a)
            cilindro_visco_elastico.presion_dispersada()
        print('¿Quieres simular otra dispersion?: s/n')
        opcion=input('Opcion: ')

    elif dispersor=='2':
        a=float(input('Ingrese el radio externo del cilindro:'))
        
        b_a=float(input('Ingrese la razon entre el radio interno y el externo: '))
        
        print('Elija el tipo de material para el nucleo:',
              '\n1. ACERO \n2. LATON \n3. POLIMERO18 \n4. POLIMERO19',
              '\n5. POLIMERO20')
    
        compuesto_nucleo=convertidor_eleccion(input('Opcion:'))
    
        print('Elija el tipo de material para la cubierta:',
              '\n1. ACERO \n2. LATON \n3. POLIMERO18 \n4. POLIMERO19',
              '\n5. POLIMERO20')
    
        compuesto_cubierta=convertidor_eleccion(input('Opcion:'))
        
        print('Ingrese el las frecuencias de la onda de ultrasonido:',
              '\nseparadas por comas. Ej: 1,2,3')
        ka=input('Frecuencias: ').split(',')
        
        for i in ka:
            cilindro_compuesto=scattering_cylinder_visco_class.cilindro_compuesto(compuesto_cubierta,compuesto_nucleo,medio, float(i), a,b_a)
            cilindro_compuesto.presion_dispersada()
        print('¿Quieres simular otra dispersion?: s/n')
        opcion=input('Opcion: ')
        
    else:
        opcion='n'
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 11:14:05 2024

@author: ICTEA
"""

# calculate magnetic fields (in tesla, T)
# with the Biot-Savart law
# the output is points (coordinates x,y,z) and B(1,B2,B3) for the 3 cases (B1,B2,B3)
import numpy as np
import wire
import biotsavart
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


df= pd.read_csv('./efecto_pos3.csv', header=None, sep=';',names=['x', 'y', 'z', 'Bx', 'By', 'Bz', 'Btotal', 'Bx_T', 'By_T', 'Bz_T'])
df = df.drop(['Btotal','Bx_T', 'By_T', 'Bz_T'], axis=1)
print(df)

puntos = df[['x', 'y', 'z']].values.tolist()

campo = df[['Bx', 'By', 'Bz']].values.tolist()



Angulo = np.arctan2(df['y'], df['x'])


#Definicion de constantes


I_rectangular = 69
I_circulares = 82


Altura = 72
pos_espira_rectangular = np.arange(-32,33,2)
radio_cilindro = 34
radio_tapas = 26
radio_tapas2 = 20
pos_espira_circular = np.arange(-33.1,34.2,2.4)
PMTs_vertical = np.arange(-33.4,32.5,0.707)
ceros = np.zeros(len(PMTs_vertical))


Media = []
points2=[]
PMTs_top=[]
PMTs_bottom=[]
Coordenadas = []
Bperp_paredes=[]
Optimizacion=[]
Longitud=[]





#Programa principal

def Espiras(puntos, campo, Angulo):

	# rectangular loops I=1 A

	
    w1r = wire.Wire(path=wire.Wire.RectangularPath(dx=np.sqrt(radio_cilindro**2 - pos_espira_rectangular[0]**2)*2, dy=Altura), discretization_length=0.1, current=I_rectangular).Rotate(axis=(1, 0, 0), deg=-90).Translate([0,pos_espira_rectangular[0], 0]).Translate([0,0, 0.5])
    sol = biotsavart.BiotSavart(wire=w1r)

    Longitud.append(np.sqrt(radio_cilindro**2 - pos_espira_rectangular[0]**2)*4+Altura*2)




    for i in range(1, len(pos_espira_rectangular)-1):
        lado_espira_rect = np.sqrt(radio_cilindro**2 - pos_espira_rectangular[i]**2)
        
        if i ==3:
            w2r = wire.Wire(path=wire.Wire.GondolaIDPath(dx=lado_espira_rect*2, dy=Altura, lado=1.91, semieje=0.5, puerta=0.5), discretization_length=0.1, current=I_rectangular).Rotate(axis=(1, 0, 0), deg=-90).Translate([0,pos_espira_rectangular[i],0]).Translate([0,0, 0.5])
            sol.AddWire(w2r)
            Longitud.append(np.sqrt(radio_cilindro**2 - pos_espira_rectangular[i]**2)*4+Altura*2)
            continue
        
        
        if i ==8:
            w2r = wire.Wire(path=wire.Wire.GondolaIDPath(dx=lado_espira_rect*2, dy=Altura, lado=10.46, semieje=1.98, puerta=1.8), discretization_length=0.1, current=I_rectangular).Rotate(axis=(1, 0, 0), deg=-90).Translate([0,pos_espira_rectangular[i],0]).Translate([0,0, 0.5])
            sol.AddWire(w2r)
            Longitud.append(np.sqrt(radio_cilindro**2 - pos_espira_rectangular[i]**2)*4+Altura*2)
            continue
        
        if i ==16:
            w2r = wire.Wire(path=wire.Wire.CompuertaPath(dx=lado_espira_rect*2, dy=Altura, r = 0.36), discretization_length=0.1, current=I_rectangular).Rotate(axis=(1, 0, 0), deg=-90).Translate([0,pos_espira_rectangular[i],0]).Translate([0,0, 0.5])
            sol.AddWire(w2r)
            Longitud.append(np.sqrt(radio_cilindro**2 - pos_espira_rectangular[i]**2)*4+Altura*2)
            continue
            
        if i ==24:
            w2r = wire.Wire(path=wire.Wire.GondolaIDPath(dx=lado_espira_rect*2, dy=Altura, lado=10.46, semieje=1.98, puerta=1.8), discretization_length=0.1, current=-I_rectangular).Rotate(axis=(1, 0, 0), deg=-90).Rotate(axis=(0, 0, 1), deg=180).Translate([0,pos_espira_rectangular[i],0]).Translate([0,0, 0.5])
            sol.AddWire(w2r)
            Longitud.append(np.sqrt(radio_cilindro**2 - pos_espira_rectangular[i]**2)*4+Altura*2)
            continue
        
        if i ==29:
            w2r = wire.Wire(path=wire.Wire.GondolaIDPath(dx=lado_espira_rect*2, dy=Altura, lado=1.91, semieje=0.5, puerta=0.5), discretization_length=0.1, current=-I_rectangular).Rotate(axis=(1, 0, 0), deg=-90).Rotate(axis=(0, 0, 1), deg=180).Translate([0,pos_espira_rectangular[i],0]).Translate([0,0, 0.5])
            sol.AddWire(w2r)
            Longitud.append(np.sqrt(radio_cilindro**2 - pos_espira_rectangular[i]**2)*4+Altura*2)
            continue
        
        else:
        
            w2r = wire.Wire(path=wire.Wire.RectangularPath(dx=lado_espira_rect*2, dy=Altura), discretization_length=0.1, current=I_rectangular).Rotate(axis=(1, 0, 0), deg=-90).Translate([0,pos_espira_rectangular[i],0]).Translate([0,0, 0.5])
            sol.AddWire(w2r)
            Longitud.append(np.sqrt(radio_cilindro**2 - pos_espira_rectangular[i]**2)*4+Altura*2)

    w2r = wire.Wire(path=wire.Wire.RectangularPath(dx=np.sqrt(radio_cilindro**2 - pos_espira_rectangular[0]**2)*2, dy=Altura), discretization_length=0.1, current=I_rectangular).Rotate(axis=(1, 0, 0), deg=-90).Translate([0,pos_espira_rectangular[-1], 0]).Translate([0,0, 0.5])
    sol.AddWire(w2r)
    Longitud.append(np.sqrt(radio_cilindro**2 - pos_espira_rectangular[0]**2)*4+Altura*2)
	
	
	# circular loops I=1 A


    w9c = wire.Wire(path=wire.Wire.PeriodicCable(radius=radio_cilindro, span=3.0, sag=0.26, pts_per_span=50), discretization_length=0.1, current=5*I_circulares).Translate([0,0,36.5])
    sol.AddWire(w9c)    

	
	
    w1c = wire.Wire(path=wire.Wire.PeriodicCable(radius=25, span=3.0, sag=0.26, pts_per_span=50), discretization_length=0.1, current=I_circulares).Translate([0,0,-35.5])
    sol.AddWire(w1c)

    w2c = wire.Wire(path=wire.Wire.PeriodicCable(radius=26, span=3.0, sag=0.26, pts_per_span=50), discretization_length=0.1, current=I_circulares).Translate([0,0, 36.5])
    sol.AddWire(w2c)
    Longitud.append(2*np.pi*27)



    w18c =  wire.Wire(path=wire.Wire.PeriodicCable(radius=radio_cilindro, span=3.0, sag=0.26, pts_per_span=50), discretization_length=0.1, current=5*I_circulares).Translate([0,0,-35.5])
    sol.AddWire(w18c)
    Longitud.append(2*3*np.pi*radio_cilindro)




    for j in range(len(pos_espira_circular)):
        w17c = wire.Wire(path=wire.Wire.PeriodicCable(radius=radio_cilindro, span=3.0, sag=0.26, pts_per_span=50), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[j]])
        sol.AddWire(w17c)
        Longitud.append(2*np.pi*radio_cilindro)

    	



    B1 = sol.CalculateB(points=puntos)*(10**7)

    Bx=[]
    By=[]
    Bz=[]
    B_perp=[]



    for l in range(len(puntos)):
        Bx.append(B1[l,0]+campo[l][0])
        By.append(B1[l,1]+campo[l][1])		
        Bz.append(B1[l,2]+campo[l][2])

        if puntos[l][2] == 32.9 or puntos[l][2] == -32.9:
            B_perp.append(np.sqrt((B1[l,0]+campo[l][0])**2+(B1[l,1]+campo[l][1])**2))
        else:
            B_perp.append(np.sqrt(((B1[l,0]+campo[l][0])*np.sin(Angulo[l])-(B1[l,1]+campo[l][1])*np.cos(Angulo[l]))**2+(B1[l,2]+campo[l][2])**2))
		

		
	# make figure with loops in 3D
    
    PMTs_malos = 0

    for alfa in range(len(B_perp)):
        Coordenadas.append([puntos[alfa][0], puntos[alfa][1], puntos[alfa][2], Bx[alfa], By[alfa], Bz[alfa], B_perp[alfa], 1 if puntos[alfa][2] == 32.9 else 3 if puntos[alfa][2] == -32.9 else 2])
        if np.abs(B_perp[alfa]) > 100:
            PMTs_malos = PMTs_malos+1
            #Coordenadas.append([puntos[alfa]])

			

    return PMTs_malos, B_perp, Coordenadas
	
	

#Resultados

Desviaciones = []

PMTs_malos, B_perp, Coordenadas = Espiras(puntos, campo, Angulo)
Media_total = np.sum(B_perp)/len(B_perp)
for i in range(len(B_perp)):
	Desviaciones.append((B_perp[i]-Media_total)**2)
	
Desviacion_est = np.sqrt(np.sum(Desviaciones)/len(B_perp))
print(Desviacion_est)
print('La media del campo magnetico es de', Media_total)

Porcentaje = PMTs_malos*100/len(puntos)
print('El numero de PMTs malos total es', PMTs_malos)
print('El porcentaje de PMTs malos es', Porcentaje)

df = pd.DataFrame(Coordenadas, columns=['x', 'y', 'z', 'Bx', 'By', 'Bz', 'Bp', 'faceid'])

# Guardar en un archivo Excel
df.to_csv('./Archivo_3.csv', index=False, sep=',')



#Histograma

intervalos = np.arange(0, 210, 5) #calculamos los extremos de los intervalos
plt.hist(x=B_perp, bins=intervalos, color='#F2AB6D', rwidth=0.85)
plt.xlabel("Remaining magnetic field perpendicular to PMT (mG)")
plt.ylabel("Number of PMTs")
plt.ylim(0,6500)
plt.xticks(np.arange(0, 200, 25))
plt.title("$B_{perp}$ distribution for all the PMTs")
textstr = '\n'.join((
    r'$\mu=%.2f\ \mathrm{mG}$' % (Media_total, ),
    r'$\sigma=%.2f\ \mathrm{mG}$' % (Desviacion_est, ),
    r'Prop. excess=%.2f\ \%%' % (Porcentaje, )
))

    
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

plt.text(140,2000, textstr, fontsize=10, bbox=props)
plt.show()
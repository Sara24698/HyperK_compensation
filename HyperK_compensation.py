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



radio_PMT_tapas = 31.95
radio_PMT= 32.4




Media = []
points2=[]
PMTs_top=[]
PMTs_bottom=[]
PMTs_paredes = []
Coordenadas = []
Longitud=[]


#Definición de puntos
x = np.arange(-31.815, 32, 0.707)
y = np.arange(-31.815, 32, 0.707)

Angulo_vuelta = np.arange(0, 6.315, 0.0218)
z=np.arange(-32.522, 32.523, 0.707)
Angulo = np.tile(Angulo_vuelta, len(z))



for r in range(len(x)):
	for h in range(len(y)):
		points2.append([x[r], y[h]])


for g in range(len(points2)):
	distancia = np.sqrt(points2[g][0]**2+points2[g][1]**2)
	if distancia <=32.01:
		PMTs_top.append([points2[g][0], points2[g][1], 32.9])
		PMTs_bottom.append([points2[g][0], points2[g][1], -32.9])


for i in range(len(z)):
	for j in range(len(Angulo_vuelta)):
		PMTs_paredes.append([radio_PMT*np.cos(Angulo_vuelta[j]), radio_PMT*np.sin(Angulo_vuelta[j]), z[i]])
          

#Función principal
def Sistema_compensacion(puntos, visualize = False, elipticas = False):

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

    	

	
    if elipticas == True:
        #Parte de abajo
        w1e = wire.Wire(path=wire.Wire.EllipticalPath(rx=9.05, ry=19.05, pts=20), discretization_length=0.1, current=-95).Rotate(axis=(0,0,1), deg=90).Translate([0,0,-35.5]).Translate([0,-24.7,0])
        sol.AddWire(w1e)
        Longitud.append(2*np.pi*np.sqrt((9.05**2+19.05**2)/2))

        
        w2e = wire.Wire(path=wire.Wire.EllipticalPath(rx=14.765, ry=24.77, pts=20), discretization_length=0.1, current=-190).Rotate(axis=(0,0,1), deg=90).Translate([0,0,36.5]).Translate([0,18.7,0])
        sol.AddWire(w2e)
        Longitud.append(2*np.pi*np.sqrt((14.765**2+24.77**2)/2))

        w3e = wire.Wire(path=wire.Wire.EllipticalPath(rx=17.8, ry=27.8, pts=20), discretization_length=0.1, current=-120).Rotate(axis=(0,0,1), deg=90).Translate([0,0,36.5]).Translate([0, 15,0])#I=-148
        sol.AddWire(w3e)
        Longitud.append(2*np.pi*np.sqrt((20.48**2+30.49**2)/2))

        w4e = wire.Wire(path=wire.Wire.EllipticalPath(rx=9.05, ry=19.05, pts=20), discretization_length=0.1, current=50).Rotate(axis=(0,0,1), deg=90).Translate([0,0,-35.5]).Translate([0,24.7,0])
        sol.AddWire(w4e)
        Longitud.append(2*np.pi*np.sqrt((9.05**2+19.05**2)/2))
        

        #Parte de arriba

        w5e = wire.Wire(path=wire.Wire.EllipticalPath(rx=17.8, ry=27.8, pts=20), discretization_length=0.1, current=-160).Rotate(axis=(0,0,1), deg=90).Translate([0,0,36.5]).Translate([0, 15,0])#I=150
        sol.AddWire(w5e)
        Longitud.append(2*np.pi*np.sqrt((20.48**2+30.49**2)/2))

        w6e = wire.Wire(path=wire.Wire.EllipticalPath(rx=14.765, ry=24.77, pts=20), discretization_length=0.1, current=-190).Rotate(axis=(0,0,1), deg=90).Translate([0,0,36.5]).Translate([0,18.7,0])
        sol.AddWire(w6e)
        Longitud.append(2*np.pi*np.sqrt((14.765**2+24.77**2)/2))

        w7e = wire.Wire(path=wire.Wire.EllipticalPath(rx=9.05, ry=19.05, pts=20), discretization_length=0.1, current=-25).Rotate(axis=(0,0,1), deg=90).Translate([0,0,36.5]).Translate([0,24.7,0])
        sol.AddWire(w7e)
        Longitud.append(2*np.pi*np.sqrt((9.05**2+19.05**2)/2))
	
    
    
    if visualize == True:
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection='3d')
        sol.mpl3d_PlotWires(ax)
        plt.show()


    B = sol.CalculateB(points=puntos)*(10**7)

    Bx=[]
    By=[]
    Bz=[]
    B_perp=[]



    if len(puntos) == len(PMTs_bottom):
        for q in range(len(puntos)):
            Bx.append(B[q,0])
            By.append(B[q,1]+303)		
            Bz.append(B[q,2]-366)
            B_perp.append(np.sqrt(B[q,0]**2+(B[q,1]+303)**2))
            #Bx = [B[q,0] for q in range (len(puntos))]
            #By = [B[q,1]+303 for q in range (len(puntos))]
            #Bz = [B[q,2]-366 for q in range (len(puntos))]
            

		
    else:
        for l in range(len(puntos)):
            Bx.append(B[l,0])
            By.append(B[l,1]+303)		
            Bz.append(B[l,2]-366)
            B_perp.append(np.sqrt((B[l,0]*np.sin(Angulo[l])-(B[l,1]+303)*np.cos(Angulo[l]))**2+(B[l,2]-366)**2))
            
		

		
	# make figure with loops in 3D
    
    PMTs_malos = 0

    for alfa in range(len(B_perp)):
        Coordenadas.append([puntos[alfa][0], puntos[alfa][1], puntos[alfa][2], Bx[alfa], By[alfa], Bz[alfa], B_perp[alfa], 1 if puntos[alfa][2] == 32.9 else 3 if puntos[alfa][2] == -32.9 else 2])
        if np.abs(B_perp[alfa]) > 100:
            PMTs_malos = PMTs_malos+1

			

    return PMTs_malos, B_perp, Coordenadas
	
	


#Extracción de parámetros

PMTs_malos_top, B_perp_top, Coordenadas_top = Sistema_compensacion(PMTs_top, visualize=True, elipticas=False)
PMTs_malos_bottom, B_perp_bottom, Coordenadas_bottom= Sistema_compensacion(PMTs_bottom, visualize=False, elipticas=False)
PMTs_malos_paredes, B_perp_paredes, Coordenadas_paredes = Sistema_compensacion(PMTs_paredes, visualize=False, elipticas=False)	
	
	



#Resultados

Desviaciones = []
B_perp_total = B_perp_bottom+B_perp_top+B_perp_paredes
Media_total = np.sum(B_perp_total)/len(B_perp_total)
for i in range(len(B_perp_total)):
	Desviaciones.append((B_perp_total[i]-Media_total)**2)
	
Desviacion_est = np.sqrt(np.sum(Desviaciones)/len(B_perp_total))
PMTs_malos_total = PMTs_malos_paredes+PMTs_bottom+PMTs_malos_top


print(f"La media total es {Media_total:.2f} ± {Desviacion_est:.2f}")
print(f"La media del campo magnético es de {Media_total}")
print(f"El número de PMTs malos en top es {PMTs_malos_top}")
print(f"El número de PMTs malos en bottom es {PMTs_malos_bottom}")
print(f"El número de PMTs malos en las paredes es {PMTs_malos_paredes}")

print(f"El número total de PMTs malos es {PMTs_malos_total}")
print(f"El número de PMTs en la pared es {len(z) * len(Angulo)}, en cada una de las tapas {len(PMTs_top)}, y en total en el detector hay {len(B_perp_total)}")
print(f"El porcentaje de PMTs malos es {PMTs_malos_total * 100 / len(B_perp_total):.2f}%")


#df = pd.DataFrame(Coordenadas, columns=['x', 'y', 'z', 'Bx', 'By', 'Bz', 'Bp', 'faceid'])

# Guardar en un archivo Excel
#df.to_excel('./Gondolas.xlsx', index=False)



#Histograma

intervalos = np.arange(0, 210, 5) #calculamos los extremos de los intervalos
plt.hist(x=B_perp_total, bins=intervalos, color="#080049", rwidth=0.85)
plt.xlabel("Remaining magnetic field perpendicular to PMT (mG)")
plt.ylabel("Number of PMTs")
plt.ylim(0,6500)
plt.xticks(np.arange(0, 200, 25))
plt.title("$B_{perp}$ distribution for all the PMTs")
textstr = '\n'.join((
    r'$\mu=%.2f\ \mathrm{mG}$' % (Media_total, ),
    r'$\sigma=%.2f\ \mathrm{mG}$' % (Desviacion_est, ),
    r'Prop. excess=%.2f%' % (PMTs_malos_total * 100 / len(B_perp_total), )
))

    
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

plt.text(140,2000, textstr, fontsize=10, bbox=props)
plt.show()

"""

#Puntos malos en el cilindro

Coord_x=[]
Coord_y=[]
Coord_z=[]

x_top=[]
y_top=[]
z_top=[]

x_bottom=[]
y_bottom=[]
z_bottom=[]




for indice in range(len(Coordenadas_paredes)):
	Coord_x.append(Coordenadas_paredes[indice][0])	 
	Coord_y.append(Coordenadas_paredes[indice][1])	 
	Coord_z.append(Coordenadas_paredes[indice][2])


for indice_top in range(len(Tapa_superior[3])):
	x_top.append(Tapa_superior[3][indice_top][0])	 
	y_top.append(Tapa_superior[3][indice_top][1])	 
	z_top.append(Tapa_superior[3][indice_top][2])


for indice_bottom in range(len(Tapa_inferior[3])):
	x_bottom.append(Tapa_inferior[3][indice_bottom][0])	 
	y_bottom.append(Tapa_inferior[3][indice_bottom][1])	 
	z_bottom.append(Tapa_inferior[3][indice_bottom][2])	


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')


ax.scatter(x_bottom, y_bottom, z_bottom, label ='PMTs bottom')
plt.plot()
"""

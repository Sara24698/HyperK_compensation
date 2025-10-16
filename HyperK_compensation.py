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
import sqlite3




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


def posiciones(geomagnetico=True, rebar=False):

    if geomagnetico==True:

        campo = [0, 303, -366]
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

        return PMTs_top, PMTs_bottom, PMTs_paredes, Angulo, campo
    
    if rebar==True:
        df= pd.read_csv('./efecto_pos3.csv', header=None, sep=';',names=['x', 'y', 'z', 'Bx', 'By', 'Bz', 'Btotal', 'Bx_T', 'By_T', 'Bz_T'])
        df = df.drop(['Btotal','Bx_T', 'By_T', 'Bz_T'], axis=1)
        puntos = df[['x', 'y', 'z']].values.tolist()
        campo = df[['Bx', 'By', 'Bz']].values.tolist()
        Angulo = np.arctan2(df['y'], df['x'])
        

        return puntos, campo, Angulo

def rotacion_campo(Angulos_rotacion, campo, geomagnetico = True, rebar=False):
    Angulos_rad = np.deg2rad(Angulos_rotacion)

    if geomagnetico==True:
        Bx_desviado = -campo[1] * np.sin(Angulos_rad)
        By_desviado = campo[1] * np.cos(Angulos_rad)
        Bz_desviado = campo[2]

        return Bx_desviado, By_desviado, Bz_desviado
    
    if rebar == True:
        campo = np.array(campo)  # forma (N,3)
        Bx_desviado = []
        By_desviado = []
        Bz_desviado = []

        for theta in Angulos_rad:
            Rz = np.array([
                [np.cos(theta), -np.sin(theta), 0],
                [np.sin(theta),  np.cos(theta), 0],
                [0, 0, 1]
            ])
            campo_rotado = campo @ Rz.T  # (39800,3)
            Bx_desviado.append(campo_rotado[:,0][:, np.newaxis])  # shape (39800,1)
            By_desviado.append(campo_rotado[:,1][:, np.newaxis])
            Bz_desviado.append(campo_rotado[:,2][:, np.newaxis])

        # Apilamos horizontalmente
        Bx_desviado = np.hstack(Bx_desviado)  # (39800, n_angles)
        By_desviado = np.hstack(By_desviado)
        Bz_desviado = np.hstack(Bz_desviado)
        print(np.shape(Bx_desviado))
        
        return Bx_desviado, By_desviado, Bz_desviado


#Función principal

def Sistema_compensacion(puntos, Angulo, campo, Angulos_rotacion, visualize = False, elipticas = False, geomagnetico=False, rebar=True):

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

    Bx_desviado, By_desviado, Bz_desviado = rotacion_campo(Angulos_rotacion, campo, geomagnetico, rebar)
    B = sol.CalculateB(points=puntos)*(10**7)

    

    parametros = []
    for i in range(len(Angulos_rotacion)):
        Bx = B[:,0] + Bx_desviado[:, i]
        By = B[:,1] + By_desviado[:, i]
        Bz = B[:,2] + Bz_desviado[:, i]
        B_perp=[]

        for q in range(len(puntos)):
            if len(puntos) == len(PMTs_bottom):
                B_perp.append(np.sqrt(Bx[q]**2 + By[q]**2))
            else:
                B_perp.append(np.sqrt(
                    (Bx[q]*np.sin(Angulo[q]) - By[q]*np.cos(Angulo[q]))**2 +
                    Bz[q]**2
                ))
                
                      
        # make figure with loops in 3D
        
        PMTs_malos = 0
        Coordenadas_i = []

        for alfa in range(len(B_perp)):
            Coordenadas_i.append([
                puntos[alfa][0], puntos[alfa][1], puntos[alfa][2],
                Bx[alfa], By[alfa], Bz[alfa], B_perp[alfa],
                1 if puntos[alfa][2] == 32.9 else 3 if puntos[alfa][2] == -32.9 else 2
            ])
            if np.abs(B_perp[alfa]) > 100:
                PMTs_malos += 1

        parametros.append({
            "angulo": Angulos_rotacion[i],
            "PMTs_malos": PMTs_malos,
            "B_perp": B_perp,
            "Coordenadas": Coordenadas_i
        })

    return parametros

def export_efficiency_data(parametros, nombre_base='Resultados'):
    """
    Exporta las coordenadas y campos B_perp de todos los ángulos a archivos CSV separados.
    Cada archivo se llama: <nombre_base>_angulo_<valor>.csv
    """
    conn = sqlite3.connect('./Compensated_field.db')
    cursor = conn.cursor()
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS Coordenadas (
        angulo REAL,
        x REAL,
        y REAL,
        z REAL,
        Bx REAL,
        By REAL,
        Bz REAL,
        Bp REAL,
        faceid INTEGER
    )
    """)
    conn.commit()


    for i in range(len(parametros)):
        angulo = parametros[i]["angulo"]

        # Extraer coordenadas de ese ángulo
        Coordenadas = np.array(parametros[i]["Coordenadas"])

        # Crear DataFrame
        df = pd.DataFrame(Coordenadas, columns=['x', 'y', 'z', 'Bx', 'By', 'Bz', 'Bp', 'faceid'])

        # Guardar CSV (uno por ángulo)
        nombre_archivo = f'./{nombre_base}_angulo_{angulo:.1f}.csv'
        df.to_csv(nombre_archivo, index=False, sep=',')

        df['angulo'] = angulo

        # Insertar en SQL
        df.to_sql('Coordenadas', conn, if_exists='append', index=False)

    conn.close()




def export_resultados(Resultados, mode, nombre_archivo='Resultados'):
    """
    Exporta un único archivo Excel con los resultados por ángulo.
    Cada fila corresponde a un ángulo de rotación.
    """
    if mode =='geomagnetico':
        df = pd.DataFrame(Resultados, columns=[
            'Angle (º)',
            'Mean (mG)',
            'Standard Deviation (mG)',
            'Exc. top',
            'Exc. bottom ',
            'Exc. paredes',
            'Prop. Excess (%)'
        ])
    
    if mode == 'rebar':
        df = pd.DataFrame(Resultados, columns=[
            'Angle (º)',
            'Mean (mG)',
            'Standard Deviation (mG)',
            'Total PMTs with excess',
            'Prop. Excess (%)'
        ])

    df.to_excel('./'+nombre_archivo+'.xlsx', index=False)


def hist(B_perp_total, Media_total, Desviacion_est, PMTs_malos_total, nombre_base, angulo):
    """
    Guarda un histograma de B_perp_total para un ángulo específico.
    """
    intervalos = np.arange(0, 210, 5)  # intervalos de la histograma
    plt.figure(figsize=(8,5))
    plt.hist(x=B_perp_total, bins=intervalos, color="#080049", rwidth=0.85)
    plt.xlabel("Remaining magnetic field perpendicular to PMT (mG)")
    plt.ylabel("Number of PMTs")
    plt.ylim(0, 6500)
    plt.xticks(np.arange(0, 200, 25))
    plt.title(f"$B_{{perp}}$ distribution for all PMTs (Angle = {angulo}°)")
    
    # Cuadro de estadísticas
    textstr = '\n'.join((
        r'$\mu=%.2f\ \mathrm{mG}$' % (Media_total, ),
        r'$\sigma=%.2f\ \mathrm{mG}$' % (Desviacion_est, ),
        r'Prop. excess=%.2f\ \%%' % (PMTs_malos_total * 100 / len(B_perp_total), )
    ))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    plt.text(140, 2000, textstr, fontsize=10, bbox=props)
    
    # Guardar archivo con nombre que incluye el ángulo
    nombre_archivo = f'./{nombre_base}_angulo_{angulo:.1f}.png'
    plt.savefig(nombre_archivo)
    plt.close()

	

#Resultados
def resultados(Angulos_rotacion, mode, export_ef_data=True, histogram=True, export_results=True):

    if mode == 'geomagnetico':
        PMTs_top, PMTs_bottom, PMTs_paredes, Angulo, campo = posiciones(geomagnetico=True, rebar=False)
        parametros_top = Sistema_compensacion(PMTs_top, Angulo, campo, Angulos_rotacion, visualize=False, elipticas=False, geomagnetico=True, rebar=False)
        parametros_bottom= Sistema_compensacion(PMTs_bottom, Angulo, campo, Angulos_rotacion, visualize=False, elipticas=False, geomagnetico=True, rebar=False)
        parametros_paredes = Sistema_compensacion(PMTs_paredes, Angulo, campo, Angulos_rotacion, visualize=False, elipticas=False, geomagnetico=True, rebar=False)	
        
        parametros = parametros_top + parametros_bottom + parametros_paredes


        Resultados = []  # Aquí guardaremos una fila por cada ángulo
        
        # Recorremos todos los ángulos
        for i, angulo in enumerate(Angulos_rotacion):
            # Extraemos los resultados correspondientes a ese ángulo
            B_perp_top = parametros_top[i]["B_perp"]
            B_perp_bottom = parametros_bottom[i]["B_perp"]
            B_perp_paredes = parametros_paredes[i]["B_perp"]

            PMTs_malos_top = parametros_top[i]["PMTs_malos"]
            PMTs_malos_bottom = parametros_bottom[i]["PMTs_malos"]
            PMTs_malos_paredes = parametros_paredes[i]["PMTs_malos"]

            # Combinamos todo
            B_perp_total = np.concatenate((B_perp_top, B_perp_bottom, B_perp_paredes))

            # Calculamos las estadísticas
            Media_total = np.mean(B_perp_total)
            Desviacion_est = np.std(B_perp_total)
            PMTs_malos_total = PMTs_malos_top + PMTs_malos_bottom + PMTs_malos_paredes
            Prop_malos = PMTs_malos_total * 100 / len(B_perp_total)

            # Añadimos los resultados de este ángulo
            Resultados.append([
                angulo,
                Media_total,
                Desviacion_est,
                PMTs_malos_top,
                PMTs_malos_bottom,
                PMTs_malos_paredes,
                Prop_malos
            ])



            print(f"La media total es {Media_total:.2f} ± {Desviacion_est:.2f}")
            print(f"El número de PMTs malos en top es {PMTs_malos_top}")
            print(f"El número de PMTs malos en bottom es {PMTs_malos_bottom}")
            print(f"El número de PMTs malos en las paredes es {PMTs_malos_paredes}")

            print(f"El número total de PMTs malos es {PMTs_malos_total}")
            print(f"El porcentaje de PMTs malos es {PMTs_malos_total * 100 / len(B_perp_total):.2f}%")

            if histogram==True:
                hist(B_perp_total, Media_total, Desviacion_est, PMTs_malos_total, nombre_base='Histograma', angulo=angulo)

        if export_ef_data:
            export_efficiency_data(parametros, nombre_base='Datos')


        if export_results:
            export_resultados(Resultados, nombre_archivo='Compensacion')
    
    if mode == 'rebar':
        puntos, campo, Angulo = posiciones(geomagnetico=False, rebar=True)
        parametros= Sistema_compensacion(puntos, Angulo, campo, Angulos_rotacion, visualize=False, elipticas=False, geomagnetico=False, rebar=True)
        Resultados = []  # Aquí guardaremos una fila por cada ángulo
        
        # Recorremos todos los ángulos
        for i, angulo in enumerate(Angulos_rotacion):
            # Extraemos los resultados correspondientes a ese ángulo
            B_perp_total= parametros[i]["B_perp"]
            PMTs_malos_total = parametros[i]["PMTs_malos"]


            # Calculamos las estadísticas
            Media_total = np.mean(B_perp_total)
            Desviacion_est = np.std(B_perp_total)
            Prop_malos = PMTs_malos_total * 100 / len(B_perp_total)

            # Añadimos los resultados de este ángulo
            Resultados.append([
                angulo,
                Media_total,
                Desviacion_est,
                PMTs_malos_total,
                Prop_malos
            ])



            print(f"La media total es {Media_total:.2f} ± {Desviacion_est:.2f}")
            print(f"El número total de PMTs malos es {PMTs_malos_total}")
            print(f"El porcentaje de PMTs malos es {PMTs_malos_total * 100 / len(B_perp_total):.2f}%")

            if histogram==True:
                hist(B_perp_total, Media_total, Desviacion_est, PMTs_malos_total, nombre_base='Histograma', angulo=angulo)

        if export_ef_data:
            export_efficiency_data(parametros, nombre_base='Datos')


        if export_results:
            export_resultados(Resultados, mode, nombre_archivo='Compensacion')


    
resultados([0,3], mode='rebar', export_ef_data=True, histogram=True, export_results=True)


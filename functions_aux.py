#FUNCIÓ PER LLEGIR ELS FITXERS
def lector_fitxers(ruta):
    """
    Donat un determinat fitxer (la seva ruta) corresponent a una certa variable d'una certa execució del codi GR-EM-KG el llegeix i converteix en un np array de dimensió 3.

    input:
        ruta: string, ruta al fitxer desitjat
    output:
        dades: array de dimensió 3 [temps (t), variable (y), posició (x)]
    """
    import numpy as np

    #declaram dues llistes, una per es blocs totals (els anirem acumulant) i una per es bloc actual
    blocs_total = []
    bloc = []
    t = 0.0
    
    with open(ruta, 'r') as file:
        for line in file:
            #netejam cada linia
            line = line.strip()

            # botam línies buides
            if not line:
                continue  
                
            # detectam les linies que ens indiquen el temps
            if "Time =" in line:
                #si ja teniem dades de un temps anterior, les guardam
                if len(bloc) > 0:
                    blocs_total.append(bloc)
                    bloc = []
                
                # extreim el valor numèric del temps
                line_clean = line.replace('"', '')
                t = float(line_clean.split('=')[1].strip())

            # detectam les linies de dades x | y         
            else:
                columnes = line.split()
                if len(columnes) == 2:
                    col1 = float(columnes[0])
                    col2 = float(columnes[1])
                    # composam totes les dades concretes: [t, col1, col2]
                    bloc.append([t, col1, col2])
                    
        # quan acabam de llegir l'arxiu, s'ha de guardar el darrer bloc
        if len(bloc) > 0:
            blocs_total.append(bloc)
            
    # Convertimos la lista de listas en un array 3D puro de NumPy
    return np.array(blocs_total)

#FUNCIÓ PER EXTREURE SA INFO DE CONVERGÈNCIA
def info_converg(variable, ruta, constraint=0, scaler=3):
    """
    Donada una certa variable per una simulació que es troba en una certa ruta, aquesta funció s'encarrega de retornar les variables de x(posicions), t(temps), i les respectives a les mètriques de comparació que li pertoquin segons el valor de constraint 
    
    input:
        variable : nom complet de la variable, exemple: evkg_repi_c4
        ruta : ruta a on es troba la simulacio
        constraint : selector de forma de la comparació, =0 per evkg, =1 per constraint...
        scaler: factor de escala utilitzat al test de convergència, per defecte és 3, però pot ser 1.5
    """
    import numpy as np

    
    dades_conv1 = lector_fitxers(ruta + "_" + variable + "_convres1.dat")
    dades_conv2 = lector_fitxers(ruta + "_" + variable + "_convres2.dat")
    dades_conv3 = lector_fitxers(ruta + "_" + variable + "_convres3.dat")

    x = dades_conv1[0,:,1] #extreim la línia de posició, que serà igual per a tots els temps i les dues variables
    times = dades_conv1[:,0,0] #extreim els temps, que seran igual per les dues variables

    if scaler==1.5:
            dades_conv4 = lector_fitxers(ruta + "_" + variable + "_convres4.dat")
            x_1 = dades_conv1[0,:,1] #en aquest cas cada un dels eixos x serà diferent per a cada convres
            x_2 = dades_conv2[0,:,1]
            x_3 = dades_conv3[0,:,1]
            x_4 = dades_conv4[0,:,1]

            x = [x_1, x_2, x_3, x_4]

            times_1 = dades_conv1[:,0,0] #idem per als t's
            times_2 = dades_conv2[:,0,0]
            times_3 = dades_conv3[:,0,0]
            times_4 = dades_conv4[:,0,0]

            times = times_1

            #el fet que després representem els x1 amb results1, x2 amb results 2 no es problema que tenguin diferents shapes perque domes ens importen les constraints, que no es mesclen les convres entre sí

            #ara bé, pels temps, necessitam homogeneitzar els eixos de les quatre convres

            results_1 = dades_conv1[:,:,2]/1.5**12 
            results_2 = dades_conv2[:,:,2]/1.5**8
            results_3 = dades_conv3[:,:,2]/1.5**4
            results_4 = dades_conv4[:,:,2]
            
            return x, times, results_1, results_2, results_3, results_4


    if constraint==0:
    #cas de les variables repi, rephi, impi, imphi
        results_1 =(dades_conv1[:,:,2] - dades_conv2[:,:,2])/3**4
        results_2 = dades_conv2[:,:,2] - dades_conv3[:,:,2]
        
        return x, times, results_1, results_2

    if constraint==1:
    #cas de les variables ham, grr, etc.
     #això es lo que mos ha dit n'alex, crec que es incorrecte
        results_1 = dades_conv1[:,:,2]/3**8 
        results_2 = dades_conv2[:,:,2]/3**4
        results_3 = dades_conv3[:,:,2]

    
        return x, times, results_1, results_2, results_3

#FUNCIÓ PER PINTAR EL FOTOGRAMA, USADA AMB INTERACT DE IPYWIDGETS
def pintar_fotograma(variable, constraint, x, t,
                    results_1, results_2, results_3=None, results_4=None,
                    i=0, scaler=3):
    """ Funció pintar_fotograma s'encarrega de pintar el fotograma corresponent a cada pas temporal per a després emprar-lo amb ipywidgets interact"""
    import matplotlib.pyplot as plt
            
    plt.figure(figsize=(8, 4))

    if scaler==1.5:
        x_1 = x[0]
        x_2 = x[1]
        x_3 = x[2]
        x_4 = x[3]

        plt.plot(x_1, results_1[i, :], "-", label="first part")
        plt.plot(x_2, results_2[i, :], "--", label="second part")
        plt.plot(x_3, results_3[i, :], "-", label="third part")
        plt.plot(x_4, results_4[i, :], "--", label="forth part")
        plt.ylim(-0.02,0.02)
    if scaler==3:
        if constraint==0: #evolution variables
            plt.plot(x, results_1[i, :], "-", label="first part")
            plt.plot(x, results_2[i, :], "--", label="second part")
            plt.ylim(-0.00002, 0.00002)
        if constraint==1: #constraints
            plt.plot(x, results_1[i, :], "-", label="first part")
            plt.plot(x, results_2[i, :], "--", label="second part")
            plt.plot(x, results_3[i, :], "-", label="third part")
            plt.ylim(-0.001,0.001)


    #plt.ylim(-0.1,1)
    plt.xlim(0,1)
    plt.title("Variable " + variable + "    "f"Tiempo: {t[i]}")
    plt.legend()
    plt.show()







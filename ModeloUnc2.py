### MODELO TERMODINAMICO DETERMINISTA DE IMPACTO DE UNCINARIAS ###
### HOOKWORM BURDEN ON A SINGLE HOST: A THERMODYNAMIC MODEL ###

import random;
import math;

# Pedirle a usuario valores de parametros

def init():
	input = raw_input("Inserta tus parametros en formato \"Ein,Eout,Eunc,DporCiclo,Lon,Mamb,Pinf,Fert,Tdes\": "); # guarda input
	global param;
	param = map(float, input.split(",")); # guarda parametros
	if len(param) != 9:
		print "Ups. Numero de parametros equivocado.";
		init();

init();

# PARAMETROS

Ein = param[0];				# Energia que entra diariamente al sistema @param (2000kcal, adjust to http://en.wikipedia.org/wiki/File:World_map_of_Energy_consumption_2001-2003.svg)
Eout = param[1];			# Energia gastada por metabolismo (5.44-7.53 MJ, Life: The Science of Biology. Ed. William K. Purves and David Sadava. 7th ed. New York: Freeman, 2004: 962. http://hypertextbook.com/facts/2009/VickieWu.shtml)
Eunc = param[2]; 			# Energia perdida por uncinaria por dia @param ( -1 kcal/(dia*uncinaria), Briscoe, J. (1979). The quantitative effect of infection on the use of food by young children in poor countries. The American journal of clinical nutrition,32(3), 648-676.)
DporCiclo = int(param[3]);	# Numero de dias en un ciclo del programa
Lon = param[4];				# Longevidad promedio en dias de una uncinaria @param
Mamb = param[5];			# Probabilidad de que un huevo o uncinaria en ambiente muera en un dia dado @param
Pinf = param[6];			# Probabilidad de infeccion por una uncinaria en el ambiente @param
Fert = param[7];			# Cantidad de huevos que pone uncinaria por dia
Tdes = param[8];			# Tiempo en dias de desarrollo promedio de un huevo de uncinaria hasta llegar al intestino @param http://www.phac-aspc.gc.ca/lab-bio/res/psds-ftss/ancylostoma-duodenale-eng.php

# VARIABLES
Enet = (Ein-Eout)*DporCiclo;	# Energia total en el sistema en el ciclo
Unc = 1; 						# Numero de uncinarias en hospedero 
EunNet = 0; 					# Energia perdida por uncinarias en un dia
Mue = 0;						# Numero de uncinarias que mueren en un dia
Inf = 1;						# Numero de uncinarias nuevas que infectan al hospedero en un dia
Mimu = 0;						# Probabilidad de que el sistema inmune mate una uncinaria
Pmna = 0;						# Probabilidad de que una uncinaria dentro de un humano muera en un dia dado.  
UncAmb = 0;						# Numero de uncinarias en el ambiente
Hue = 0;						# Numero de huevos puestos en un dia
mod = 0;						# Modifica el parametro de Pinf para crear un negative feedback loop
Tbirth = int(round(Tdes*1.0/DporCiclo));		# Tiempo en ciclos de incubacion de huevos y desarrollo a L3 infectivas 

# CONTAINERS
ciclos = [ [Enet, Unc, EunNet, Mue, Inf, Mimu, UncAmb, Hue, 2] ];
				# Contiene todos los ciclos. ciclos[0] es caso base, ciclos[1] es primera iteracion. Ultima posicion del array contiene numero real del anyo. 

# STOCHASTIC FCTN (random number check)
def prob(p):
	r = random.random();
	resp = False;
	if r <= p:
		resp = True;
	return resp;
	
print "ciclo", 0, "Ein-Eout", (Ein-Eout)*DporCiclo, "Enet", Enet, "mod", mod;
# MAIN LOOP
for i in range(1,6400): # itera un periodo de desarrollo. Empieza con 1 para igualarse a indice de ciclos.
	
	Hue = int(round(Unc*0.5*random.gauss(Fert*DporCiclo, 2500/math.sqrt(DporCiclo)))); # Asumimos que la mitad de las uncinarias son hembras y todas ponen huevos. @param
	
	#if i>Tbirth:
	#	UncAmb = ciclos[i-Tbirth][7]; # nacen las uncinarias de huevos puestos hace un ciclo.
	UncAmb = 0;
	max = 14/DporCiclo;
	for k in range(1,max+1):
		if i>Tbirth+k:
			UncAmb += ciclos[i-Tbirth-k][7]; # siguen vivas las uncinarias de huevos puestos hace dos semanas.
	
	UncAmb = int(round(UncAmb*(1-Mamb))); # ambiente mata uncinarias en ambiente
	
	Inf = 0;
	p = 0.000001;
	if i>int(round(9125/DporCiclo)):
		p = Pinf*0.000001;#/round(1+i*DporCiclo/(365.0)); # cada anyo disminuye probabilidad de infeccion
	for j in range(0, UncAmb):
		if prob(p): 
			Inf += 1; # sets infecciones nuevas 
	
	#Mimu = (0.25*Enet/2092000 + 0.25); # sets probabilidad de muerte por sistema inmune @param
	Mue = 0;
	if Enet>0:
		mod = ((1.0*(Ein-Eout)*DporCiclo/Enet)**(120))/10;
		Pmna = (1/Lon)*mod*DporCiclo; # Negative feedback loop
	else:
		Pmna = 1;
	for j in range(0, Unc):
		if prob(Pmna):
			Mue += 1; # sets muertes 
	#Mue = (Unc*(Pmna+Mimu)); # sets muertes 
		
	Unc = Unc + Inf - Mue; # sets numero de uncinarias en hospedero
	
	EunNet = Unc*Eunc*DporCiclo; # sets energia total perdida por uncinarias
	Enet = Ein*DporCiclo - Eout*DporCiclo - EunNet; # sets energia total ese ciclo
	
	ciclos.append( [Enet, Unc, EunNet, Mue, Inf, Mimu, UncAmb, Hue, 2+i*DporCiclo/365.0] ); # guarda datos de este ciclo (al final, tiempo en anyos)
	print "ciclo", i, "Ein-Eout", (Ein-Eout)*DporCiclo, "Enet", Enet, "mod", mod, "Pmna", Pmna;

# Exportar archivo datos
data = open("data.csv", "w"); # crea archivo para guardar datos en csv. SI EXISTE ARCHIVO, SOBREESCRIBE.
data.write("Ciclo,Enet,Unc,EunNet,Mue,Inf,Mimu,UncAmb,Hue,Anyo\n");
for i in range(0,len(ciclos)):
	data.write(str(i) + "," + ",".join(map(str, ciclos[i])) + "\n"); # pasa elementos de ciclos de int a str, concatena dia como string, escribe cada dia en una linea
data.close(); # cierra archivo

print "done.";
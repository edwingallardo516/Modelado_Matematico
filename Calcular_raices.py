import numpy as np
import sympy as sym # para matematicas simbolicas 
from sympy.plotting import plot # Libreria para graficar
from matplotlib import pyplot as plt
import math

# Desarrollar el polinomio
# Al calcular r se forma el siguiente polinomio
# (pi)^2 r^4 + (pi)^2 h^2 r^2 - s^2 = 0

# Parametros de la función
r = sym.Symbol('r')
s = 750
h = 9
fx = (math.pi)**2 *r**4 + (math.pi)**2 *h**2 *r**2 - (s**2) # Función
fx_prima = fx.diff(r) # Derivada de la función
fx_prima2 = fx_prima.diff(r) # Derivada segunda de la función


# 1. Graficamos la función para ver donde se encuentra la raíz 

p_fx = plot(fx,(r,13,16),title='Raíz de la función',show=True)

# 2. Metodos para calcular las raices
ite = 10 # Número de iteraciones a utilizar 

# 2.1. Método de Bisección ------------------------------------------------------------------------------
# Utilizando el procedimiento de la figura 5.5 del libro Metodos numéricos para ingenieros 5a ed. pag. 124

# Paso 1:
# Según la grafica el rango de valores donde está la raíz es 13-16, por tanto:
xi = 13
xu = 16
#Comprobación:
f_xi_xu = fx.subs(r,xi)*fx.subs(r,xu)
print("  1. MÉTODO DE BISECCIÓN  ")

# Paso 2:
xr_anterior = 16 # Guarda el valor de la raíz aproximado anterior para hacer el calculo de errores
error1 = np.zeros(ite) # Vector para guardar los errores
for i in range(0,ite):
    xr = (xi + xu)/2

    f_xi_xr = fx.subs(r,xi)*fx.subs(r,xr) 
    if(f_xi_xr < 0):
        xu = xr
    if(f_xi_xr > 0):
        xi = xr
    if(f_xi_xr == 0):
        i = ite # termina el proceso
    
    # Calculo de errores
    # Utilizando la ecuación 5.2 del libro Metodos numéricos para ingenieros 5a ed. pag. 126
    Erel = abs((xr - xr_anterior)/xr)*100
    xr_anterior = xr 
    error1[i] = Erel
    print("niter= ", i+1, "xi= ",round(xr,6), "f(xi)= ",round(fx.subs(r,xr),6), "Erel(%)= ", round(Erel,6)) # Presentación de resultados

# 2.2. Método Regula Falsi -------------------------------------------------------------------------------------
 # Utilizando la ecuación 5.7 del libro Metodos numéricos para ingenieros 5a ed. pag. 132
 # Se utiliza el mismo intervalo que en el caso anterior 
xi = 13
xu = 16

print("\n"); print("  2. MÉTODO REGULA FALSI  ")
xr_anterior = 16 # Guarda el valor de la raíz aproximado anterior para hacer el calculo de errores
error2 = np.zeros(ite) # Vector para guardar los errores
for i in range(0,ite):
    xr = xu - (fx.subs(r,xu)*(xi-xu))/(fx.subs(r,xi) - fx.subs(r,xu))

    f_xi_xr = fx.subs(r,xi)*fx.subs(r,xr)
    if(f_xi_xr < 0):
        xu = xr
    if(f_xi_xr > 0):
        xi = xr
    if(f_xi_xr == 0):
        i = ite # termina el proceso

    # Calculo de errores
    # Utilizando la ecuación C6.2.6 del libro Metodos numéricos para ingenieros 5a ed. pag. 150
    Erel = abs((xr - xr_anterior)/xr)*100
    xr_anterior = xr
    error2[i] = Erel
    print("niter= ", i+1, "xi= ",round(xr,6), "f(xi)= ",round(fx.subs(r,xr),6), "Erel(%)= ", round(Erel,6)) # Presentación de resultados

# 2.3. Método Newton - Raphson --------------------------------------------------------------------------------
# utilizando la ecuación 6.6 del libro Metodos numéricos para ingenieros 5a ed. pag. 149

# Según la grafica, el valor aproximado de la raíz es 14, pero por probar la eficacia del método
# lo empezaremos en 20
xi = 20
ite = 5 # Para este caso el numero de iteraciones suficientes es 5
Erel_anterior = 1 #Guarda el valor del error anterior 

print("\n"); print("  3. MÉTODO NEWTON - RAPHSON  ")
error3 = np.zeros(ite) # Vector para guardar los errores
for i in range(0,ite):
    xr = xi - fx.subs(r,xi)/fx_prima.subs(r,xi)
    xi = xr

    # Calculo de errores
    # Utilizando la ecuación 5.2 del libro Metodos numéricos para ingenieros 5a ed. pag. 126
    Erel = abs(((-fx_prima2.subs(r,xr)) / (2 * fx_prima.subs(r,xr)))) * (Erel_anterior**2)
    Erel_anterior = Erel
    error3[i] = Erel*100
    print("niter= ", i+1, "xi= ",round(xr,8), "f(xi)= ",round(fx.subs(r,xr),6), "Erel(%)= ", round(Erel*100,10)) # Presentación de resultados

# 2.4. Método de la secante  --------------------------------------------------------------------------------
# utilizando la ecuación 6.6 del libro Metodos numéricos para ingenieros 5a ed. pag. 149

xi = 20
xi_1 =19
Erel_anterior = 1
ite = 7 # Para este caso el numero de iteraciones suficientes es 7
print("\n"); print("  4. MÉTODO DE LA SECANTE  ")
error4 = np.zeros(ite) # Vector para guardar los errores
for i in range(0,ite):
    xr = xi - (fx.subs(r,xi)*(xi_1 - xi)) / (fx.subs(r,xi_1) - fx.subs(r,xi))
    xi = xi_1 #
    xi_1 = xr

    # Calculo de errores
    # Utilizando la ecuación 5.2 del libro Metodos numéricos para ingenieros 5a ed. pag. 126
    Erel = abs(((-fx_prima2.subs(r,xr)) / (2 * fx_prima.subs(r,xr)))) * Erel_anterior**2
    Erel_anterior = Erel
    error4[i] = Erel*100
    print("niter= ", i+1, "xi= ",round(xr,8), "f(xi)= ",round(fx.subs(r,xr),6), "Erel(%)= ", round(Erel*100,10)) # Presentación de resultados

# GRAFICAS ---------------------------------------------------------------------------
 
# Generamos vector de iteraciones
ite1_2 = np.array(range(1,11,1)) # iteraciones para los metodos 1 y 2
ite3 = np.array(range(1,6)) # iteraciones para el método 3
ite4 = np.array(range(1,8)) # iteraciones para el método 4

plt.plot(ite1_2,error1,'r',label='BISECCIÓN')
plt.plot(ite1_2,error2,'b',label='REGULA FALSI')
plt.plot(ite3,error3,'g',label='NEWTON - RAPHSON')
plt.plot(ite4,error4,'y',label='SECANTE')

plt.xlabel('N° de Iteraciones')
plt.ylabel('Error Relativo (%)')
plt.title('Error de cada método')
plt.legend(loc=1)
plt.grid()
plt.show()


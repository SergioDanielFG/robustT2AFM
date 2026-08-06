# Notas metodológicas

Decisiones del proyecto que no se deducen leyendo el código, y los datos que
las sostienen. Cada apartado responde a una pregunta que ya se ha planteado
más de una vez, casi siempre en la dirección equivocada.

---

## 1. El enmascaramiento vuelve sordo al método clásico, no nervioso

Es lo más fácil de invertir de todo el método, y conviene tenerlo por escrito
porque la intuición natural apunta al revés.

La intuición equivocada dice: «si la Fase 1 está contaminada, la carta clásica
saltará a cada rato». La cadena real va en sentido contrario:

1. La contaminación de Fase 1 infla la covarianza combinada `Sp`. En el caso
   Tennessee Eastman su determinante es **16.14 veces** el de la matriz de
   referencia robusta `Sw`.
2. Un límite de control calculado sobre una covarianza inflada es **más
   ancho**.
3. Un límite más ancho avisa **menos veces**, no más.

La contaminación no vuelve nerviosa a la carta clásica. La vuelve sorda. Lo
que hay que ilustrar en cualquier figura o texto es **detección perdida**.

### Los datos

El método clásico produjo **cero falsas alarmas** en todos los escenarios
medidos en este proyecto:

- Tabla 9 del artículo.
- Los cinco fallos de la Tabla 10.
- Las 20 composiciones del análisis de estabilidad.
- Las 12 semillas y las 200 réplicas del escenario de ejemplo del paquete.

No «pocas». Cero. No existe ningún escenario medido en el que el método
clásico produzca falsas alarmas.

### La consecuencia práctica

Ninguna figura, página de ayuda ni viñeta del paquete debe sugerir que el
método clásico genera más falsas alarmas. Y el cero no es una carencia del
escenario elegido que haya que subsanar buscando otro: no hay que salir a
buscar un escenario con falsas alarmas del clásico, porque no existe. Esa
búsqueda es precisamente la intuición invertida del principio.

Esta es la razón de ser de `plot_method_comparison()`, y por eso su pie de
figura empieza por el mecanismo y no por un recuento.

---

## 2. Colorear por verdad y colorear por decisión

La Figura 6 del artículo colorea los lotes por **verdad**: «In-control batch»
frente a «Faulty batch». Esa elección es la que hace el argumento visible, y es
legítima porque el conjunto Tennessee Eastman **viene etiquetado**: se sabe qué
lotes llevan fallo porque el dataset lo dice. La imagen resultante es el
argumento entero: veinte puntos rojos, casi todos por debajo del límite en el
panel clásico. Lotes que se sabe defectuosos y una carta que calla.

Una carta de control en producción no puede colorear así, y no por una
limitación técnica: **qué lotes llevan fallo es justo lo que se quiere
averiguar**. Si esa información estuviera disponible, la carta sobraría.

De ahí la regla del paquete:

- `plot_control_chart()` colorea siempre por **decisión** —por encima o por
  debajo del límite— y nunca usa etiquetas de contaminación. Es la carta de
  producción.
- `plot_method_comparison()` colorea por **verdad** solo si quien la llama le
  aporta esa verdad, mediante el argumento `faulty`. Sin él vuelve a colorear
  por decisión. Pedir el coloreado por verdad sin aportarla es un error, no un
  valor por defecto silencioso.

El coloreado por verdad es un instrumento de **validación**, no de operación.
Sirve para enseñar por qué un método detecta y otro no sobre datos donde la
respuesta se conoce de antemano; no para vigilar un proceso.

Nota histórica: durante un tiempo el paquete no pudo reproducir su propia
Figura 6, porque la única función de carta coloreaba por decisión por diseño y
la figura publicada salía de un script externo. `plot_method_comparison()` con
`faulty` cierra esa brecha, y lo hace sin relajar la regla anterior.

---

## 3. Por qué el gráfico de pesos AFM no colorea por la línea 1/K

`plot_afm_weights()` dibuja la referencia en el peso uniforme `1/K` pero no
colorea las barras según caigan por debajo o por encima. La tentación es
evidente y la razón para resistirla es aritmética.

Los pesos AFM son inversos y normalizados:

```
w_k = (1/lambda1_k) / sum_i (1/lambda1_i)
```

Suman 1 por construcción, luego su media es exactamente `1/K`. Y la
distribución de `1/lambda1` es asimétrica a la derecha, así que su mediana
queda por debajo de su media. Consecuencia: **más de la mitad de los lotes cae
por debajo de `1/K` incluso con una Fase 1 limpia**.

Medido sobre los escenarios del paquete, con K = 30: con 2 lotes realmente
contaminados quedan **17 de 30** por debajo de `1/K`, y con 6 lotes
contaminados también **17 de 30**. Colorear por ese umbral pintaría más de
media calibración sana como sospechosa, y comunicaría algo falso con mucha
firmeza.

Las barras van por tanto en un color neutro, y el resalte queda en manos de
quien llama, mediante `highlight_lowest`, que es una decisión de presentación
y no un juicio del paquete.

Conviene recordar además qué mide el peso: dispersión interna, vía el primer
autovalor. No mide posición. Un lote trasladado en bloque, sin cambiar de
forma, conserva su peso intacto. La ponderación protege `Sw`; no protege el
centro de referencia `mu_r`. En el escenario base, de los seis pesos más bajos
cuatro corresponden a lotes que sí llevan atípicos y dos a lotes limpios: el
peso ordena dispersión, no clasifica lotes.

---

## 4. La cita de la descomposición del T²

La referencia correcta para la descomposición es:

> Mason, R. L., Tracy, N. D., & Young, J. C. (1995). Decomposition of T² for
> Multivariate Control Chart Interpretation. *Journal of Quality Technology*,
> 27(2), 99–108. doi:10.1080/00224065.1995.11979573

Hay dos confusiones recurrentes que conviene dejar zanjadas.

**El orden de autores.** Existen dos artículos distintos y cercanos, con orden
de firma distinto, y ambos son correctos en su sitio:

- Tracy, N. D., Young, J. C., & Mason, R. L. (1992). *Multivariate control
  charts for individual observations*. JQT 24(2), 88–95. Es el de los límites,
  y es el que cita el paquete en la parte clásica.
- Mason, R. L., Tracy, N. D., & Young, J. C. (1995). El de la descomposición.

No hay que «corregir» ninguno hacia el otro. El seguimiento práctico, con el
esquema de cómputo, es Mason, Tracy & Young (1997), JQT 29(4), 396–406.

**La paginación.** La lista de referencias de Montgomery (6.ª ed.) da 109–119
para el artículo de 1995. Es una errata. Los metadatos del editor dan 99–108,
y lo que zanja la cuestión es qué ocupa las páginas 109–119 de ese mismo
número: el artículo de Chang, T. C. y Gan, F. F., *A Cumulative Sum Control
Chart for Monitoring Process Variance*, JQT 27(2), 109–119. Dos artículos no
comparten páginas, de modo que la entrada de Montgomery arrastra la
paginación del artículo contiguo.

Queda pendiente la comprobación definitiva contra la primera página del
artículo impreso; lo anterior se apoya en el depósito bibliográfico del editor
y en la incompatibilidad con el artículo vecino, que es evidencia fuerte pero
indirecta.

---

## 5. Estado de la descomposición del T² sobre `Sw`

La descomposición por variables no está implementada. No es un olvido: la
parte algebraica se traslada sin problema, y la parte distributiva no.

**Lo que sí se conserva.** La descomposición es la factorización de Cholesky
de una forma cuadrática. Para cualquier matriz simétrica definida positiva y
cualquier vector de desviación, la forma se parte en J sumandos no negativos
que suman exactamente el total. Eso vale para `Sw` igual que para la
covarianza muestral combinada, y se comprueba numéricamente sin margen de
duda. Se conservan también la dependencia del orden, y la identidad que da el
término condicional cuando la variable va la última:

```
cond_j = ((Sw^-1 d)_j)^2 / (Sw^-1)_jj  =  T² − T²(todas menos j)
```

donde `d = sqrt(I) * (media del lote − mu_r)`.

**Lo que no se conserva.** Las distribuciones de referencia por término, que
son las que darían p-valores. Las distribuciones exactas de Mason, Tracy y
Young suponen una covarianza Wishart e independiente de la observación
evaluada. `Sw` no lo es, por tres motivos que se acumulan: cada covarianza de
lote sale de un MCD que retiene m\* = 13 de I = 20 observaciones mediante una
selección no lineal; los pesos son función de los datos, así que `Sw` no es
siquiera una combinación lineal fija de matrices Wishart; y sus grados de
libertad efectivos son desconocidos.

**Sobre el umbral por variable.** Existe y es estándar: el corte
`chi²_{alpha,1}` que proponen Runger, Alt y Montgomery (1996), recogido en
Montgomery, Sección 11.3, junto a la ecuación (11.22). Con alpha = 0.01 vale
6.63. No hay que inventar nada.

Lo que hay que decir es más fino: ese corte está derivado suponiendo **Sigma
conocida**. Aquí la covarianza está estimada, y estimada por una vía —MCD por
lote más ponderación dependiente de los datos— cuyo comportamiento
distributivo no está caracterizado en esta configuración. El corte sigue
siendo utilizable como **contorno de referencia**, con una calidad de
aproximación desconocida; lo que no puede es presentarse como una región de
confianza exacta ni traducirse en un p-valor por variable.

> Runger, G. C., Alt, F. B., & Montgomery, D. C. (1996). Contributors to a
> multivariate statistical process control chart signal. *Communications in
> Statistics - Theory and Methods*, 25(10), 2203-2213.

**Lo que sí puede reportarse**, todo álgebra exacta y sin distribución de por
medio: la desviación con signo en unidades de desviación típica, el término
incondicional de cada variable, y el término condicionado a todas las demás.
La lectura útil es el cruce de los dos últimos: alto y alto señala a la
variable como origen; alto e incondicional pero bajo el condicional indica una
variable arrastrada por la correlación; bajo y alto es una violación de
relación, una variable cerca de su objetivo pero fuera de línea dada el resto.

Tres cautelas numéricas medidas sobre el escenario del paquete:

- Los términos incondicionales **no son porcentajes del T²**. Para un lote con
  T² = 117.25 suman 316.21, porque las variables están correlacionadas a
  aproximadamente 0.58.
- El reparto aditivo depende del orden: hay J! órdenes válidos, y en el mismo
  lote la primera variable del orden se lleva el 77.4 % simplemente por ir
  primera.
- La descomposición hereda el sesgo de `mu_r`. Si la Fase 1 traía lotes
  desplazados, el centro de referencia está sesgado y el reparto de culpa se
  hace contra un objetivo equivocado, con total aplomo. La ponderación protege
  `Sw`, no `mu_r`.

**El problema abierto de la elipse.** La idea era una elipse bivariante sobre
la pareja de variables de mayor contribución, centrada en `mu_r` y con la
submatriz 2×2 de `Sw/I`, dividida por I porque lo que se dibuja es la media
del lote. Medida sobre el escenario del paquete, esa elipse con alpha = 0.001
—contorno nominal del 99.9 %— deja **fuera 5 de los 30 centros MCD de Fase 1**,
es decir cubre el 83.3 %. No es un artefacto de la pareja elegida: según cuál
se tome, la cobertura va de 25 a 28 de 30. Y no se arregla bajando alpha: con
0.0027 siguen siendo 25 de 30.

La causa es que `Sw/I` describe la dispersión de una media de lote dentro del
núcleo en control, y no recoge la variación entre lotes ni la contaminación de
la propia Fase 1. Dibujar los centros históricos sobre esa elipse mostraría
cinco lotes de referencia fuera de la referencia. Antes de implementar la
elipse hay que decidir qué covarianza usa; la especificación anterior no
sirve.

---

## 6. Lotes de Fase 1 de tamaño desigual

El límite de control supone un tamaño de lote común I, a través de
`m* = round(I·h)`. La calibración en sí no lo necesita: MCD, los pesos, `Sw` y
`mu_r` se calculan lote a lote y admiten tamaños distintos sin problema. El
límite sí.

En datos de planta los lotes desiguales son lo normal: una tanda se detiene
antes, otra pierde mediciones descartadas. Cuando ocurre, no existe un I
exacto, y hay que elegir uno.

**Se usa el tamaño cuyo `m*` iguala la media de los `m*_k`.** Los grados de
libertad reales del estimador de covarianza con lotes desiguales son
`Σ_k (m*_k − 1)`; igualando esa suma con los `K(m*−1)` de la fórmula publicada
sale `m* = media(m*_k)`. De ahí:

```
m_objetivo = round( mean( round(I_k · h) ) )
I_phase1   = round( m_objetivo / h )
```

**Cuidado, no es lo mismo que redondear la media de los tamaños.** Las dos
formas se confunden con facilidad y difieren de verdad:

- `mean(round(I_k · h))` es la que se deriva de los grados de libertad.
- `round(mean(I_k) · h)` es la que sale de promediar tamaños y convertir
  después.

Sobre 20 000 vectores aleatorios de 30 tamaños entre 10 y 25, las dos
discrepan en el **18.6 %** de los casos. Un ejemplo comprobable a mano, con
h = 0.67 y tamaños 11, 11, 12, 20, 20, 20:

```
m*_k = 7, 7, 8, 13, 13, 13   ->  media 10.167  ->  m* = 10  ->  I_phase1 = 15
round(mean(I_k) · h) = round(round(15.667) · 0.67) = 11     ->  I        = 16
```

En la peor divergencia encontrada (30 lotes de 10 a 25) el `m*` cambia de 11 a
12 y el UCL de 19.81864 a 19.74987, un 0.347 % de diferencia. Con lotes
iguales las dos formas coinciden siempre, así que nada de lo ya medido se
mueve.

El código usa la primera, la derivada. El redondeo va y vuelve sin pérdida:
`round(round(m/h) · h) = m` para todo m entre 5 y 30.

**No se usa el mínimo**, que sería lo «conservador» por instinto. El UCL
decrece de forma monótona con I: con K = 30, J = 4 y h = 0.67, I = 14 da
20.0097, I = 18 da 19.7499, I = 20 da 19.6929, I = 25 da 19.5374. El mínimo da
el límite **más ancho**, y un límite más ancho avisa menos. Sería prudente
frente a las falsas alarmas, que es el riesgo que en este método no existe, a
costa de la detección, que es el riesgo real. Ver el apartado 1.

La diferencia de magnitud es pequeña, alrededor del 1.6 % entre I = 14 e
I = 20. Lo que importaba no era el tamaño del error sino que fuera arbitrario:
la implementación anterior tomaba el tamaño del **primer** lote, y como los
lotes se recorren en orden de aparición, **reordenar las filas del data.frame
cambiaba el límite** sin cambiar un solo dato.

---

## 7. Por qué el ejemplo usa la semilla 20260425

Toda la documentación del paquete —ejemplos, viñeta y README— usa un único
escenario, el de la configuración base del artículo:

```r
simulate_batch_process(
  K1 = 30, K2 = 20, I = 20, J = 4, rho = 0.6,
  outlier_batches_F1 = 6, outlier_rate = 0.20, outlier_shift = 4,
  prop_ooc_F2 = 0.5, shift_ooc = 1.0,
  seed = 20260425
)
```

El escenario anterior tenía un defecto de fondo: ambos métodos marcaban los
mismos 6 de 20 lotes y la razón de determinantes era 1.02. Quien ejecutase el
ejemplo por omisión no veía para qué sirve el método.

**La elección de semilla es deliberada y no es la más favorable.** Sobre 200
réplicas de esta configuración, el paquete promedia 0.824 de detección para el
método propuesto y 0.208 para el clásico. La semilla 20260425 da 8 de 10 y 2
de 10, la réplica individual más próxima a ambas medias a la vez. Las semillas
20260417 y 20260419 dan 10 de 10 al método propuesto: son el resultado más
favorable y por eso están descartadas. No conviene «mejorar» el ejemplo
cambiando a ellas.

**Los conteos del ejemplo no son tasas, y no miden lo mismo que el artículo.**
El artículo publica 0.863 y 0.270 promediando 2000 repeticiones, y las mide
igualando antes ambos métodos a un mismo ARL0 mediante el cuantil empírico de
la distribución nula de cada uno sobre 5000 lotes bajo control. El ejemplo
aplica a cada método su límite operativo. Son cantidades parecidas y no
idénticas: el mismo orden y la misma separación son esperables, la coincidencia
al decimal no. La diferencia entre 0.824 y 0.863, y entre 0.208 y 0.270, queda
explicada por esa hipótesis pero **no está comprobada**, y no debe intentarse
cuadrarla ajustando el ejemplo.

Con delta = 1.0 y alpha = 0.001, además, lo único que separa a los dos métodos
en este escenario es la potencia de detección. Ni en las 12 semillas ni en las
200 réplicas hubo una sola falsa alarma en ninguno de los dos, lo cual
concuerda con el apartado 1 y no es una limitación del ejemplo.

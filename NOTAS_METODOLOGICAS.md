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

## 3. El vocabulario: dos parejas, y por qué son dos

El artículo usa dos parejas de términos y hace bien, porque nombran dos cosas
distintas: «fuera de control» es lo que **dice la carta**, y «lote con fallo»
es lo que **el lote es**. Decisión frente a verdad, que es la misma distinción
que sostiene el diseño de las cartas (apartado 2).

El problema no era tener dos parejas sino tener cinco formas de nombrarlas.
El paquete fija estas:

| Idea | Pareja | Título de la leyenda |
|---|---|---|
| Veredicto de la carta | **Out of control / In control** | `Chart verdict` |
| Condición del lote | **Faulty batch / Fault-free batch** | `Batch condition` |

**Por qué «In control» y no «Under control».** «Under control» es un calco del
español. La literatura de control estadístico dice *in control* y *out of
control*, sin excepción en Montgomery. Además la consola del paquete ya decía
«out of control» e «in control» en `print` y en `summary`, así que esta
elección no obligó a reescribir nada de texto.

**Por qué «Fault-free» y no «Normal».** «Lote en condiciones normales» se
traduce de forma natural por *normal batch*, y en este artículo concreto eso
es una ambigüedad cara: el texto está lleno de distribuciones normales.
*Fault-free* es el complemento exacto de *faulty*, no contiene la palabra
*control* y describe justo lo que ocurre en el Tennessee Eastman, donde el
lote lleva un fallo inyectado o no lo lleva.

**Por qué pueden convivir.** Tres condiciones, y las tres se cumplen: nunca
aparecen en la misma leyenda, porque el color codifica o una cosa o la otra y
nunca las dos; no comparten ninguna palabra discriminante, solo el sustantivo
«batch»; y el título de la leyenda nombra cuál de las dos se está mirando. Esa
tercera condición es la que de verdad lo resuelve: quien ve dos figuras
seguidas no tiene que preguntarse si el mismo color significa lo mismo, porque
la leyenda se lo dice.

### El tercer nivel: la columna `Status` de los datos

La misma decisión hubo que aplicarla un nivel más abajo. `simulate_batch_process()`
generaba la columna `Status` con los valores «Under Control» / «Out of
Control», y esa columna es **verdad de campo**, no veredicto: dice qué lotes
se generaron con fallo.

El síntoma era el ejemplo de la documentación, que construía el argumento
`faulty` así:

```r
faulty <- afm_phase2$Batch[afm_phase2$Status == "Out of Control"]
```

Un argumento llamado `faulty`, que es verdad, filtrando por una etiqueta que
nombra un veredicto. El ejemplo enseñaba exactamente la confusión que las
figuras acababan de eliminar.

`Status` pasa por tanto a **`Faulty` / `Fault-free`**, la misma pareja que la
leyenda de la verdad. Se hizo **antes de publicar los datos**, y ese momento
importa: una vez que los `.rda` estén en Zenodo o en CRAN, cambiar los valores
de una columna rompe el código de quien ya los use.

Comprobado al regenerar con la misma semilla: **las columnas `Var1` a `Var4`
son idénticas byte a byte** —`identical()` sobre los objetos serializados,
diferencia máxima 0—, el identificador de lote no se mueve, y las tres anclas
siguen en 19.69285, 13 y 19.46440, con el T² máximo en 44.8193. No fue un
cambio de datos sino de rótulos.

En la misma pasada se cerró `ContaminationType`, que tenía el nivel `OOC` en
Fase 2: la misma jerga de veredicto, un nivel más abajo. Pasa a **`Shifted`**,
el mismo nombre que ya usaba Fase 1, y los niveles quedan en
`Clean | Outliers | Shifted`.

El motivo es que **es literalmente la misma línea de código**: Fase 1 hace
`mu + shift_contam * sigma_vec` y Fase 2 `mu + shift_ooc * sigma_vec`, y
ambas siguen con el mismo `gen_batch(I, mu_k, Sigma)`. Mismo desplazamiento de
la media, misma covarianza; solo cambia la magnitud, que es un parámetro y no
una categoría.

Sí hay una diferencia real entre los dos casos, pero **es de papel y no de
mecanismo**: en Fase 1 un lote desplazado es contaminación que el método debe
absorber, en Fase 2 es la señal que debe detectar. Ese papel lo determina por
completo la columna `Phase`, así que codificarlo otra vez en
`ContaminationType` era repetir información, y repetirla con dos palabras
distintas inducía a creer que eran dos fenómenos.

Cada columna responde ahora a una sola pregunta: `ContaminationType` a cómo se
estropeó el lote, `Phase` a dónde está, `Status` a si lleva fallo.

### PENDIENTE DE LA TRADUCCIÓN DEL ARTÍCULO

La Figura 6 rotula hoy «In-control batch» / «Faulty batch». Al traducir el
artículo, **esa leyenda debe pasar a «Fault-free batch» / «Faulty batch»**, no
a la traducción literal.

El motivo: la Figura 6 colorea por **verdad**, y «In-control batch» usaría
para la verdad la misma expresión que el resto del artículo usa para el
**veredicto**. La misma frase nombraría dos cosas distintas en figuras
contiguas, que es exactamente la confusión que la separación de vocabularios
evita. Es un cambio de dos palabras y es fácil que se pierda en la traducción,
de ahí que quede anotado aquí.

---

## 4. Por qué el gráfico de pesos AFM no colorea por la línea 1/K

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

**Y por eso su línea de referencia no es burdeos.** En las dos cartas de
control la línea es un **límite**: cruzarla dispara una alarma, y el burdeos
es el color de la alarma en todo el paquete. Aquí `1/K` **no es un límite**:
cruzarlo no significa nada, como acaban de demostrar los 17 de 30. Se dibuja
en pizarra `#2C3E50` para decir que es una referencia y no un umbral. Los dos
colores distinguen dos clases de línea; no es un descuido, y está escrito en
el `@details` de las tres funciones.

Conviene recordar además qué mide el peso: dispersión interna, vía el primer
autovalor. No mide posición. Un lote trasladado en bloque, sin cambiar de
forma, conserva su peso intacto. La ponderación protege `Sw`; no protege el
centro de referencia `mu_r`. En el escenario base, de los seis pesos más bajos
cuatro corresponden a lotes que sí llevan atípicos y dos a lotes limpios: el
peso ordena dispersión, no clasifica lotes.

---

## 5. La cita de la descomposición del T²

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

## 6. Estado de la descomposición del T² sobre `Sw`

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

## 7. Lotes de Fase 1 de tamaño desigual

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

## 8. Por qué el ejemplo usa la semilla 20260425

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

### Decisiones visuales de `plot_afm_weights`

**Color de la línea de referencia.** Pizarra (#2C3E50), no el burdeos (#A02D31)
de las líneas de límite en `plot_control_chart` y `plot_method_comparison`.
Allí la línea es un límite: cruzarla dispara una alarma. Aquí 1/K es la media
aritmética de los pesos por construcción, y cruzarla no significa nada.
Pintarla de burdeos insinuaría un umbral que no existe.

**Por qué no se colorea por 1/K.** Medido con K = 30 sobre `afm_phase1`:
17 de los 30 lotes caen por debajo de 1/K = 0.0333, y solo 6 llevan atípicos.
Colorear por esa línea marcaría como sospechosa a más de media calibración
sana. El mismo recuento sale con 2 lotes contaminados y con 6.

**Ejes transpuestos.** Los identificadores van en el eje y para que queden
horizontales y legibles con los K = 30 lotes del artículo. Ordenar por peso
ya descarta el orden de corrida, así que no se pierde nada al transponer.

**Rango observado de los pesos.** Sobre `afm_phase1`: de 0.0065 (F1_B17) a
0.0900 (F1_B04), un factor de casi 14 entre extremos con solo 6 lotes
ligeramente contaminados.
### Decisiones de diseño de `plot_method_comparison`

**Por qué la figura no lleva pie.** El mensaje central es que la carta clásica
no se vuelve nerviosa bajo contaminación, se vuelve sorda. Esa frase es el
punto entero de la figura y no se imprime debajo: un pie que nadie lee es peor
que ninguno, porque las líneas que sí importan se saltan con él.

**Por qué `scale = "ratio"` es el defecto.** Con estadísticos crudos hacen
falta ejes verticales independientes, y entonces las alturas de los dos
paneles no son comparables. Eso exige una advertencia, y esta figura no tiene
pie donde ponerla. Dividir por el límite propio elimina la necesidad de la
advertencia en vez de repetirla. Con `scale = "T2"` la línea de aviso reaparece
automáticamente: el pie existe solo cuando la figura necesita defenderse de su
propia lectura natural.

**Verificado sobre `afm_phase2`.** El mismo lote F2_B01 aparece en 1,83 del
límite en el panel robusto y en 0,87 en el clásico. F2_B04 igual: 1,79 frente
a 0,87. El conjunto de alarmas del clásico está contenido en el del robusto,
sin excepciones. Es la demostración del enmascaramiento sin una sola fórmula.
### Nombres retirados de la documentación de `ucl_F_adjusted`

**"Way 1".** La ayuda llamaba así al límite analítico, citándolo como si el
artículo usara ese nombre. No lo usa: la Sección 2.5 dice "límite analítico"
y "límite operativo". El nombre venía de una etapa anterior en la que se
numeraban las vías candidatas (analítica, K_eff, bootstrap). Al descartarse
las otras, el nombre perdió sentido y quedó como fósil. Sustituido por la
referencia a la Ecuación (8).

**Hardin-Rocke.** La ayuda declaraba que esa corrección no se aplica. Hardin y
Rocke no aparecen en el artículo ni en su bibliografía; la mención venía del
documento de aclaración del método. Una negación sin referencia y sin
explicación no informa: el lector no sabe qué es lo que no se aplica ni por
qué debería importarle. Retirada.

**Los dos alphas.** `mcd_alpha` (fracción de retención MCD, 0.67) y `alpha`
(tasa de falsas alarmas, 0.001) son parámetros sin relación con nombres
parecidos. Las firmas están publicadas y no se tocan; el aviso va en la ayuda
de ambos argumentos.
### `I_phase1` con lotes de tamaño desigual

**Derivación.** Con lotes desiguales los grados de libertad reales de la
estimación de covarianza son Σ_k (m*_k − 1). Igualando eso con los K(m* − 1)
de la fórmula publicada sale m* = media de los m*_k, y de ahí
I_phase1 = round(m* / h).

**Por qué NO es `round(mean(sizes))`.** Parece lo mismo y no lo es: las dos
formas discrepan en el 18,6 % de 20 000 vectores aleatorios probados. Ejemplo
comprobable a mano con tamaños 11, 11, 12, 20, 20, 20 y h = 0,67: los m*_k son
7, 7, 8, 13, 13, 13, media 10,167, m* = 10, I_phase1 = 15. La otra forma da 11
e I = 16. No lo "simplifiques" de vuelta.

**Por qué no se usa el mínimo.** Sería lo conservador por instinto, pero da el
límite más ancho, y un límite más ancho avisa menos. Sería prudente frente a
las falsas alarmas, que es el riesgo que en este método no existe, a costa de
la detección, que es el riesgo real.

**Orden de aparición.** `unique(data$Batch)` toma el orden de aparición de las
filas. Antes, reordenar las filas cambiaba el UCL sin cambiar un solo dato.
Hay un test que lo impide.

**Verificado sobre `afm_phase1`.** El centro robusto sale en
0,053 / 0,014 / 0,036 / 0,011 con 6 de 30 lotes contaminados. El centro clásico
de `hotelling_classical_calibrate` sobre los mismos datos sale en
0,167 / 0,142 / 0,157 / 0,140: 4,6 veces más lejos del cero verdadero.
### El `sample()` de `simulate_batch_process`

**El bug.** `sample(x, n)` muestrea de `1:x` cuando `x` tiene longitud 1, en
vez de tomar el propio `x`. Es una comodidad histórica de R que se convierte
en error en cuanto un vector se queda con un elemento. Aquí ocurría con
`available_for_outliers`: con `prop_contam_F1 = 0.97` y `K1 = 30` quedan 29
lotes desplazados y uno disponible, y entonces se contaminaba con atípicos un
lote distinto del elegido, dejando las etiquetas `Status` y
`ContaminationType` mal asignadas sin ningún aviso.

**La corrección.** `x[sample.int(length(x), n)]`, que es el idioma recomendado
en la propia ayuda de `sample`.

**Comprobado que no mueve los datasets.** Con `prop_contam_F1 = 0`, que es la
configuración de `afm_phase1`, ninguna rama previa consume el generador y
`available_for_outliers` llega con los 30 lotes, así que las dos formas
consumen la aleatoriedad igual. Los 24 tests de `test-data.R`, que regeneran
los datos desde la llamada con semilla y los comparan con lo distribuido,
siguen en verde.

**La semilla del ejemplo.** `seed = 20260417` coincide con una de las dos
semillas descartadas para el escenario del artículo por dar 10 de 10. Aquí el
escenario es otro (K1=30 con 2 lotes con atípicos y 7% desplazados, K2 con 30%
OOC), así que no es el mismo caso, pero el número coincide y puede confundir a
quien lea ambas cosas.

---

## 9. Qué hace de verdad cada entrada inválida, medido

Las guardas de validación se añadieron sobre una lista de fallos supuestos.
Al medirlos antes de escribir el código, tres de los cuatro resultaron ser
distintos de lo que se creía, y en dos casos peores. Se dejan aquí medidos
porque la intuición sobre ellos ha sido errónea, y porque explican por qué
Fase 1 y Fase 2 necesitan guardas distintas y no una compartida.

### Una columna de texto no hace fallar a `covMcd`: la calibración miente

Se suponía que una variable no numérica hacía fallar a `covMcd` con un mensaje
interno de `robustbase`. **No falla en absoluto.** `data.matrix()` convierte la
columna de texto en códigos de factor, y los niveles se ordenan
alfabéticamente, no numéricamente. Sobre `F1_B01`:

```
Var2 real:          -0,840  -0,258   0,560   0,616   2,979
lo que covMcd ve:        7       4      13      14      20
```

La calibración termina sin avisar y devuelve `mu_r["Var2"] = 10,32` y
`Sw["Var2","Var2"] = 44,09`, frente a ~0,01 y ~1,4 en las columnas honestas.
Todo finito, todo verosímil, todo falso. Es el mismo modo de fallo que motiva
el mensaje de auto-detección de variables: un resultado creíble y equivocado,
que es peor que un error.

### El texto solo es silencioso en Fase 1, y la razón importa

En Fase 2 la misma columna de texto **sí** aborta, con `'x' debe ser numérico`.
La diferencia es qué función toca los datos: Fase 1 pasa por `covMcd`, que
coacciona vía `data.matrix()`; Fase 2 pasa por `colMeans()`, que se niega.

De ahí que las dos fases no lleven la misma guarda por simetría estética, sino
porque fallan distinto: en Fase 1 la comprobación de tipo evita un desastre
silencioso, y en Fase 2 solo mejora un mensaje que ya era ruidoso.

### Un NA en Fase 2 no oculta un lote: destruye el recuento entero

`T2` sale NA, de ahí `is_ooc` sale NA, y como el conteo es `sum(mon$is_ooc)`,
**la suma completa sale NA en vez de un número**. No se pierde la alarma de un
lote: se pierde el total. Ocurría por igual en las dos gemelas, robusta y
clásica, y llegaba impreso al usuario por el camino real,
`run_afm_mcd(compare_classical = TRUE)` → `summary()`.

Este es el motivo de que la guarda de Fase 2 se aplicase a las dos gemelas,
mientras que el rechazo de J = 1 se aplicó solo a la principal: J = 1 es una
entrada absurda de resultado visible, y un NA aparece solo en producción con
efecto invisible.

### Un NA en Fase 1 sí lo absorbe `covMcd`

Es el único de los cuatro que resultó ser **menos** grave de lo que se creía.
Inyectando un NA en `afm_phase1`, `mu_r` sale idéntico a cuatro decimales
(0,0535 / 0,0137 / 0,0356 / 0,0113 en ambos casos): `covMcd` lo absorbe sin
propagarlo. La guarda se mantiene por prevención y por simetría con la gemela
clásica, que ya la tenía, no porque hoy produzca números malos.

---

## 10. Deudas conocidas del paquete

Cosas que están mal o incompletas a sabiendas, con la razón por la que se
dejaron así. No son descubrimientos pendientes de investigar: están medidas y
decididas, y lo que falta es hacerlas.

### `hotelling_classical_calibrate` acepta J = 1 y su gemela no

`calibrate_afm_mcd` rechaza una sola variable con un mensaje que remite a un
gráfico de Shewhart univariante. La gemela clásica no: acepta J = 1 y devuelve
un resultado degenerado.

La asimetría es deliberada por ahora. La función clásica es la implementación
de referencia de la Ecuación (1) del artículo, no una herramienta de
producción, y su resultado degenerado con J = 1 es visible, no silencioso. Aun
así queda coja, y conviene igualarla antes de publicar o documentar por qué no.

### Dos mensajes de error dicen qué falló pero no qué hacer

Son `"The following 'variables' are not numeric: ..."` y `"Batch 'X' contains
non-finite values (NA/NaN/Inf)."`. Incumplen la convención que sigue el resto
del paquete, que es decir qué hacer y no solo qué falló.

Están hoy en las cuatro funciones de calibración y monitoreo porque se
copiaron literales desde `hotelling_classical_calibrate` al añadir las guardas,
buscando que el mismo fallo diera el mismo mensaje en los cuatro sitios.

**Arreglarlos en una sola pasada sobre las cuatro, no de una en una.** El valor
de haberlos copiado literales es precisamente la uniformidad; una corrección
parcial deja una gemela mejor que la otra y destruye lo único que ese copiado
compraba.

### Faltan tests dedicados de dos funciones clásicas

`test-hotelling_classical.R` cubre `hotelling_classical_monitor` y la
delegación desde `run_afm_mcd`. `hotelling_classical_calibrate` y
`hotelling_classical_ucl` no tienen tests propios: se ejercitan de rebote desde
`test-batch_col.R` y `test-run_afm_mcd.R`, lo cual está bien para lo que esos
archivos prueban pero no sustituye a los suyos.

De las dos, la más expuesta es `hotelling_classical_ucl`, porque produce un
ancla publicada (19,464401) y es la que menos cobertura tiene en proporción a
lo que garantiza.

### La guarda de valores no finitos vive dentro del bucle a propósito

Al vectorizar los dos `monitor_*` (los T² se acumulan en vectores y el data
frame se monta una sola vez, en vez de un `rbind` por iteración) quedó la
tentación evidente: sacar la comprobación de finitud a una pasada previa sobre
todo el data frame, que es más rápida y parece equivalente.

**No es equivalente, y la diferencia no rompe nada visible.** El error nombra
un lote, y con la guarda dentro del bucle ese lote es **el primero inválido en
orden de aparición**. Una pasada previa sobre el data frame completo nombraría
el primero por fila, que no es el mismo cuando las filas no vienen agrupadas
por lote. El resultado sigue siendo un error correcto sobre datos inválidos;
solo cambia cuál de los lotes malos se señala, que es justo lo que un
ingeniero usa para ir a buscar el problema.

Hay un test que lo fija reordenando las filas para que cambie el orden de
aparición y comprobando que cambia el lote nombrado. Es el test que detecta
esta optimización si alguien la hace por instinto.

### Las dos gemelas tratan distinto el lote de una sola observación

`hotelling_classical_monitor` tiene una rama para `I == 1`
(`as.numeric(subset_batch[1, ])` en vez de `colMeans()`), y `monitor_afm_mcd`
no la tiene.

La asimetría es real y está pinada por un test, pero no está explicada en el
código ni justificada en el artículo. Con una sola observación el T² de lote
deja de ser el estadístico que el método define —la media de lote es la propia
observación y no hay promediado que reduzca la varianza—, así que lo razonable
sería que ninguna de las dos lo aceptara en silencio, no que una lo trate
aparte. Queda pendiente decidir cuál de las dos formas es la correcta y
alinearlas; hasta entonces, el test impide que la rama desaparezca por
descuido al refactorizar.

### `Rplots.pdf` en la raíz no es un defecto del paquete

`run_examples()` lo escribe en la raíz del proyecto cuando algún ejemplo abre
un dispositivo gráfico, y el siguiente `check()` lo marca con un NOTE:
`Non-standard file/directory found at top level`. Se borra y se vuelve a
comprobar. **No añadirlo a `.Rbuildignore`**, que ocultaría una señal legítima
de `check()`.

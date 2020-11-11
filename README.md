Respositorio de Tesis de Magister Gabriel Esteban Zamora Munoz


## DELINEACIÓN DE ZONAS DE MANEJO Y ROTACIÓN DE CULTIVOS CON RESTRICCIONES DE ADYACENCIA USANDO COLUMN AND ROW GENERATION

**Tesis para optar al grado de Magíster en ciencias de la Ingeniería Industrial y al título de Ingeniero Civil Industrial**





# Resumen

En este trabajo se aborda el problema integrado de delineación de zonas de manejo y planificación de cultivos aplicado en agricultura orgánica. Se busca encontrar un plan de
cultivos que maximice el rendimiento de un terreno agrícola a la vez que se tiene cuidado del impacto medioambiental. Para ello se genera una partición del predio conformada por
zonas de manejo rectangulares relativamente homogéneas en sus propiedades del suelo y se asigna un plan de cultivos que considera la interacción que existe entre distintas especies vegetales cultivadas en zonas adyacentes. Esto ayuda a disminuir el uso de recursos en control de plagas, además de favorecer la salud del suelo.

Con esto en mente, se propone un modelo de programación lineal binario para el problema integrado de delineación de zonas de manejo y planificación de cultivos con
restricciones de adyacencia. En este modelo, se maximizan los ingresos generados por el plan de cultivos a la vez que se respeta el criterio de homogeneidad del suelo en las zonas
de cultivo y las condiciones de adyacencia entre distintas especies.

Las dimensiones del modelo propuesto crecen de forma combinatorial al aumentar el número de divisiones que se admiten en la partición. Por esto, se desarrolla una estrategia
algorítmica basada en Column-and-row generation with column-dependent-rows para resolver el problema. La estrategia de descomposición involucra un problema maestro que
asegura la homogeneidad de las zonas y permite que el plan de cultivos cumpla con los criterios de adyacencia; además de un subproblema que propone zonas que mejoran el
valor de la función objetivo del problema maestro.

La estrategia algorítmica fue implementada en JuMP, un lenguaje de modelado para optimización alojado en Julia. Los resultados del conjunto de instancias ensayadas muestran
la relevancia del uso de la estrategia de descomposición a la hora de resolver el problema.

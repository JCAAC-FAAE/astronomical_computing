[Journal of Computational Astronomy & Astronomical Computing (JCAAC)](https://jcaac.org)

# SECTION: Astronomical Computing

[`No. 1, Nov. 2024, 57–63`](https://federacionastronomica.es/index.php/the-journal/archive/135-contents/613-astronomical-computing)

## Presentación. Cálculo del Día Juliano
> En esta sección se desarrollarán algoritmos que cubren desde las operaciones más básicas, como obtener el día Juliano en este primer capítulo, hasta otros más complejos como calcular la posición de los cuerpos celestes y eventos astronómicos como eclipses y muchos otros. El lenguaje de programación elegido para presentar el código es Java, pero se espera que el lector escriba el código en su lenguaje preferido, aprovechando las ventajas del repositorio de software para evitar errores en la transcripción de las operaciones matemáticas.

## Presentation. Computing the Julian Date
> In this section we will develop algorithms to cover from the most basic calculations, like obtaining the Julian day in this first chapter, to more complex ones that will allow us to calculate the positions of the celestial bodies and astronomical events like eclipses and many others. The programming language chosen to present the code is Java, but the readers are expected to write this code in their preferred language, taking advantage of the software repository to avoid errors in the transcription of the mathematical operations.

<br><br>

---

[`No. 2, Mar. 2025, 87–103`](https://federacionastronomica.es/index.php/the-journal/archive/)

## Cálculo de ángulos de orientación de la Tierra y transformaciones de coordenadas
> En esta entrega de la sección discutimos el cálculo de la diferencia TT-UT1, que resulta básico para cualquier cálculo de fenómenos, ya que relaciona el tiempo en el que se expresan las efemérides planetarias con el tiempo asociado a la rotación de la Tierra y, por tanto, a la posición de un observador sobre la misma. El Tiempo Sidéreo Aparente, que es esencialmente el ángulo de rotación terrestre, se relaciona de manera directa con UT1 y con la ecuación de los equinoccios. Además, se discute el cálculo de los parámetros de precesión y nutación, y de la oblicuidad de la eclíptica, necesarios para diversas y necesarias transformaciones entre sistemas de coordenadas astronómicas, que se tratan al final del artículo.


## Computation of Earth Orientation Angles and Coordinate Transformations
> In this article we discuss the calculation of the TT-UT1 difference, which is fundamental for any calculation of astronomical phenomena, as it relates the time in which planetary ephemerides are expressed to the time associated with the Earth's rotation and, therefore, to the position of an observer on its surface. Apparent Sidereal Time, which is essentially the Earth's rotation angle, is directly related to UT1 and the equation of the equinoxes. Additionally, the calculation of precession and nutation parameters, as well as the obliquity of the ecliptic, is discussed, as these are necessary for various essential transformations between astronomical coordinate systems, which are addressed at the end of the article.

<br><br>

---

[`No. 3, Nov. 2025, 65–78`](https://federacionastronomica.es/index.php/the-journal/archive/contents/706-calculo-de-efemerides)

## Cálculo de efemérides
> En este número veremos cómo reducir las coordenadas de un objeto para obtener la posición del mismo en el
cielo para un observador situado sobre la superficie terrestre. Veremos el caso más sencillo en el que suponemos
que ya tenemos la posición eclíptica geocéntrica del objeto referida al equinoccio medio de la fecha, de manera
que no es necesario aplicar la corrección de precesión, vista en entregas anteriores. También se introducirá cierta
cantidad de código adicional para tener bien organizado los resultados de los cálculos, incluyendo los resultados
de operaciones potencialmente iterativas, como obtener los instantes de salida y puesta del objeto. En el caso de
cuerpos estáticos, como las estrellas, esta iteración no es necesaria.

## Ephemeris Computation
> In this article we see how to reduce the coordinates of an object to obtain its position in the sky for an observer
located on the Earth’s surface. We will examine the simplest case, in which we assume that we already have the
geocentric ecliptic position of the object referred to the mean equinox of the date, so it is not necessary to apply
the correction for precession, as discussed in previous articles in this Section. We also introduce some additional
code to keep the calculation results well organized, including the results of potentially iterative operations, such as
obtaining the rise and set times of the object. In the case of static bodies, such as stars, this iteration is not necessary.

<br><br>

---

[`No. 4, Dec. 2025, 45–60`](https://www.federacionastronomica.es/index.php/the-journal/archive/contents/751-efemerides-precisas-del-sol-y-la-luna)

## Efemérides precisas del Sol y la Luna
> En este número presentamos algoritmos para el cálculo preciso de la posición del Sol y la Luna para un observador en la Tierra, utilizando el procedimiento de reducción explicado en el número anterior. Para alcanzar la mejor consistencia posible, respecto de la última integración numérica disponible del JPL, se han modificado algoritmos clásicos, de manera notable en el caso de la Luna. Las discrepancias máximas que se obtienen son de 2" para la posición geocéntrica del Sol, y en torno a 15" para la de la Luna, durante al menos 4000 años alrededor del año 2000. También se presentan cambios que se han introducido en ficheros publicados con anterioridad, para corregir algunos problemas.

## Accurate Ephemeris of the Sun and the Moon

> In this article we will see how to compute the accurate positions of the Sun and the Moon for any observer on Earth, using the reduction procedures presented in a previous article in this Section. To achieve the maximum consistency, with respect to the latest JPL numerical integration, some classic algorithms has been modified, especially for the Moon. The maximum discrepancies found are of 2" for the geocentric posicion of the Sun, and around 15" for that of the Moon, during at least 4000 years around the year 2000. Some corrections to files published previously are also presented, solving different problems identified.

---

[`No. 5, Mar. 2026, xx–xx`](https://federacionastronomica.es/index.php/the-journal/archive/contents/751-efemerides-precisas-del-sol-y-la-luna)

## Efemérides precisas mediante integración numérica. Aplicación a los planetas y al asteroide Apophis

> En esta entrega introduciremos un método extraordinariamente preciso para el cálculo de las efemérides
de planetas y cuerpos menores: la integración numérica. Veremos una implementación simplificada que
resulta en un código eficiente y de extensión reducida, que gracias a la gran velocidad de los ordenadores
actuales, puede ejecutarse con Java en menos de un segundo. Otros métodos tradicionales, utilizados aún
hoy en la mayoría de las situaciones, son mucho menos precisos y más extensos. Por supuesto, también
son más veloces, pero en la mayoría de las situaciones esa rapidez ya no es relevante en la práctica.
Una gran ventaja de este método es que permite calcular con gran precisión la posición de asteroides
y cometas arbitrarios, sin aumentar significativamente la extensión del código. Analizaremos en detalle
los casos de Apophis y 2024 YR4. Apophis es un objeto potencialmente peligroso para la Tierra, con un
tamaño de unos 400 m, que pasará muy cerca de la Tierra el 13 de abril de 2029. Por su parte, 2024 YR4
podría colisionar con la Luna en 2032.

## Accurate ephemerides with numerical integration. Application to the planets and the asteroid Apophis
> In this article we will introduce an extraordinary accurate method to compute the ephemerides of planets
and minor bodies: the numerical integration. We will see a simplified implementation that results in an
efficient and short code, that due to the high speed of current computers, can be executed with Java in less
than a second. Other traditional methods, still used today in most situations, are less accurate and longer.
Of course, they are also faster, but in most situations that speed is no longer relevant in practice.
An important advantage of this method is that it allows to compute, with great accuracy, the position
of any asteroid or comet, without increasing significantly the extension of the code. We will analyze in
detail the cases of Apophis and 2024 YR4. Apophis is a potentially dangerous body for the Earth, with a
size of around 400 m, that will pass very close to the Earth on April, 13, 2029. On the other hand, 2024
YR4 may collide with the Moon in 2032.


<!-- plantilla para incluir siguientes números en el mismo repositorio

<br><br>

---

[`No. 1, Nov. 2024, 57–63`](https://federacionastronomica.es/index.php/the-journal/archive/135-contents/613-astronomical-computing)

## Título artículo
> Texto.

## Article Title
> Text.
-->

<br><br>

---

:copyright: **Tomás Alonso Albi** 2024-2026. Licensed under the EUPL.



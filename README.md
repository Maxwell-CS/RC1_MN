# RC1 — Puntos de operación de un sistema de bombeo

Resolución de Casos 01 · Métodos Numéricos · UTEC · 2026-2
**Sección 9 — Grupo 3**

| Integrante | Participación |
|---|---|
| Maxwell Lupo Gregorio Collazos Solis (líder) | 100% |
| Pablo Nicolás Dávila Wong | 100% |
| Sebastian Manuel Avila Torres | 100% |
| Matias Fernando Castillo Quincho | 100% |
| Melvín Billy León Capcha | 100% |

## Contenido

```
Plantilla_RC1_2026-2_Grupo3.mlx   entregable (abrir con MATLAB R2023b+)
enunciado/                        PDF del caso y plantilla en blanco
src/                              código de las celdas, en texto plano
```

El `.mlx` se ejecuta completo desde un entorno limpio (`clearvars` en la
primera celda). Los archivos de `src/` son el mismo código en `.m` para
poder revisarlo y hacerle diff sin abrir MATLAB.

## Problema

Una bomba centrífuga eleva agua hacia un tanque. El punto de operación es
la intersección entre la altura que entrega la bomba y la que exige la
instalación. Como las pérdidas dependen del factor de fricción y éste del
caudal, el sistema es no lineal y acoplado:

```
F1(Q,f) = H0 + aQ - bQ² - Δz - C·B(f)·Q² = 0     (energía, en metros)
F2(Q,f) = f^(-1/2) + 2·log10 A(Q,f)     = 0     (Colebrook, adimensional)
```

con `C = 8/(gπ²D⁴)`, `B(f) = fL/D + K_T` y `A(Q,f) = ε/(3.7D) + 2.51/(Re√f)`.

## Resultados

Dos equilibrios en el dominio turbulento:

| | Q [m³/s] | f | Hb = Hs [m] | V [m/s] | Re | F'(Q) | estabilidad |
|---|---|---|---|---|---|---|---|
| P1 | 0.0076743 | 0.021776 | 32.081 | 0.244 | 48 856 | +222.92 | inestable |
| P2 | 0.0543246 | 0.016197 | 35.231 | 1.729 | 345 841 | −221.75 | **estable** |

Comparación de métodos:

| Método | Semilla / intervalo | Iter. | Eval. F | Iter. internas | Estado |
|---|---|---|---|---|---|
| Bisección | [0.00765, 0.00781] | 15 | 17 | ~125 | converge a P1 |
| Bisección | [0.05419, 0.05435] | 15 | 17 | ~125 | converge a P2 |
| Newton | (0.005, 0.030) | 5 | 11 | 0 | converge a P1 |
| Newton | (0.060, 0.020) | 5 | 11 | 0 | converge a P2 |
| Punto fijo | (0.010, 0.030) | 299 | 598 | 0 | converge a P2 |
| Punto fijo | (0.060, 0.020) | 307 | 614 | 0 | converge a P2 |
| Punto fijo | (0.005, 0.030) | — | — | — | dominio inválido en k=0 |

Newton alcanza ambas raíces con convergencia cuadrática. El punto fijo
sólo alcanza P2: en P1 el radio espectral de `JG` es 11.43 > 1, de modo
que la raíz existe pero la iteración la repele. **La estabilidad de una
iteración y la estabilidad física de un equilibrio son cosas distintas.**

Sensibilidad (peor caso lineal, con Δε = 10% y ΔD = 1%):

| | ∂Q/∂ε | ∂Q/∂D | ΔQ/Q | dominante |
|---|---|---|---|---|
| P1 | +0.2269 | −0.008551 | 0.236 % | D (94 %) |
| P2 | −26.790 | +0.34739 | 1.501 % | D (85 %) |

Domina el diámetro porque interviene en cuatro dependencias (`C ~ D⁻⁴`,
`L/D`, `ε/D` y `Re`), mientras la rugosidad sólo entra en `A` con efecto
logarítmico.

Barrido de la altura estática, Δz ∈ [28, 36] m:

- **Δz < 30.186 m** — una solución turbulenta
- **30.186 < Δz < 34.590 m** — dos soluciones
- **Δz > 34.590 m** — ninguna

Los dos cambios tienen causas distintas: en Δz = 34.590 m las ramas se
funden en un cruce **tangente** (`F'(Qc) = 0` en `Qc = 0.030975`), que una
malla de 0.1 m no puede localizar porque `F` no cambia de signo; en
Δz = 30.186 m la raíz menor **sale del dominio** por el borde `Re = 4000`.

## Métodos

Bisección con iteración interna de Colebrook, Newton para el sistema
acoplado con Jacobiano analítico (paso obtenido de `J·s = −F`, sin formar
`J⁻¹`) y punto fijo simultáneo con análisis de normas y radio espectral.
Las funciones del laboratorio del curso (`newton2v`, `newtonVec`,
`pfijos`) van adjuntas íntegramente en la última celda del `.mlx`.

## Alcance

Modelo didáctico. No sustituye la verificación de cavitación (NPSH),
eficiencia, potencia al eje ni los límites operativos del fabricante.

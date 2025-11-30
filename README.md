# Predator-prey model for investment py implementation
A Python implementation of a 3-dimensional Lotka-Volterra model applied to stock market acquisitions which was already in matlab. Includes numerical simulations of stability and Hopf Bifurcation analysis. Developed as an undergraduate thesis project




# Predator-Prey for Investment Dynamics

## Overview
This repository is the codebase  from my **undergraduate final year project**. The project explores the application of biological dynamical systems to financial markets, modeling company shares as an ecosystem.

Based on the research by **Addison, Bhatt, and Owen (2016)**, I implemented a modified three-dimensional Lotka-Volterra model where large investors act as "predators" and smaller competing companies as "prey." The study investigates how varying financial parameters (such as risk premium and price volatility) influence market stability.

## Key Features
* **ODEs Implementation:** Solves the non-linear system of differential equations using `scipy.integrate.solve_ivp`.
* **Stability Analysis:** Simulates both stable equilibrium and unstable limit cycle scenarios based on specific parameter sets.
* **Bifurcation Analysis:** Numerically reproduces the **Hopf Bifurcation** diagram to identify the exact threshold where market stability breaks down.
* **Visualization:** Generates publication-quality plots using `Matplotlib` and `Seaborn`.

## The Model
The system is defined by three differential equations representing:
1.  **$X_1, X_2$:** Relative share prices of two competing "prey" companies.
2.  **$Y$:** Relative share price of the "predator" (investing company).

### Governing Equations
The model is governed by the following system of three non-linear ordinary differential equations:

$$
\begin{aligned}
\frac{dX_1}{dt} &= s_1 X_1 \left(1 - \frac{X_1}{K_1}\right) - m_{12} X_1 X_2 - v_1 X_1 Y r_1 \\
\frac{dX_2}{dt} &= s_2 X_2 \left(1 - \frac{X_2}{K_2}\right) - m_{21} X_1 X_2 - v_2 X_2 Y r_2 \\
\frac{dY}{dt}   &= -\mu Y + c_1 v_1 X_1 Y r_1 + c_2 v_2 X_2 Y r_2
\end{aligned}
$$

Where the predator switching functions ($r_1$ and $r_2$) for $n=1$ are defined as:

$$
r_1(X_1, X_2, Y) = \frac{X_1}{X_1 + a_2 X_2 + b_2 Y}, \quad r_2(X_1, X_2, Y) = \frac{X_2}{X_2 + a_1 X_1 + b_1 Y}
$$


The predator utilizes a "switching function" to invest in the prey species that is most abundant (profitable), creating complex non-linear dynamics.
 
 

### 2. Hopf Bifurcation
A numerical analysis of the conversion rate parameter ($c_2$). The diagram in code clearly shows the transition from a stable fixed point to a limit cycle (fork) at the critical value $c_2 \approx 0.228$.

 
## Requirements
* Python 3.x
* NumPy
* SciPy
* Pandas
* Matplotlib
* Seaborn

 

# Numerical methods computerlaboration 1 --- SF1544

This computer laboration is one of three which are part of the Numerical methods course taken on the engineering physics program at KTH university.
The laboration is examined with an on-site examination, and an assessment of our results and code-documentation (Matlab code) in the form of a submission. 
This readme-file is a shorter version of the actual instructions of the laboration that we had answered. The answers are briefly dockumented in the code along with the required plots

Special thanks to my laborationpartner [Mohammad Kanjo](https://github.com/MohiMad).

## Introduction
This project explores a quarter-car suspension system using numerical methods along with knowledge in mechanics,
differential equations, linear algebra, and multivariable calculus to structure and solve ordinary differential equations (ODEs) 
associated with the system's dynamics. 

A vehicle's suspension system connects its chassis with springs and shock absorbers, which is critical for comfort and road handling. 
that connect a vehicle to its wheels. The specific system that we have analysed is the following, also refered to as a quarter-car modell,
which is the suspentionsystem connected to one wheel:\\

<img width="232" alt="image" src="https://github.com/user-attachments/assets/08a20a89-d879-4c76-9d47-8b6bd1f4044a" />

The goal is to analyze the stability and accuracy of the suspension system, and optimize parameters for improved ride comfort and safety.

### Parameters

As observed by the image above we have in order: <br>
$z1, m1, k1, c1$: Chassis displacement, mass, spring constant and damping constant. <br>
$z2, m2, k2, c2$: Wheel displacement, mass, related spring and damping constant.

Furthermore, h is the vertical position of the road surface, which is given by

$$
h(t) = \begin{cases} 
      \frac{H}{2}(1 - cos(\frac{2 \pi u t}{L})), & t\leq L/u \\
      0 & t\geq L/u 
   \end{cases}
$$

The values of the parameters are given by the follwoing table.
| Parameter | Value     | Units  |
|-----------|-----------|--------|
|   $m_1$   | 475       | kg     |
|   $m_2$   | 53        | kg     |
|   $k_1$   | 5400      | N/m    |
|   $k_2$   | 135000    | N/m    |
|   $c_1$   | 310       | Ns/m   |
|   $c_2$   | 1200      | Ns/m   |
|   $u$     | 65/3.6    | m/s    |
|   $H$     | 0.24      | m      |
|   $L$     | 1         | m      |



To describe the equations governing the suspentionsystem we use analytical mechanics to:
1. Choose our generalized paramters z1 and z2 and their respective expressions.
2. Write expressions for the mechanical energy T and potential energy V.
3. Derrive L = T - V with respect to both of our generalized paramters, as well as their **first** derrivatives.
4. Write the related Euler-lagrange equations.

The Euler-lagrange equations give us the final expressions.

$$
m_1 \ddot{z_1} + c_1(\dot{z_1} - \dot{z_2}) + k_1(z_1 - z_2) = 0
$$

$$
m_2 \ddot{z_2} + c_1(\dot{z_2} - \dot{z_1}) + c_2(\dot{z_2} - \dot{h}) + k_1(z_2 - z_1) + k_2(z_2 - h) = 0
$$

We can finally answer the following the questions related to the laboration.

## Question 1
Write the equations governing the movement of the equations as a first order differential equation with respect to time, in the form:

$$
\frac{d\mathbf{v}}{dt} = \mathbf{A}\mathbf{v} + \mathbf{g}(t)
$$

where

$$
\mathbf{v} = [z_1, z_2, \dot{z_1}, \dot{z_2}]^T
$$.

## Question 2
Write code that solves the differential equation problem from question 1 with the help of: <br>
- Forward Euler
- Matlabs built-in function ODE-45

Assume that the vertical velocities and positions are zero initially.

a) Plot the numerical solution for z1 and z2 as a function of time using ode45 with relative tolerance 10−6. What will be the maximum deviation for z1 and z2 and what does this mean for the passengers?

b) Visualize the size of the time steps for ode45 and the size-change over time. Explain the changes in the time steps over time by studying the plots of the numerical solution and the time steps together.


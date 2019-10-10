Sample Matlab and Mathematica code for effective velocity and diffusivity calculation

This code accompanies the preprint "Renewal reward perspective on linear switching diffusion systems". 

The Mathematica sample code provides the calculation of effective transport properties (velocity and diffusivity) using the proposed renewal rewards framework for a cooperative model with at most 2 kinesin-1 motors (i.e., with three dynamic states or behaviors). The inputs include the probability transition matrix and the vectors of moments of rewards (state durations and spatial displacements) for each state. The final section of the code implements the calculation procedure for a range of load forces of the model motor.

The Matlab sample code provides calculation of effectve velocity and diffusivity in tug-of-war models of cargo driven by maximumtwo kinesin-1 and maximum two conventional dynein motors (i.e., switching between nine dynamic states or behaviors). The main code is "tug_war.m", and the supporting Matlab functions calculate the probability transition matrix, the rates of binding and unbinding of the motors to microtubules, the cargo forces and velocities. This sample code allows for diffusion in the detached state and can be easily adapted to model tug-of-war interactions of motors with different kinetic parameters.



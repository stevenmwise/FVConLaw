# FVConLaw

This MATLAB code is designed to approximate solutions of the inviscid Burger equation.

It is a modified version of a code developed by Justin Dong (https://jdongg.github.io),
released as numCL (https://github.com/jdongg/numCL). In particular, the numerical
fluxes and the exact solution functions were coded by Justin. He has a very nice blog
post on his website describing the basic ideas of numerical methods for conservation
laws.

This is a teaching code that allows students to explore various solver combinations,
including SSP-RK-2, SSPRK-3 time stepping, linear and WENO5 function reconstruction
procedures, and several different numerical flux functions.
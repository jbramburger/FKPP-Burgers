This repository contains the MATLAB files to reproduce the data and figures from "The speed of traveling waves in a FKPP-Burgers system" by Jason J. Bramburger and Christopher Henderson (2020).

Computations use YALMIP version R20190425 to translate the SOS constraints into semidefinite programs and these are solved using Mosek version 9.0. These programs are publicly available at:

YALMIP: https://yalmip.github.io/download/
Mosek: https://www.mosek.com/downloads/

To reproduce data the scripts assume that YALMIP files are stored in the folder 'YALMIP-master' and Mosek files in the folder 'mosek'. 

In all scripts the parameter 'd' represents the degree of the auxiliary polynomial, denoted H in the manuscript. Refer to the appendix of the manuscript for the inequalities that are being implemented to obtain the bounds in Table 1. The scripts in this repository do the following:

- Inviscid_Upper: Bounding proceedure for the upper bound on the minimum wave speed of the inviscid FKPP-Burgers equation.

- Inviscid_Lower: Bounding proceedure for the lower bound on the minimum wave speed of the inviscid FKPP-Burgers equation.

- Viscous_Upper: Bounding proceedure for the upper bound on the conjectured minimum wave speed of the viscous FKPP-Burgers equation.

- Viscous_Lower: Bounding proceedure for the lower bound on the conjectured minimum wave speed of the viscous FKPP-Burgers equation.

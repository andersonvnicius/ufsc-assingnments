# Welcome to finte element analysis

In this section there are the codes that i've used to solve simple finite element analysis problems during graduation

## Spring Elements

The most basic mechanical element there is, it only accepts traction or compression loads along it's length, it's displacement is given by hooke's law

F = k*x

where:
F: load applied
k: stiffness array of the element
x: displacement caused by the application of the load

The stiffness aray of an sping element is of the form [[k, -k],[-k, k]]


## Beam elements

TBD

## Truss elements

TBD


# Requirements

Python 3

numpy

matplotlib

pandas


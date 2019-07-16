#pragma once
// Version 1 alpha demo
// This file related to aCRD/Main.cpp
// This file defines several parameters
// Further dll Prototype for parameters
// The main input memebers use the following units
// Atomic radius : nm
// Crystal structure : nm
// Diffusivity: m^2/s - nm^2/s
// Time : s
// Energy : keV

#ifndef PARAMETER_H
#define PARAMETER_H

// Unit in nm
#define THRES_1 1.00 // Trigger of E02 during protecting, with another HPD
#define THRES_2 0.50 // Trigger of E10 during protecting, with another walker. 4 x single radius
#define THRES_3 1.50 // Minimum distance in E11 releasing, within final walkers
// Unit in s
#define LIMIT_1 5.0e-8 // E10 Holdup destination time, stamp for buildstamp + LIMIT_1

// This is only for first initialization estimation
#define INITIAL_EVENTS_TIME 0.001

#define FE_RADIUS 0.126
#define FE_BCC_A 0.28665 // Assume to be single vacancy
#define INTERST_MODIFIER 1.00
#define VACANCY_MODIFIER 1.25

#define INIDOMAIN_MODIFIER 6.0 // For sqrt ( Factor * Dv * t). This will trigger initial FE 900, need to change TH2

#endif
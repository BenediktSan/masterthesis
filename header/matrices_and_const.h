#include <iostream>
#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include <iomanip>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <complex>
#include<string>
#include<sys/types.h>
#include<sys/stat.h>


#include "../header/functions.h"


using cplx = std::complex<double>;
using namespace std;
using Eigen::MatrixXd;
using Observable = Eigen::MatrixXcd;

// Constants

cplx cplx_i = cplx(0,1.);


Observable build_local_Splus(const float j);


const double conversion_factor = 1;

//Lande faktoren dabei recht willkürlich als 2
const double hbar = 1.05457182 * pow(10., -34.);
const double const_e = 1.602176634 * pow(10., -19.);
const double mu_theo = 9.2740100783; //in 10^-24 J/T


//this one is used
//const double mu_b =  mu_theo * pow(10.,-24.) * 1./hbar  ; //hbar/µs*T
const double mu_b = mu_theo;
//Constants for coupling constand d


//vacuum permeability
const double mu_0 =1.25663706127 * pow(10,-6);//* pow(10,4);//kg*Angstrom/( C^2)      // 1.25663706127 *  pow(10.,-6.); // N/A^2

//gyromagnetic ratio of fluor
//const double gamma_i = 25.181  ; // 10^7 rad / (T s)

//constants from book 
const double gamma_i = 251.662 * pow(10,6) ;// rad / µsT   40.077;//40.069244;//40.077 ;// 1/(T µs)
//##### IMPORTNANT NOTE: The following lattice constant was caclulated for room temperature
//However all experiments were done at approcimately 4K, which still hlds for infinite temperature in spin direction distribution
//Therefore I have to correct the lattice constant in the post processing of the data
//The correct value at 4K is 2.72325 Angstrom
const double lattice_const = 2.7315; // 5.451;// 2.724;//5.451; //in Angstrom

//B field used 
//is alloed to be not normalized. All calculations only unse B_field_direction.norm() ( for example coupling constant)
//needs to be changed. therefore no const

//############ THIS LINE IS IMPORTANT ######################
Eigen::Vector3cd B_field_direction(0,0,1);



//Prefactor for that is used for deriving T:max in calculation of the correlation.
// is only allowed to be <=0.5
const double T_max_prefactor = 0.5; 

//seed for RNG
//seed is being kept. Generator hast to be scratched due to use of MPI

double seed = 5;
mt19937 gen(seed); // Mersenne Twister generator initialized with seed





const Observable SIGMA_0{   { cplx(1.,0.), cplx(0.,0.) }, 
                                    { cplx(0.,0.), cplx(1.,0.) }    };

//const Observable SIGMA_X{   { cplx(0.,0.), cplx(1.,0.) }, 
//                                    { cplx(1.,0.), cplx(0.,0.) }    };
//
//const Observable SIGMA_Y{   { cplx(0.,0.), cplx(0.,-1.) }, 
//                                    { cplx(0.,1.), cplx(0.,0.)  }   };
//
const Observable SIGMA_Z{   { cplx(1./2.), cplx(0.,0.)  },                    
                                    { cplx(0.,0.), cplx(-1./2.,0.) }   };

Observable S_p = build_local_Splus(1./2.);
Observable S_m = S_p.adjoint();


const Observable SIGMA_X = 1./2. * ( S_p + S_m);

const Observable SIGMA_Y = 1./(2. * cplx_i) * ( S_p - S_m);


// MAtrices for Spin S = 1

const Observable SIGMA_0_S1{   { cplx(1.,0.), cplx(0,0.), cplx(0.,0.)  }, //For S = 1
                            { cplx(0,0.), cplx(1.,0.), cplx(0,0.)  },
                            { cplx(0.,0.), cplx(0,0.), cplx(1.,0.)  } };

//const Observable SIGMA_X_S1{   { cplx(0.,0.), cplx(sqrt(2.),0.), cplx(0.,0.)  }, //For S = 1
//                            { cplx(sqrt(2.),0.), cplx(0.,0.), cplx(sqrt(2.),0.)  },
//                            { cplx(0.,0.), cplx(sqrt(2.),0.), cplx(0.,0.)  } };
//
//const Observable SIGMA_Y_S1{   { cplx(0.,0.), cplx(0.,-1 *sqrt(2.)), cplx(0.,0.)  }, //For S = 1
//                            { cplx(0.,sqrt(2.)), cplx(0.,0.), cplx(0.,-1* sqrt(2.))  },
//                            { cplx(0.,0.), cplx(0.,sqrt(2.)), cplx(0.,0.)  } };
//                        
const Observable SIGMA_Z_S1{   { cplx(1.,0.), cplx(0.,0.), cplx(0.,0.)  }, //For S = 1
                            { cplx(0.,0.), cplx(0.,0.), cplx(0.,0.)  },
                            { cplx(0.,0.), cplx(0.,0.), cplx(-1.,0.)  } };


Observable S_p1 = build_local_Splus(1);
Observable S_m1 = S_p1.adjoint();


const Observable SIGMA_X_S1 = 1./2. * ( S_p1 + S_m1);

const Observable SIGMA_Y_S1 = 1./(2. * cplx_i) * ( S_p1 - S_m1);
                        

// NOTE: to run the supplementary material, use Jupyter notebook. See README

// Front propagation based overhang filer
//
// Part of supplementary material of "Overhang control in topology optimization: 
// a comparison of continuous front propagation-based and discrete layer-by-layer
// overhang control", E. van de Ven, R. Maas, C. Ayas, M. Langelaar, F. van Keulen,
// 2019
//
// Code by Emiel van de Ven, 2020
//
//
// Disclaimer:                                                             
// The author reserves all rights but does not guarantee that the code is  
// free from errors. Furthermore, the author shall not be liable in any   
// event caused by the use of the program.       
//                         
// emiel@emielvandeven.nl

#include <iostream>
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>
#include <limits>
#include <queue>
#include <utility>
#include <algorithm>
#include <chrono>
#include <functional>
#include "MatrixXX.h"

class FPAMFilter {
    public:     
        // constants:    
        double kappa; // parameter that influences length of black-white transisiton (EQ. 6)
        double v_void; // void speed in front. Only relevant for speedFunctionType=0 (EQ. 11)
        double a_oh; // overhang angle. For square elements this can only be pi/4
        double vx; // propagation speed perpendicular on build direction

        int nelx; // number of elements in x direction
        int nely; // number of elements in y direction
        int nv;   // total number of elements
        int max_dep; // number of elements the arrival time for an element is dependent on
        int speedFunctionType; // 0: old, 1: improved

        MatrixXX C_p; // projection matrix
        MatrixXX build_dir; // build direction. For this script it can only be [0,1]

        std::vector<int> diff_mat, dm_x, dm_y; // arrays to calculate indices of neighbouring elements

        // vectors for front propagation:
        std::vector<bool> Accepted; // array which determines if an element is Accepted (true) or not (false)
        std::vector<double> T; // array with arrival time for each element
        // priority queue to quickly get the element with the lowest arrival time. The STL priority queue does not support
        // updates of values already in the queue, so its would be faster to implement a custom binary tree. However this is very convenient.
        std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<std::pair<double, int>> > pqueue; 

        // vectors for sensitivity calculations:
        std::vector<int> order; // the order in which elements are Accepted (sensitivities are calculated in opposite order)
        std::vector<int> prev;  // lists for each element on which two other elements its arrival time is based
        std::vector<double> dep_der; // partial derivative of arrival time to the arrival time of the two elements its arrival time is based on
        std::vector<double> dep_rho; // partial derivative of arrival time to the density value of the element

        // initialize the object: allocate memory and set constant values
        FPAMFilter(int _nelx, int _nely, int _speedFunctionType, double _kappa, double _v_void): 
            nelx(_nelx), nely(_nely), speedFunctionType(_speedFunctionType), kappa(_kappa), v_void(_v_void) 
        {
            nv = nelx*nely;
            std::cout << speedFunctionType << std::endl;
            // constants
            a_oh = M_PI/4.0; 
            vx = 1.0/std::tan(a_oh);
            build_dir = MatrixXX {2,1,{0,1}};
            C_p =  MatrixXX(2,2,true) - build_dir.dot(build_dir.T());

            // the 8 neighbouring elements of element 'ind' can be found with ind+diff_mat[i]
            diff_mat = {-nely+1, 1, nely+1, nely, nely-1, -1, -nely-1, -nely};
            // dm_x[i] and dm_y[i] can be added to the element coordinates to get the 
            // coordinates of the neighbouring element corresponding to ind+diff_mat[i]
            dm_x = {-1, 0, 1, 1, 1, 0, -1, -1};
            dm_y = {1, 1, 1, 0, -1, -1, -1, 0};

            max_dep = 2; // in 2D, the arrival time of an element is calculated from the arrival time of 2 other elements
            Accepted.resize(nv,false); 
            T.resize(nv);
            dep_rho.resize(nv);
            prev.resize(nv*max_dep);
            dep_der.resize(nv*max_dep);
            order.resize(nv);
        }

        ~FPAMFilter() {
        }

        // apply AM filter to given input field xin
        void evaluate(double *xin, double *xout) {
            // initialize all elements
            for (int i=0; i<nv; i++) {
                Accepted[i] = false; // all nodes start as not accepted
                T[i] = std::numeric_limits<double>::max(); // set distance to infinite
            }

            // initialize arrival times at base plate. Numbering is top to bottom, left to right
            for (int i=nely-1; i<nv; i+=nely) {
				switch (speedFunctionType) {
					case 0:
						T[i] = -log(xin[i])/(kappa*log(2)); // set arrival time according to EQ.
						dep_rho[i] = -1.0/(kappa*log(2)*xin[i]); // derivative of T to xin
						break;
					case 1:
						T[i] = 1.0-xin[i]; // set arrival time according to EQ.
						dep_rho[i] = -1.0; // derivative of T to xin
						break;
				}
				
                
                pqueue.push(std::make_pair(T[i], i)); // add element to queue
                for (int j=0; j<max_dep; j++) prev[i*max_dep+j] = -1; // base plate elements do not depend on other nodes
            }

            // start front propagation loop 
            int c = 0;
            while (pqueue.size()>0) {
                int min_index = pqueue.top().second; // get element with lowest arrival time
                pqueue.pop(); // remove element from queue

                // if an element is already Accepted, skip it. As the arrival time of a node can be calculated from multiple other 
                // elements, its added to the queue multiple times because the STL priority_queue does not support updating entries 
                if (Accepted[min_index]) continue;

                Accepted[min_index] = true; // add to Accepted
                order[c] = min_index; // add to order.
                c++;

                int n_x = min_index / nely; // get x index of element
                int n_y = min_index % nely; // get y index of element
                
                // loop over 8 neighbours of min_index to find a target element to update
                for (int i=0; i<8; i++) {
                    int target = min_index + diff_mat[i];
                    // check if target is inside mesh, and not yet Accepted
                    if (n_x+dm_x[i]>=0 && n_x+dm_x[i]<nelx && n_y+dm_y[i]>=0 && n_y+dm_y[i]<nely && Accepted[target]==false) {
                        // find common neighbours of min_index and target
                        for (int j=0; j<8; j++) {
                            int neighbour2 = min_index + diff_mat[j];
                            // check if neighbour2 is next to target
                            if (std::any_of(diff_mat.begin(), diff_mat.end(), [target,neighbour2](int dm){return dm + neighbour2 == target;})) {
                                // check if neighbour2 is inside mesh, and Accepted
                                if (n_x+dm_x[j]>=0 && n_x+dm_x[j]<nelx && n_y+dm_y[j]>=0 && n_y+dm_y[j]<nely && Accepted[neighbour2]==true) {
                                    // update target from min_index and neighbour2
                                    if (updateArrivalTime(target, min_index, neighbour2, xin[target])) {
                                        // if the update resulted in a lower arrival time, add to queue
                                        pqueue.push(std::make_pair(T[target], target));
                                    }
                                }
                            }
                        }
                    }
                }
            }
            // post-processing: calculate densities from arrival times according to EQ
            for (int i=0; i<nv; i++) {
                double coordZ = (nely - 1) - (i % nely); // coordZ = 0 at base plate
				switch (speedFunctionType) {
					case 0:
						xout[i] = pow(2,-kappa*(T[i]-coordZ));
						break;
					case 1:
						xout[i] = 1.0-(T[i]-coordZ);
                        if (!isfinite(xout[i])) std::cout << "xout[i] not finite" << std::endl;
						break;
				}
            }
        }

        // calculate sensitivities
        void sens(double *sensIn, double *sensOut) {
            std::vector<double> lmbda(nv,0.0); // set up adjoint vector

            // loop over elements in reverse order as in which they were accepted
            for (int i=order.size()-1; i>=0; i--) {
                int index = order[i]; // get element index
                double coordZ = (nely - 1) - (index % nely); // calculate coordZ
				switch (speedFunctionType) {
					case 0:
						lmbda[index] += sensIn[index]*(-kappa*log(2)*pow(2,-kappa*(T[index] - coordZ))); // update adjoint
						break;
					case 1:
						lmbda[index] += -1.0*sensIn[index]; // update adjoint
						break;
				}
                sensOut[index] = dep_rho[index]*lmbda[index]; // calculate sensitivity value
                if (!isfinite(sensOut[index])) std::cout << "sensOut[index] not finite" << std::endl;
                // update adjoint of elements from which the arrival time of current element is calculated
                for (int j=0; j<max_dep; j++) {
                    if (prev[index*max_dep+j]>=0) {
                        lmbda[prev[index*max_dep+j]] += dep_der[index*max_dep+j]*lmbda[index];
                    }
                }
            }
        }

        // calculate the arrival time of target from source0 and source1 with known arrival time
        // This function is executed multiple times per elements, and is the computational bottleneck
        // of the algorithm. Make it efficient!
        bool updateArrivalTime(int target, int source0, int source1, double rho) {
            double p_min = 10.0; // factor for smooth minimum

            // create coordinate vectors for target and source elements
            // the y-coordinate is opposite to the ordering: it points bottom to top
            MatrixXX x1 {2,1, {double(source0 / nely), double((nely - 1) - source0 % nely)}};
            MatrixXX x2 {2,1, {double(source1 / nely), double((nely - 1) - source1 % nely)}};
            MatrixXX x3 {2,1, {double(target / nely), double((nely - 1) - target % nely)}};

            // vector 'gamma' determines the locations that are probed along the inverval 
            // from source0 to source1, fromt which target is updated. Normally these
            // locations can be determined exactly, but here we use a brute force for
            // brevity. 3 locations is sufficient in this case
            std::vector<double> gamma;
            int ns = 3;
            for (int i=0; i<=ns; i++) {
                gamma.push_back(i/ns);
            }

            bool updated = false;
            // calculate arrival time of target from locations determined by gamma
            for (int i=0; i<gamma.size(); i++) { 
                double alt_T, xiip, min_val, speed;

                MatrixXX x_dest = x1*gamma[i] + x2*(1-gamma[i]); // the location between x1 and x2 from which target is updated (EQ)
                double T_dest = T[source0]*gamma[i] + T[source1]*(1.0-gamma[i]); // interpolated arrival time at x_dest (EQ)
                MatrixXX d = x3 - x_dest; // vector from x_dest to target
                //double dr = std::max(std::abs(build_dir.dotVV(d)), C_p.dot(d).normL2() / vx); // speed in the direction d (EQ)
                double dr = std::max(std::abs(d.dat[0]), abs(d.dat[1]) / vx);

                // the speed scaling function g (EQ 10) for the old and improved implementation
                switch (speedFunctionType) {
                    case 0: // old version: linear scaling with denstiy according to EQ 11
                            speed = v_void + (1.0 - v_void)*rho;
                            alt_T = dr / speed + T_dest; // new arrival time
                            break;
                    case 1: // improved version: accordin to EQ 33 (slightly rewritten)
                            xiip = 1.0-(T_dest-x_dest.dat[1]);
                            min_val = (xiip-rho)/(1.0+std::exp(-p_min*(xiip-rho)));
			    speed = dr/(min_val+dr);
                            alt_T = (min_val+dr) + T_dest; // new arrival time
                            break;
                }     
               
                
                // if new arrival time is lower than old arrival time, update
                if (alt_T < T[target]) {
                    updated = true;
                    T[target] = alt_T;

                    // set dependency of target to source0 and source1
                    prev[target*max_dep+0] = source0;
                    prev[target*max_dep+1] = source1;

                    // calcualte partial derivatives according to speed function that was used
                    switch (speedFunctionType) {
                        case 0: {
                                    double dvdrho = -(1.0 - v_void) / (speed*speed);
                                    dep_der[target*max_dep+0] = gamma[i];
                                    dep_der[target*max_dep+1] = 1.0-gamma[i];
                                    dep_rho[target] = dr * dvdrho;
                                }
                                break;
                        case 1: 
								{
                                    double dvdspeed = 1.0;
									double dspeeddmin_val = 1.0/dr;
									double dmin_valdvals = 1.0/(std::exp(-p_min*(xiip-rho)) + 1.0) + (p_min*(xiip-rho)*std::exp(-p_min*(xiip-rho)))/(std::pow(std::exp(-p_min*(xiip-rho)) + 1.0,2.0));
									double dvalsdrho = -1.0;
									double dvalsdxiip = 1.0;
                                    double dxiipdu = -1.0;

                                    dep_der[target*max_dep+0] = gamma[i] + dr*dvdspeed*dspeeddmin_val*dmin_valdvals*dvalsdxiip*dxiipdu*gamma[i];
                                    dep_der[target*max_dep+1] = 1.0-gamma[i] + dr*dvdspeed*dspeeddmin_val*dmin_valdvals*dvalsdxiip*dxiipdu*(1.0-gamma[i]);
                                    dep_rho[target] = dr*dvdspeed*dspeeddmin_val*dmin_valdvals*dvalsdrho;
                                }
                                break;
                    }
                }
            }

            return updated;
        }
};

// define functions to interface with Python
extern "C" {
#ifdef _WIN32
#define EXPORTED __declspec(dllexport) 
#else
#define EXPORTED
#endif
    EXPORTED FPAMFilter* FPAMFilter_new(int _nelx, int _nely, int _speedFunctionType, double _kappa, double _v_void){ return new FPAMFilter(_nelx, _nely, _speedFunctionType, _kappa, _v_void); }
    EXPORTED void FPAMFilter_delete(FPAMFilter *fpamfilter) {delete fpamfilter;}
    EXPORTED void evaluate(FPAMFilter *fpamfilter, double *xin, double *xout){ fpamfilter->evaluate(xin, xout); }
    EXPORTED void sens(FPAMFilter *fpamfilter, double *sensIn, double *sensOut){ fpamfilter->sens(sensIn, sensOut); }

}

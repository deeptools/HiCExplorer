#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <set>
#include <map>
#include <math.h>
#include <algorithm>
#include <Eigen/Dense>
#include <limits>
namespace py = pybind11;

class kr_balancing{
     public:
      kr_balancing(std::vector<std::vector<double> > & input_matrix){
        A.resize(input_matrix.size(), input_matrix[0].size());
        for (int i = 0; i < input_matrix.size(); ++i){ //TODO how to avoid this for loop?!
            A.row(i) = Eigen::VectorXd::Map(&input_matrix[i][0], input_matrix[0].size());
        }
        e.resize(A.rows());
        e.setOnes();
        x0 = e;
      }
      void outter_loop(){
        double stop_tol = tol*0.5;
        double eta = etamax;
        x = x0;
        double rt = std::pow(tol,2);
        v = x.cwiseProduct(A*x);
        rk = Eigen::VectorXd::Constant(v.rows(),1) - v; //Vecotr of 1 - v
        rho_km1 = rk.conjugate().transpose()*rk;
        rho_km2 = rho_km1;
        double rout = rho_km1(0);
        double rold = rout;
        if(fl == 1) std::cout<< 'intermediate convergence statistics is off'<<std::endl;
        while(rout > rt){//outer itteration
          i=i+1; k=0; y=e;
          innertol = std::max(std::pow(eta,2)*rout,rt);
          inner_loop();
          x=x.cwiseProduct(y);
          v=x.cwiseProduct(A*x);
          rk = Eigen::VectorXd::Constant(v.rows(),1) - v;
          rho_km1 = rk.conjugate().transpose()*rk;
          rout = rho_km1(0);
          MVP = MVP + k + 1;
          //Update inner iteration stopping criterion.
          double rat = rout/rold; rold = rout; double res_norm = std::sqrt(rout);
          double eta_o = eta; eta = g*rat;
          if(g*std::pow(eta_o,2) > 0.1){
            eta = std::max(eta,g*std::pow(eta_o,2));
          }
          eta = std::max(std::min(eta,etamax),stop_tol/res_norm);
          if(fl == 1){
            res.push_back(res_norm);
          }
         }
         std::cout << "output"<<std::endl;
         Eigen::MatrixXd Ax = A.array().colwise() * x.array();
         xTAx = Ax.array().rowwise() * x.transpose().array();
      //   const static Eigen::IOFormat tsvFormat(-2, 1, "\t","\n");
      //   std::ofstream file("output.tsv");
      //   file << xTAx.format(tsvFormat);
      }
      void inner_loop(){
        while (rho_km1(0) > innertol){ // Inner itteration (conjugate gradient method)
          k++;
          if(k == 1){
            //Eigen::VectorXd Z = rk.cwiseProduct(v. cwiseInverse());
            Z = rk.cwiseQuotient(v);
            p=Z;
            rho_km1 = rk.conjugate().transpose()*Z;
          }else{
            Eigen::VectorXd beta=rho_km1.cwiseQuotient(rho_km2);
            p = Z + (beta(0)*p);

          }
          //Update search direction efficiently.
          Eigen::VectorXd w=x.cwiseProduct(A*(x.cwiseProduct(p))) + v.cwiseProduct(p);
          double alpha = rho_km1.cwiseQuotient(p.conjugate().transpose()*w)(0);
          Eigen::VectorXd ap = alpha * p;
          //Test distance to boundary of cone.
          Eigen::VectorXd ynew = y + ap;

          if(ynew.minCoeff() <= delta){
            if(delta == 0) break;
            Eigen::Matrix<bool, Eigen::Dynamic , Eigen::Dynamic> ind_helper = (ap.array()<0);
            Eigen::VectorXi ind = Eigen::VectorXi::LinSpaced(ind_helper.size(),0,ind_helper.size()-1);
            ind.conservativeResize(std::stable_partition(
              ind.data(), ind.data()+ind.size(), [&ind_helper](int i){return ind_helper(i);})-ind.data());
            Eigen::VectorXd y_copy = y;
            for(i = 0; i < ind.size(); i++){ //TODO proper masking? //ind_vec.unaryExpr(x);
               y_copy(ind(i)) = (delta - y_copy(ind(i)))/ap(i);
            }
            double gamma = y_copy.minCoeff();
            y = y + gamma*ap;
            break;
          }
          if(ynew.minCoeff() >= Delta){
            Eigen::Matrix<bool, Eigen::Dynamic , Eigen::Dynamic> ind_helper = (ynew.array() > Delta);
            Eigen::VectorXi ind = Eigen::VectorXi::LinSpaced(ind_helper.size(),0,ind_helper.size()-1);
            ind.conservativeResize(std::stable_partition(
              ind.data(), ind.data()+ind.size(), [&ind_helper](int i){return ind_helper(i);})-ind.data());
            Eigen::VectorXd y_copy = y;
            for(i = 0; i < ind.size(); i++){ //TODO proper masking? //ind_vec.unaryExpr(x);
                 y_copy(ind(i)) = (Delta - y_copy(ind(i)))/ap(i);
            }
            double gamma = y_copy.minCoeff();
            y = y + gamma*ap;
            break;
          }

          y = ynew;
          rk = rk - (alpha*w); rho_km2 = rho_km1;
          Z = rk.cwiseQuotient(v); rho_km1 = rk.conjugate().transpose()*Z;

        }//End of the inner 'while'
      }
      const Eigen::MatrixXd & get_output(){
        return xTAx;
      }
     private:
       std::vector<double> res;
       unsigned int fl = 0; //0 = on , 1 = off
       unsigned int Delta = 3;
       double delta = 0.1;
       double tol = 1e-6;
       double g = 0.9;
       double etamax = 0.1;
       Eigen::VectorXd x0;
       Eigen::VectorXd e;
       Eigen::MatrixXd A;
       Eigen::VectorXd rho_km1;
       Eigen::VectorXd rho_km2;
       unsigned int k;
       Eigen::VectorXd y;
       Eigen::VectorXd p;
       Eigen::VectorXd Z;
       double innertol;
       unsigned int i = 0; //Outer itteration count
       unsigned int MVP = 0;
       Eigen::VectorXd v;
       Eigen::VectorXd x;
       Eigen::VectorXd rk;
       //Output:
       Eigen::MatrixXd xTAx;
};


PYBIND11_MODULE(KRBalancing, m) {
  py::class_<kr_balancing>(m, "kr_balancing")
  //.def_buffer([](kr_balancing &m));
    .def(py::init<std::vector<std::vector<double> >& >())
    .def("outter_loop", &kr_balancing::outter_loop)
    .def("get_output",&kr_balancing::get_output, py::return_value_policy::reference_internal);
  //   return py::array(xTAx.size(), xTAx.data());
    //m.doc() = "pybind11 example plugin"; // optional module docstring

    //m.def("add", &add, "A function which adds two numbers");
}

//c++ -O3 -Wall -I . -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) KRBalancing.cpp -o KRBalancing$(python3-config --extension-suffix)

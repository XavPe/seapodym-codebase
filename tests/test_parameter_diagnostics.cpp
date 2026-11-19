#include <admodel.h>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include "ParameterDiagnostics.h"

int main()
{
    gradient_structure structure;
    int count=0;
    for(int hybrid:{0,1}) {
        gradient_structure::Hybrid_bounded_flag=hybrid;
        for(double x:{-.9,-.5,0.,.2,.9}) {
            for(double lower:{-20.,0.,5.}) {
                double upper=lower+30.;
                independent_variables v(1,1);v(1)=x;
                dvar_vector variables(v);
                dvariable penalty=0.;
                dvariable theta=boundp(variables(1),lower,upper,penalty);
                double physical=value(theta);
                dvector derivative(1,1);gradcalc(1,derivative);
                double jac=seapodym_diagnostics::parameter_jacobian(x,lower,upper);
                if(std::fabs(jac-derivative(1))>1e-9*std::max(1.,std::fabs(jac))) {
                    std::cerr<<"Jacobian mismatch hybrid="<<hybrid<<" x="<<x<<" calculated="<<jac<<" AD="<<derivative(1)<<"\n";
                    return 1;
                }
                if(physical!=0.) {
                    double cv=seapodym_diagnostics::physical_cv(physical,.5,jac);
                    assert(std::fabs(cv-std::fabs(derivative(1))*.5/std::fabs(physical))<1e-9);
                }
                ++count;
            }
        }
    }
    gradient_structure::Hybrid_bounded_flag=0;
    const double pi=std::acos(-1.);
    double cv=seapodym_diagnostics::physical_cv(10.,.5,seapodym_diagnostics::parameter_jacobian(0.,0.,20.));
    assert(std::fabs(cv-pi/4)<1e-9);
    assert(std::isnan(seapodym_diagnostics::physical_cv(0.,1.,1.)));
    assert(std::isnan(seapodym_diagnostics::physical_cv(1.,std::numeric_limits<double>::quiet_NaN(),1.)));
    std::cout<<"jacobian_AD_cases="<<count<<" physical_cv="<<cv<<" zero_estimate_and_invalid_SE=NaN\n";
}

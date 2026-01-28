#include "roots.hpp"
#include <cmath>
const double tolerance = 1e-6;

bool bisection(std::function<double(double)> f,
               double a, double b,
               double *root)
{    
        
        
        if (f(a) * f(b) > 0) {
            return false;
        }
            
        for (int i = 0; i < 1e6; ++i) {
            double mid = (a+b)/2;
            
            
            if (abs(f(mid)) < tolerance){
                *root = mid;
                return true;
            }
                
            if (f(a) * f(mid) < 0){
                b = mid;
                
            } else{
                a = mid;
                
            }

            }

            *root = (a + b)/2;
            return true;
        }

bool regula_falsi(std::function<double(double)> f,
                  double a, double b,
                  double *root)
{
    
    if (f(a) * f(b) > 0){
        return false;
    }
    for ( int i = 0;  i < 1e6; i++){
        double mid = a - (f(a)*(b - a))/(f(b) - f(a));
        
        if (abs(f(mid)) < tolerance){
            *root = mid;
        }

        if (f(a) * f(mid) < 0){
                b = mid;
                
            } else{
                a = mid;
                
            } 

    }
    *root = a - (f(a)*(b - a))/(f(b) - f(a));
    return true;
}

bool newton_raphson(std::function<double(double)> f,
                    std::function<double(double)> g,
                    double a, double b, double c,
                    double *root)
{
    
    for ( int i = 0;  i < 1e6; i++){
        c -= f(c)/g(c);
        
        if (abs(f(c)/g(c)) < tolerance){
            *root = c;
            return true;
        }
    }

    return false;
}

bool secant(std::function<double(double)> f,
            double a, double b, double c,
            double *root)
{
    double d = b;
    for (int i = 0;  i < 1e6; i++){
        double denom = f(d)-f(c);
    
        double dn = d-f(d)*(d-c)/denom;
        c = d;
        d = dn;

        if (abs(d-c) < tolerance){
            *root = d;
            return true;
        }


    }
    return false;
}



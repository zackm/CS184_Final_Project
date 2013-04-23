#include "FluidSimulation.h"


double Wpoly6(double r_square) 
{
    static double coefficient = 315.0/(64.0*PI*pow((double)H, 9));
    static double h_square = H*H;

    return coefficient*pow(h_square-r_square, 3);
}

//double Wpoly6(double *r)
//{
//    static double coefficient = 315.0/(64.0*PI*pow((double)H, 9));
//    static double h_square = H*H;
//
//    double r_square = DotProduct(r,r);
//
//    return coefficient*pow(h_square-r_square, 3);
//}

void Wpoly6Grad(double *r, double r_square, double *gradient)
{
    static double coefficient = -6.0 * 315.0/(64.0*PI*pow((double)H, 9));
    static double h_square = H*H;

    double weight = coefficient*pow(h_square-r_square, 2);

    for(int i = 0; i < 3; i++) gradient[i] = weight * r[i];
}


double Wpoly6Laplacian(double r_square)
{
    static double coefficient = 315.0/(64.0*PI*pow((double)H, 9));
    static double h_square = H*H;

    double h_square_minus_r_square = h_square - r_square;

    return coefficient*(-18.0*pow(h_square_minus_r_square, 2) + 24.0*r_square*h_square_minus_r_square);
}

//double Wpoly6Laplacian(double *r)
//{
//    static double coefficient = 315.0/(64.0*PI*pow((double)H, 9));
//    static double h_square = H*H;
//
//    double r_square = DotProduct(r,r);
//    double h_square_minus_r_square = h_square - r_square;
//
//    return coefficient*(-18.0*pow(h_square_minus_r_square, 2) + 24.0*r_square*h_square_minus_r_square);
//}


double Wspiky(double *r)
{
    static double coefficient = 15.0/(PI*pow((double)H, 6));

    double r_norm = sqrt(DotProduct(r,r));

    return coefficient*pow(H-r_norm, 3);
}



void WspikyGrad(double *r, double r_square, double *gradient) // from WspikyGrad
{
    static double coefficient = -3.0 * 15.0/(PI*pow((double)H, 6));

    double r_norm = sqrt(r_square);
    double h_minus_r = H - r_norm;

    double weight = coefficient*pow(h_minus_r, 2) / r_norm;

    for(int i = 0; i < 3; i++) gradient[i] = weight * r[i];
}


//void WspikyGrad(double *r, double *gradient)
//{
//    static double coefficient = -3.0 * 15.0/(PI*pow((double)H, 6));
//
//    double r_norm = sqrt(DotProduct(r,r));
//    double h_minus_r = H - r_norm;
//
//    double weight = coefficient*pow(h_minus_r, 2) / r_norm;
//
//    for(int i = 0; i < 3; i++) gradient[i] = weight * r[i];
//}

double WspikyLaplacian(double *r)
{
    static double coefficient = 15.0/(PI*pow((double)H, 6));

    double r_norm = sqrt(DotProduct(r,r));
    double h_minus_r = H - r_norm;

    return coefficient*6.0*(h_minus_r - pow(h_minus_r, 2) / r_norm);
}


double Wviscosity(double *r)
{
    static double coefficient  = 15.0/(2.0*PI*pow((double)H, 3));
    static double inv_h_cubic  = 1.0/pow((double)H, 3);
    static double inv_h_square = 1.0/pow((double)H, 2);

    double r_square = DotProduct(r,r);
    double r_norm   = sqrt(r_square);
    double r_cubic  = r_square * r_norm;

    return coefficient*(-r_cubic*0.5*inv_h_cubic + r_square*inv_h_square + H*0.5/r_norm - 1.0);
}

void WviscosityGrad(double *r, double *gradient)
{
    static double coefficient  = 15.0/(2.0*PI*pow((double)H, 3));
    static double inv_h_cubic  = 1.0/pow((double)H, 3);
    static double inv_h_square = 1.0/pow((double)H, 2);

    double r_square = DotProduct(r,r);
    double r_norm   = sqrt(r_square);
    double r_cubic  = r_square * r_norm;

    double weight = coefficient*(-3.0*r_norm*0.5*inv_h_cubic + 2.0*inv_h_square - H*0.5/r_cubic);

    for(int i = 0; i < 3; i++) gradient[i] = weight * r[i];
}

double WviscosityLaplacian(double r_square) // from WviscosityLaplacian
{
    static double coefficient  = 45.0/(PI*pow((double)H, 6));

    double r_norm = sqrt(r_square);

    return coefficient*(H - r_norm);
}

//double WviscosityLaplacian(double *r)
//{
//    static double coefficient  = 45.0/(PI*pow((double)H, 6));
//
//    double r_square = DotProduct(r,r);
//    double r_norm   = sqrt(r_square);
//
//    return coefficient*(H - r_norm);
//}
#include <iostream>
#include <armadillo>
#include <fstream>

using namespace arma;
using namespace std;

int main()
{
    double m,hbar,dx,dt,omega2,omega,x0;
    complex<double> alph;
    complex<double> i (0,1.0);
    int n = 30;
    int nt = 4000;
    int nit = 100;
    cx_mat A = cx_mat(n,n);
    A.fill(0);
    cx_mat Anew;
    cx_mat Aold;
    mat V(n,n);
    ofstream outfile;
    outfile.open("matrix.dat");
    vec xx = linspace(0,1,n);
    vec yy = linspace(0,1,n);
    x0=0.5;
    hbar = 1;
    m    = 1;
    dx   = 1.0/(n-1);
    dt   = 0.0001;
    omega = 30.0;
    omega2 = omega*omega; //omega squared
    alph = i*hbar*dt/(2*m*dx*dx);

    for (int i = 1; i < n-1; i++) {
        for (int j = 1; j < n-1; j++) {
             A(i,j) = (4*yy(j)*yy(j) -2)*exp( -(((xx(i) - x0)*(xx(i) - x0) + (yy(j) - x0)*(yy(j) - x0)))*m*omega/(2*hbar)) + 2*xx(i)*exp( -(((xx(i) - x0)*(xx(i) - x0) + (yy(j) - x0)*(yy(j) - x0)))*m*omega/(2*hbar));
           //1/M_PI*sin(M_PI*xx(i))*sin(2*M_PI*yy(j)) + 1/M_PI*sin(M_PI*xx(i))*sin(M_PI*yy(j)) + 1/M_PI*sin(3*M_PI*xx(i))*sin(5*M_PI*yy(j));        //exp( -(((xx(i) - m)*(xx(i) - m) + (yy(j) - m)*(yy(j) - m))) / (0.01));
            V(i,j) = 2*m*dx*dx/(hbar*hbar)*(1/2.0*omega2*(xx(i)-x0)*(xx(i)-x0) +1/2.0*omega2*(yy(j)-x0)*(yy(j)-x0));
        }
    }


    outfile << abs(A) << endl;

    Anew = A;
    Aold = A;


    // Print info.
    printf(" ------------------------------\n");
    printf("| n   = %22d | \n", n);
    printf("| nt  = %22d | \n", nt);
    printf("| nit = %22d | \n", nit);
    printf(" ------------------------------\n");
    printf("| dx  = %22.4f | \n", dx);
    printf("| dt  = %22.4f | \n", dt);
    printf(" ------------------------------\n");

    for(int k = 0;k<nt;k++) {//

        // Print progress.
        if (k % 5 == 0) {
            printf("Progress: %8.2f \% \r", 100*k/ (double)nt);
        }

        for(int l = 0;l<nit;l++){
            for(int i = 1;i<n-1;i++) {
                for(int j = 1;j<n-1;j++) {
                    Anew(i,j) = Aold(i,j) + alph*(A(i,j+1) + A(i-1,j) + A(i+1,j) + A(i,j-1));
                    Anew(i,j) = Anew(i,j)/(1.0 +(4.0+V(i,j))*alph);
                }
            }
            A = Anew;
        }
        Aold = A;
        outfile << abs(A) << endl;
    }
    outfile.close();
    return 0;
}


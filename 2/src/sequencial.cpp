#include "sequencial.h"

/*
   function [u,ITER]=PoissonGS(N, TOL)
   % Finds the steady-state solution for the temperature distribution on a
   % square plate, for particular boundary conditions
   % USES GAUSS-SEIDEL ITERATIONS
   % N is the number of grid points in each direction (including boundary
   % points)
   % TOL is the tolerance for the stopping citeria
   % IN THE GRID, USES i FOR ROW (ASCENDENT), j FOR COLUMN (LEFT-RIGHT)
   % SEE Quinn's book, pp.330-332

   % SET BOUNDARY VALUES
   % Temperature is zero at top and 100 on the other boundaries
   u(1,1:N)=100;  % lower boundary
   u(1:N,1)=100;  % left boundary
   u(1:N,N)=100;  % right boundary
   u(N,1:N)=0;    % boundary above

   w(1,1:N)=100;  % lower boundary
   w(1:N,1)=100;  % left boundary
   w(1:N,N)=100;  % right boundary
   w(N,1:N)=0;    % boundary above

   % initial values for interior points
   u(2:N-1,2:N-1)=50;

   % COMPUTE STEADY STATE SOLUTION
   DIFF=TOL+1;
   ITER=0;
   while DIFF>TOL
   for i=2:N-1
   for j=2:N-1
   w(i,j)=(w(i-1,j)+w(i,j-1)+u(i,j+1)+u(i+1,j))/4;
   end
   end
   DIFF=max(max(abs(w-u)));
   u=w;
   ITER=ITER+1;
   end

   % DISPLAYS A COLOURFUL MAP OF DISTRIBUTIONS
   u=flipud(u); image(u)
   */

void poisson_gs(const int n, int tol) {
    int diff = tol + 1;
    int iter = 0;

    Matrice<int> u(n, n);
    Matrice<int> w(n, n);

    for (int i = 0; i < n; i++) {
        u(0, i) = 100;
        w(0, i) = 100;
        u(i, 0) = 100;
        w(i, 0) = 100;
        u(i, n - 1) = 100;
        w(i, n - 1) = 100;
        u(n - 1, i) = 0;
        w(n - 1, i) = 0;
    }

    for (int i = 1; i < n - 1; i++) {
        for (int j = 1; j < n - 1; j++) {
            u(i, j) = 50;
        }
    }

    while (diff > tol) {
        for (int i = 1; i < n - 1; i++)
            for (int j = 1; j < n - 1; j++)
                w(i, j) = (w(i - 1, j) + w(i, j - 1) + u(i, j + 1) + u(i + 1, j)) / 4;

        diff = (w - u).max();
        u = w;
        iter++;
    }

    std::cout << u << std::endl;
    std::cout << w << std::endl;
}

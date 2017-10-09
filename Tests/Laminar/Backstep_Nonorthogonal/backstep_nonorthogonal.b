%-- phisical properties
  1
  1  Fluid  0.001  1.00  
%-- boundary conditions 
   4 
   1    WALL      0.0      0.0   0.0
   2    WALL      0.0      0.0   0.0
   3    INFLOW    file profile.dat      
   4    OUTFLOW   0.5      0.0   0.0
%-- initial conditions  
   1
   1   1.00  0.0   0.0

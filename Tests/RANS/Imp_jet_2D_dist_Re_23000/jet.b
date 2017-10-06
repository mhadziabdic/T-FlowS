#::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                      #
# Physical properties, Boundary and Initial conditions #
#                                                      #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::#


#------------------------#
#-  physical properties -#
#------------------------#
#  No.   Type         VISc     DENc     CONs      CAPc
   1
   1     Fluid       1.00E-4   1.00    1.40e-4    1.0 

#------------------------#
#-- boundary condition --#
#------------------------#
#  No.   Type          U       V       W       T      Kin     Eps     v_2     f22
   6
#--lower wall
   1     WallFlux   0.0    0.0    0.0  0.1     0.0     1.0e-3  0.0     0.0   
#--pipe wall
   2     WallFlux   0.0    0.0    0.0  0.0     0.0     1.0e-3  0.0     0.0
   3     PRESSURE   0.0      0.0     0.0    0.0  20.0    1.0E-2  1.0E-3  6.6E-2  1.0e-3
   4     INFLOW       File    InletProfile_zeta_Re23000.dat 
   5     SYMMETRY        0.0     0.0    0.0  20.0    1.0E-2  1.0E-3  6.6E-2  1.0e-3

#--outer outlet
   6     PRESSURE   0.0      0.0     0.0    0.0  20.0    1.0E-2  1.0E-3  6.6E-2  1.0e-3

#--sides

#------------------------#
#-- initial conditions --#
#------------------------#
#  No.                  U       V       W       T      Kin     Eps     v_2     f22
   1
   1                   -0.1     0.0     0.0     20.0   1.0E-3  1.0E-4  6.6E-4  1.0e-4

#
#
#     7---------------8------------------16
#    /|              /|                  /|
#   5---------------6------------------15 |
#   | |            /| |                 | |
#   | 3 - - - - -4/ |                   | |
#   |/           /| | |                 | |
#   1-----------2 | |                   | |
#               | | | |                 | |
#               |11-|14-----------------|12
#               |/  |/                  |/
#               9--13------------------10
#
#-------------------------------------------#
#  Nodes (cells), boundary cells and sides  #
#-------------------------------------------#
  30000 10000 90000

#----------#
#  Points  #
#----------#
16
  1   0.0             0.0             0.0
  2   2.0             0.0             0.0
  3   0.0             1.0             0.0
  4   2.0             1.0             0.0
  5   0.0             0.0             1.0            
  6   3.0             0.0             1.0            
  7   0.0             1.0             1.0            
  8   3.0             1.0             1.0            
  9   2.0             0.0            -1.0
 10   8.0             0.0            -1.0
 11   2.0             1.0            -1.0
 12   8.0             1.0            -1.0
 13   3.0             0.0            -1.0
 14   3.0             1.0            -1.0
 15   8.0             0.0             1.0
 16   8.0             1.0             1.0
#----------#
#  Blocks  #
#----------#
3
  1   21   4  11 
      1.0  1.0  1.0 
      1  2  3  4  5  6  7  8
  2   11   4  11 
      1.0  1.0  1.0 
      2  9  4 11  6 13  8 14
  3   51   4  11 
      1.0  1.0  1.0 
     13 10 14 12  6 15  8 16
#--------#
#        #  
#--------#
   0 
   0
#-----------------------#
#  Boundary conditions  #
# (it will use default) #
#-----------------------#
  2
    1     Imin
        1   2
    2     Imax
        3   3
#-----------------------#
#  Periodic boundaries  #
#-----------------------#
  3
      1    1  2  6  5
           3  4  8  7
      2    2  9 13  6
           4 11 14  8
      3   13 10 15  6
          14 12 16  8
#-------------------#
#  Copy boundaries  #
#-------------------#
   0
#------------
# Refinement
#-----------
   0
#------------
# Smoothing 
#-----------
  1
    1  x  y  z
       10  0.5
      -0.50 -0.50 -1.50  9.00  1.50  1.50

$//////////////////////////////////////////////////
$
$  18-------------19-------------20
$   |\            /|\             |\
$   | \              \              \
$   |  \        /  |  \           |  \ 
$   |  15-------------16-------------17
$   | 13|\---14    |  /|          |   |
$   |   ||\   |\     / |              | 
$   |   || \  | \  |/  |          |   |
$   |   |9- \10  \ /   |              |
$   |   | \ 11---12|   |          |   |
$   |  /|  \ |  \ |    |              |
$   |   |   \|   \||   |          |   |
$   |/  |    7----8    |              |
$   4 - | - / - - -\- -|- - - - - 6   |
$    \  |  /       5\  |           \  | 
$     \ | /          \ |            \ |
$      \|/            \|             \|
$       1--------------2--------------3
$
$//////////////////////////////////////////////////

%===========================================%
%  Nodes (cells), boundary cells and sides  %
%===========================================%
  300000 100000 900000

%===========%
%   Nodes   %
%===========%
20
%-----------%
%   Floor   %
%-----------%
 1   0.0   0.0   0.0
 2   3.0   0.0   0.0
 3   9.0   0.0   0.0
 4   0.0   3.0   0.0
 5   3.0   3.0   0.0
 6   9.0   3.0   0.0
%--------------%
%   Cylinder   %
%--------------%
 7   1.25   0.0   1.25
 8   1.75   0.0   1.25
 9   1.25   3.0   1.25
10   1.75   3.0   1.25
11   1.25   0.0   1.75
12   1.75   0.0   1.75
13   1.25   3.0   1.75
14   1.75   3.0   1.75
%----------%
%   Roof   %
%----------%
15   0.0   0.0   3.0
16   3.0   0.0   3.0
17   9.0   0.0   3.0
18   0.0   3.0   3.0
19   3.0   3.0   3.0
20   9.0   3.0   3.0
%============%
%   Blocks   %
%============%
  5 
$================
$
$   +-----+-----+
$   |\ 3 /|     |
$   | +-+ |     |
$   |1| |4|  5  |
$   | +-+ |     |
$   |/ 2 \|     |
$   +-----+-----+
$    
$================
1   21 21 21   
    4.0  1.0  1.0
    1  7  4  9 15 11 18 13 
2   21 21 21   
    1.0  1.0  4.0
    1  2  4  5  7  8  9 10
3   21 21 21   
    1.0  1.0  0.25
   11 12 13 14 15 16 18 19
4   21 21 21   
    0.25 1.0  1.0
    8  2 10  5 12 16 14 19
5   21 21 21   
    0.25 1.0  1.0
    2  3  5  6 16 17 19 20
%%%%%%%%%%%%%%%%%%%%
%   line section   %
%%%%%%%%%%%%%%%%%%%%
  8
*---- front
      1    8  12
   71    1.853554     .000000    1.146447
   72    1.880203     .000000    1.175276
   73    1.904509     .000000    1.206107
   74    1.926320     .000000    1.238751
   75    1.945503     .000000    1.273005
   76    1.961940     .000000    1.308658
   77    1.975528     .000000    1.345492
   78    1.986185     .000000    1.383277
   79    1.993844     .000000    1.421783
   80    1.998459     .000000    1.460770
    1    2.000000     .000000    1.500000
    2    1.998459     .000000    1.539230
    3    1.993844     .000000    1.578217
    4    1.986185     .000000    1.616723
    5    1.975528     .000000    1.654508
    6    1.961940     .000000    1.691342
    7    1.945503     .000000    1.726995
    8    1.926320     .000000    1.761249
    9    1.904508     .000000    1.793893
   10    1.880203     .000000    1.824724
   11    1.853553     .000000    1.853553
      2   12  11
   11    1.853553     .000000    1.853553
   12    1.824724     .000000    1.880203
   13    1.793893     .000000    1.904508
   14    1.761249     .000000    1.926320
   15    1.726995     .000000    1.945503
   16    1.691342     .000000    1.961940
   17    1.654508     .000000    1.975528
   18    1.616723     .000000    1.986185
   19    1.578217     .000000    1.993844
   20    1.539230     .000000    1.998459
   21    1.500000     .000000    2.000000
   22    1.460770     .000000    1.998459
   23    1.421783     .000000    1.993844
   24    1.383277     .000000    1.986185
   25    1.345492     .000000    1.975528
   26    1.308658     .000000    1.961940
   27    1.273005     .000000    1.945503
   28    1.238751     .000000    1.926320
   29    1.206107     .000000    1.904508
   30    1.175276     .000000    1.880203
   31    1.146447     .000000    1.853553
      3   11   7
   31    1.146447     .000000    1.853553
   32    1.119797     .000000    1.824724
   33    1.095492     .000000    1.793893
   34    1.073680     .000000    1.761249
   35    1.054497     .000000    1.726995
   36    1.038060     .000000    1.691342
   37    1.024472     .000000    1.654508
   38    1.013815     .000000    1.616723
   39    1.006156     .000000    1.578217
   40    1.001541     .000000    1.539230
   41    1.000000     .000000    1.500000
   42    1.001541     .000000    1.460770
   43    1.006156     .000000    1.421783
   44    1.013815     .000000    1.383277
   45    1.024472     .000000    1.345492
   46    1.038060     .000000    1.308658
   47    1.054497     .000000    1.273005
   48    1.073680     .000000    1.238751
   49    1.095492     .000000    1.206107
   50    1.119797     .000000    1.175276
   51    1.146447     .000000    1.146447
      4    7   8
   51    1.146447     .000000    1.146447
   52    1.175276     .000000    1.119797
   53    1.206107     .000000    1.095491
   54    1.238751     .000000    1.073680
   55    1.273005     .000000    1.054497
   56    1.308658     .000000    1.038060
   57    1.345491     .000000    1.024472
   58    1.383277     .000000    1.013815
   59    1.421783     .000000    1.006156
   60    1.460770     .000000    1.001541
   61    1.500000     .000000    1.000000
   62    1.539230     .000000    1.001541
   63    1.578217     .000000    1.006156
   64    1.616723     .000000    1.013815
   65    1.654509     .000000    1.024472
   66    1.691342     .000000    1.038060
   67    1.726995     .000000    1.054497
   68    1.761249     .000000    1.073680
   69    1.793893     .000000    1.095492
   70    1.824724     .000000    1.119797
   71    1.853554     .000000    1.146447
*---- back
      5   10  14
   71    1.853554    3.000000    1.146447
   72    1.880203    3.000000    1.175276
   73    1.904509    3.000000    1.206107
   74    1.926320    3.000000    1.238751
   75    1.945503    3.000000    1.273005
   76    1.961940    3.000000    1.308658
   77    1.975528    3.000000    1.345492
   78    1.986185    3.000000    1.383277
   79    1.993844    3.000000    1.421783
   80    1.998459    3.000000    1.460770
    1    2.000000    3.000000    1.500000
    2    1.998459    3.000000    1.539230
    3    1.993844    3.000000    1.578217
    4    1.986185    3.000000    1.616723
    5    1.975528    3.000000    1.654508
    6    1.961940    3.000000    1.691342
    7    1.945503    3.000000    1.726995
    8    1.926320    3.000000    1.761249
    9    1.904508    3.000000    1.793893
   10    1.880203    3.000000    1.824724
   11    1.853553    3.000000    1.853553
      6   14  13
   11    1.853553    3.000000    1.853553
   12    1.824724    3.000000    1.880203
   13    1.793893    3.000000    1.904508
   14    1.761249    3.000000    1.926320
   15    1.726995    3.000000    1.945503
   16    1.691342    3.000000    1.961940
   17    1.654508    3.000000    1.975528
   18    1.616723    3.000000    1.986185
   19    1.578217    3.000000    1.993844
   20    1.539230    3.000000    1.998459
   21    1.500000    3.000000    2.000000
   22    1.460770    3.000000    1.998459
   23    1.421783    3.000000    1.993844
   24    1.383277    3.000000    1.986185
   25    1.345492    3.000000    1.975528
   26    1.308658    3.000000    1.961940
   27    1.273005    3.000000    1.945503
   28    1.238751    3.000000    1.926320
   29    1.206107    3.000000    1.904508
   30    1.175276    3.000000    1.880203
   31    1.146447    3.000000    1.853553
      7   13   9
   31    1.146447    3.000000    1.853553
   32    1.119797    3.000000    1.824724
   33    1.095492    3.000000    1.793893
   34    1.073680    3.000000    1.761249
   35    1.054497    3.000000    1.726995
   36    1.038060    3.000000    1.691342
   37    1.024472    3.000000    1.654508
   38    1.013815    3.000000    1.616723
   39    1.006156    3.000000    1.578217
   40    1.001541    3.000000    1.539230
   41    1.000000    3.000000    1.500000
   42    1.001541    3.000000    1.460770
   43    1.006156    3.000000    1.421783
   44    1.013815    3.000000    1.383277
   45    1.024472    3.000000    1.345492
   46    1.038060    3.000000    1.308658
   47    1.054497    3.000000    1.273005
   48    1.073680    3.000000    1.238751
   49    1.095492    3.000000    1.206107
   50    1.119797    3.000000    1.175276
   51    1.146447    3.000000    1.146447
      8    9  10
   51    1.146447    3.000000    1.146447
   52    1.175276    3.000000    1.119797
   53    1.206107    3.000000    1.095491
   54    1.238751    3.000000    1.073680
   55    1.273005    3.000000    1.054497
   56    1.308658    3.000000    1.038060
   57    1.345491    3.000000    1.024472
   58    1.383277    3.000000    1.013815
   59    1.421783    3.000000    1.006156
   60    1.460770    3.000000    1.001541
   61    1.500000    3.000000    1.000000
   62    1.539230    3.000000    1.001541
   63    1.578217    3.000000    1.006156
   64    1.616723    3.000000    1.013815
   65    1.654509    3.000000    1.024472
   66    1.691342    3.000000    1.038060
   67    1.726995    3.000000    1.054497
   68    1.761249    3.000000    1.073680
   69    1.793893    3.000000    1.095492
   70    1.824724    3.000000    1.119797
   71    1.853554    3.000000    1.146447
%%%%%%%%%%%%%%%%%%%%%%%
%   surface section   %  -> still under construction
%%%%%%%%%%%%%%%%%%%%%%%
  0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   boundary conditions and material data   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  2
#--------#
# Inflow #
#--------#
    1     Imin 
        1   2
#---------#
# Outflow #
#---------#
    2     Imax 
        5   3
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   periodic boundaries   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
  5
  1    1  2  8  7  
       4  5 10  9
  2    2 16 12  8
       5 19 14 10
  3   11 12 16 15
      13 14 19 18
  4    1  7 11 15
       4  9 13 18
  5    2  3 17 16
       5  6 20 19
%%%%%%%%%%%%%%%%%%%%%%
%   copy boundaries  %
%%%%%%%%%%%%%%%%%%%%%%
  0
%%%%%%%%%%%%%%%%%%
%   refinement   %
%%%%%%%%%%%%%%%%%%
  1
     1  1
        1   FILL ELIPSOID
            2.3  1.5  1.5  1.5  1000.0  1.0
%%%%%%%%%%%%%%%%%%%%
%   don't smooth   %
%%%%%%%%%%%%%%%%%%%%
  0

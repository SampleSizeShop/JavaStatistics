data yUnivariate;
input y x;
datalines;
2.67 1
0.97 1
0.55 1
1.57 1
-0.12 1
-1.02 1
0.94 1
-1.25 1
0.94 1
-0.08 1
2.09 0
2.7800000000000002 0
1.7 0
1.4100000000000001 0
1.47 0
2.2 0
4.890000000000001 0
0.96 0
4.05 0
0.9199999999999999 0
;
run;

proc glm data=yUnivariate;
title "univariate";
model y = x ;
run;

data yMulti;
input y1 y2 y3 x1 x2 x3 x4 x;
datalines;
0.43999999999999995 -0.88 0.58 1 0 0 0 1
4.029999999999999 -1.88 -0.16 1 0 0 0 1
2.19 -0.07 0.67 1 0 0 0 1
2.83 -0.22 0.63 1 0 0 0 1
0.6000000000000001 -1.47 1.71 1 0 0 0 1
0.07 -0.53 -0.09 0 1 0 0 2
-0.2 -0.87 0.06 0 1 0 0 2
0.27 -1.12 -1.42 0 1 0 0 2
-0.01 0.08 -0.38 0 1 0 0 2
1.66 -0.53 -0.75 0 1 0 0 2
0.25 1.73 0.08 0 0 1 0 3
-1.05 0.61 -0.15 0 0 1 0 3
2.01 -0.54 -1.34 0 0 1 0 3
-1.06 0.16 1.88 0 0 1 0 3
-1.43 2.18 -0.39 0 0 1 0 3
-0.75 -0.19 -1.46 0 0 0 1 4
-1.08 -0.91 0.8 0 0 0 1 4
-1.16 -1.68 0.61 0 0 0 1 4
1.67 1.0 0.55 0 0 0 1 4
-0.47 -0.07 0.14 0 0 0 1 4
;
run;


data yMulti2;
input y1 y2 y3 x1 x2 x3 x4 x;
datalines;
1 2 3 1 0 0 0 1
1 2 3 1 0 0 0 1
1 2.7 3 1 0 0 0 1
1 1 3 1 0 0 0 1
1 2 3 1 0 0 0 1
1 4 2 0 1 0 0 2
1 2 2 0 1 0 0 2
1 3.5 2 0 1 0 0 2
1 4 2 0 1 0 0 2
1 4 2 0 1 0 0 2
1 4 2 0 0 1 0 3
1 6 2 0 0 1 0 3
1 4 2 0 0 1 0 3
1 5 2 0 0 1 0 3
1 4 2 0 0 1 0 3
-1 0 2 0 0 0 1 4
-1 0 1 0 0 0 1 4
-1 0 2 0 0 0 1 4
-1.2 0 1.5 0 0 0 1 4
-1 0 2 0 0 0 1 4
;
run;


proc glm data=yMulti;
title "multivariate 1";
class x;
model y1 y2 y3 = x / noint solution;
*model y1 y2 y3 = x1 x2 x3 x4 / noint solution;
contrast 'x pairwise effects' x 1 -1 0 0, x 1 0 -1 0, x 1 0 0 -1 / E;
*contrast 'x pairwise effects' x1 x2 x3 x4 1 -1 0 0, x1 x2 x3 x4 1 0
-1 0, x1 x2 x3 x4 1 0 0 -1;
manova h = x / printe printh ;
run;

/*
proc glm data=yMulti2;
title "multivariate 2";
class x;
model y1 y2 y3 = x / noint solution;
contrast 'x pairwise effects' x 1 -1 0 0, x 1 0 -1 0, x 1 0 0 -1;
manova h = x / printe printh ;
run;*/

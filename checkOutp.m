
%% this is for a case where you use a 2 x 2 grid
%%
p11 = load('../OUT/out_0');
p12 = load('../OUT/out_1');
p21 = load('../OUT/out_2');
p22 = load('../OUT/out_3');

%%-------------------- gather 
pall = [p11, p12; p21 p22];

%%-------------------- checking
t1 = sum(sum(pall))
%%-------------------- plot
surf(pall) 


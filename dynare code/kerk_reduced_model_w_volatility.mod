% A DSGE Model w 6 period gov. debt and no rigidities.
%
% Aaron Betz
% Aug. 27 2018

%----------------------------------------------------------------
% 0. Housekeeping (close all graphic windows, import fiscal params)
%----------------------------------------------------------------


close all;

%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------


var z sigma y c k i l w rk q q2 q3 q4 q5 q6 Bm1 Bm2 Bm3 Bm4 Bm5 B1 B2 B3 B4 B5 B6 uc; 

varexo uep usig;

parameters rhoz rhosig beta psi phi delta alpha gam;

%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------

gam=2.55;
rhoz=0.9;
beta    = 0.9975;
alpha   = 0.36;
delta   = 0.015;
rhosig     = 0.9;  
psi=1;
phi=0.55;


model; 
//1
z=rhoz*z(-1)+sigma*uep;

//2-7 Bond lag equations and portfolio (Also GBC, b1=0.0)
B6=B1;
B5=B1;
B4=B1;
B3=B1;
B2=B1;
B1=0.0;
Bm5=B6(-1)+B5;
Bm4=Bm5(-1)+B4;
Bm3=Bm4(-1)+B3;
Bm2=Bm3(-1)+B2;
Bm1=Bm2(-1)+B1;

//4 LoM for capital
i=k-(1-delta)*k(-1);

//5 Production function
y=k(-1)^alpha*(exp(z)*l)^(1-alpha);

//6 MPL
w=(1-alpha)*y/l;

//7 MPK
rk=alpha*(y/k(-1));

//8 marg. util. consum.
uc=phi*(c^gam+psi*(1-l)^gam)^((phi-gam)/gam)*c^(gam-1);

//9 labor f.o.c.
c^(gam-1)*w=psi*(1-l)^(gam-1);

//10 Capital Euler
uc=beta*uc(+1)*(1+rk(+1)-delta);

//11 Budget Constraint
c=w*l+rk*k(-1)-i+Bm1(-1)+q*Bm2(-1)+q2*Bm3(-1)+q3*Bm4(-1)+q4*Bm5(-1)+q5*B6(-1)-q*B1-q2*B2-q3*B3-q4*B4-q5*B5-q6*B6;

//12 Portfolio Choice

//13-18 Bond Euler
uc*q=beta*uc(+1);
uc*q2=beta*uc(+1)*q(+1);
uc*q3=beta*uc(+1)*q2(+1);
uc*q4=beta*uc(+1)*q3(+1);
uc*q5=beta*uc(+1)*q4(+1);
uc*q6=beta*uc(+1)*q5(+1);

//19 LoM for volatility
sigma=sigma(-1)*rhosig+(1-rhosig)*log(.021) + usig;

end;



%----------------------------------------------------------------
% 4. Computation%----------------------------------------------------------------

initval;

y               	=       1;
c               	=       0.8 ;
k               	=      4.5 ;
i               	=       0.2 ;
l               	=       0.33 ;
w               	=       2.4448 ;
rk              	=       0.0240 ;

q=beta;
B1=0;
B6=0; 
end;

shocks;
var uep; stderr 0.6;
%var usig; stderr 0.015;
end;

maxit_=2000; 

steady(solve_algo=0);
stoch_simul( order=3, irf=40) y c i k;
gammab=0.01;
plot_policy_fun('sigma',[log(.011),log(.021),log(0.031),log(0.041),log(.051)],'q',oo_.dr.ys,1);
%plot_policy_fun('sigma',[log(.011),log(.021),log(0.031),log(0.041),log(.051)],'q1',oo_.dr.ys,1);
%plot_policy_fun('sigma',[log(.011),log(.021),log(0.031),log(0.041),log(.051)],'q2',oo_.dr.ys,1);
%plot_policy_fun('sigma',[log(.011),log(.021),log(0.031),log(0.041),log(.051)],'q3',oo_.dr.ys,1);
%plot_policy_fun('sigma',[log(.011),log(.021),log(0.031),log(0.041),log(.051)],'q3',oo_.dr.ys,1);
%plot_policy_fun('sigma',[log(.011),log(.021),log(0.031),log(0.041),log(.051)],'q3',oo_.dr.ys,1);
%plot_policy_fun('sigma',[log(.011),log(.021),log(0.031),log(0.041),log(.051)],'q6',oo_.dr.ys,1);


%----------------------------------------------------------------
% 5. Some Results
%----------------------------------------------------------------

%statistic1 = 100*sqrt(diag(oo_.var(1:6,1:6)))./oo_.mean(1:6);
%dyntable('Relative standard deviations in %',strvcat('VARIABLE','REL. S.D.'),M_.endo_names(1:6,:),statistic1,10,8,4);

%first order rules (declaration order)
 oo_.dr.ghx(oo_.dr.inv_order_var',:)


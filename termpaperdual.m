%termpaper data
tbl = readtable('data.csv');
mh = [];
mh_d = [];
for i=1:size(tbl.MH)
    if(mod(i,3)==1)
        if(i==1)
            mh = [mh; tbl.MH(i)];
        else
            mh = [mh; mh(int16(i/3))+tbl.MH(i)];
        end
        mh_d = [mh_d; tbl.MH(i)];
    end
end
x1 = linspace(1,size(mh,1),size(mh,1));

f_mh = fit(x1.',mh,'exp1');

ka = [];
ka_d = [];
for i=1:size(tbl.KA)
    if(mod(i,3)==1)
        if(i==1)
            ka = [ka; tbl.KA(i)];
        else
            ka = [ka; ka(int16(i/3))+tbl.KA(i)];
        end
        ka_d = [ka_d; tbl.KA(i)];
    end
end
f_ka = fit(x1.',ka,'exp1');

ap = [];
ap_d = [];
for i=1:size(tbl.AP)
    if(mod(i,3)==1)
        if(i==1)
            ap = [ap; tbl.AP(i)];
        else
            ap = [ap; ap(int16(i/3))+tbl.AP(i)];
        end
        ap_d = [ap_d; tbl.AP(i)];
    end
end
f_ap = fit(x1.',ap,'exp1');

alpha = [f_mh.b; f_ka.b; f_ap.b];
beta  = [f_mh.a; f_ka.a; f_ap.a];
n = [mh(size(mh,1)); ka(size(ka,1)); ap(size(ap,1))];
M = 30000;
ldown = [0;0;0];
gamma = 1.1; %Lockdown parameter
c = ldown.*(alpha./gamma)+(1-ldown).*(alpha);

m = sym('m', [3 1]);
k = sym('k', [3 1]);
p = sym('p');
d = sym('d');
u1 = sym('u1', [3 1]);
u2 = sym('u2', [3 1]);
u3 = sym('u3', [3 1]);
u4 = sym('u4', [3 1]);

l = sym('l');
l0 = 10;

M = 10;
m0 = [1;1;1];
k0 = [1;1;1];
%u0 = [0;0;0;0;0;0;0;0;0;0;0;0];
u01 = [0;0;0];
u02 = [0;0;0];
u03 = [0;0;0];
u04 = [0;0;0];

d0 = 1;
p0 = 10;
eta = 0.001;

%DUAL Ascent on LP
for i=20:size(mh,1)
   i
   f_mh = fit(x1(1:i).',mh(1:i),'exp1');
   f_ka = fit(x1(1:i).',ka(1:i),'exp1');
   f_ap = fit(x1(1:i).',ap(1:i),'exp1');
   
   alpha = [f_mh.b; f_ka.b; f_ap.b];
   beta  = [f_mh.a; f_ka.a; f_ap.a];
   n = [mh(i); ka(i); ap(i)];
   
   gamma = 0.2; %Lockdown parameter
   mu = 10;
   b = 100;
   
   f = (-l*double(M)-double(gamma*mu*sum(beta))+l*double(gamma*sum(beta.*exp(alpha.*i))));
   g = double(diff(f,l));
   l0 = l0+eta*subs(g,l,l0);
   if l0<0 
       l0 = 0;
   elseif l0>b
       l0 = b;
   end
   double(l0)
   
   %Analytically
   del_m = sum(gamma.*beta.*exp(alpha.*(i+1)))-M;
   p = l0*delm;
   k = gamme.*beta;
   m = k.*exp(alpha.*(i+1));
end

l1 = sym('l1');
l2 = sym('l2', [3,1]);
l20 = [1;1;1];
l10 = 10;
eta = 0.01;

%DUAL Ascent on Convex Objective
for i=20:size(mh,1)
   i
   f_mh = fit(x1(1:i).',mh(1:i),'exp1');
   f_ka = fit(x1(1:i).',ka(1:i),'exp1');
   f_ap = fit(x1(1:i).',ap(1:i),'exp1');
   
   alpha = [f_mh.b; f_ka.b; f_ap.b];
   beta  = [f_mh.a; f_ka.a; f_ap.a];
   n = [mh(i); ka(i); ap(i)];
   
   gamma = 0.2; 
   mu = 10;
   b = 100;
   
   f = (-l1*double(M)+l1*double(gamma*sum(beta.*exp(alpha.*i)))) ...
       + mu*(size(n,1) + sum(log(mu.*(l1.*exp(alpha.*i)-l2+0.0000000001))))...
       - mu*gamma*sum(beta./mu.*(l1.*exp(alpha.*i)-l2));   
   g = (diff(f,l1))
   g1= (gradient(f,l2))
   l10 = l10+eta*subs(subs(g,l1,l10),l2,l20);
   l20 = l20+eta*subs(subs(g1,l1,l10),l2,l20);

   %The lower bound is not limiting here, hence need not be as strictly
   %enforced. Hence slight relaxation on l2
   if l10<0 
       l10 = 0;
   elseif l10>b
       l10 = b;
   end
   l20(l20<0) = 0;
   l20(l20>b) = b;
   double([l10; l20])
   %Assuming this is the solution, using this for the complimentary slack
   %criterion. Choosing an appropriate step size solves this for this
   %problem!!!
   %Similar to LP, with analytical relaxation of constraint for
   %simplification
   del_m = sum(gamma.*beta.*exp(alpha.*(i+1)))-M;
   p = l10*delm;
   k = gamme.*beta;
   m = k.*exp(alpha.*(i+1));
   
end

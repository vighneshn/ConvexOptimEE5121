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
   
   f = (p-mu*sum(log(k)));
   f = f+u1.'*(gamma.*beta-k);
   f = f+u2.'*(k.*exp(alpha.*i)-m);
   f = f-u3.'*m;
   f = f+u4.'*[sum(m)-M-d; b*d-p; -d];
   g1 = gradient(f,m);
   g2 = gradient(f,k);
   g3 = diff(f,d);
   g4 = diff(f,p);
   g5 = gradient(f,u1);
   g6 = gradient(f,u2);
   g7 = gradient(f,u3);
   g8 = gradient(f,u4);

   m0
   thresh = 0.02;
   %while (1>0)
    fval = double(subs(subs(subs(subs(subs(subs(subs(subs(f,m,m0),k,k0),d,d0),p,p0),u1,u01),u2,u02),u3,u03),u4,u04));
    d1 = double(subs(subs(subs(subs(subs(subs(subs(subs(g1,m,m0),k,k0),d,d0),p,p0),u1,u01),u2,u02),u3,u03),u4,u04));
    d2 = double(subs(subs(subs(subs(subs(subs(subs(subs(g2,m,m0),k,k0),d,d0),p,p0),u1,u01),u2,u02),u3,u03),u4,u04));
    d3 = double(subs(subs(subs(subs(subs(subs(subs(subs(g3,m,m0),k,k0),d,d0),p,p0),u1,u01),u2,u02),u3,u03),u4,u04));
    d4 = double(subs(subs(subs(subs(subs(subs(subs(subs(g4,m,m0),k,k0),d,d0),p,p0),u1,u01),u2,u02),u3,u03),u4,u04));
    d5 = double(subs(subs(subs(subs(subs(subs(subs(subs(g5,m,m0),k,k0),d,d0),p,p0),u1,u01),u2,u02),u3,u03),u4,u04));
    d6 = double(subs(subs(subs(subs(subs(subs(subs(subs(g6,m,m0),k,k0),d,d0),p,p0),u1,u01),u2,u02),u3,u03),u4,u04));
    d7 = double(subs(subs(subs(subs(subs(subs(subs(subs(g7,m,m0),k,k0),d,d0),p,p0),u1,u01),u2,u02),u3,u03),u4,u04));
    d8 = double(subs(subs(subs(subs(subs(subs(subs(subs(g8,m,m0),k,k0),d,d0),p,p0),u1,u01),u2,u02),u3,u03),u4,u04));
    
    alpha = 0.1;
    beta = 0.5;
    eta1 = 1;
    eta2 = 1;
    
    while  double(subs(subs(subs(subs(subs(subs(subs(subs(f,m,m0-eta1*d1),k,k0-eta1*d2),d,d0-eta1*d3),p,p0-eta1*d4),u1,u01),u2,u02),u3,u03),u4,u04)) ...
            >fval-alpha*eta1*(dot(d1,d1)+dot(d2,d2)+dot(d3,d3)+dot(d4,d4))
        eta1 = beta*eta1; %Backtracking line search
        eta1
    end
    while  double(subs(subs(subs(subs(subs(subs(subs(subs(f,m,m0),k,k0),d,d0),p,p0),u1,u01+eta2*d5),u2,u02+eta2*d6),u3,u03+eta2*d7),u4,u04+eta2*d8)) ...
            < fval+alpha*eta2*(dot(d5,d5)+dot(d6,d6)+dot(d7,d7)+dot(d8,d8))
        eta2 = beta*eta2
    end
    eta1
    eta2
    %eta1*(sum(abs(d1./m0))+sum(abs(d2./k0))+abs(d3/d0)+abs(d4/p0));
    %if(eta1*(sum(abs(d1./m0))+sum(abs(d2./k0))+abs(d3/d0)+abs(d4/p0))<thresh)
    %    break
    %end
    m0  = m0  - d1*eta1;
    k0  = k0  - d2*eta1;
    d0  = d0  - d3*eta1;
    p0  = p0  - d4*eta1;
    u01 = u01 + d5*eta2;
    u02 = u02 + d6*eta2;
    u03 = u03 + d7*eta2;
    u04 = u04 + d8*eta2;

    u01(u01<0) = 0;
    u02(u02<0) = 0;
    u03(u03<0) = 0;
    u04(u04<0) = 0;
   %end
   m0;
   M = M+d0-sum(m0);
end


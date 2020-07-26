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

tn = [];
tn_d = [];
for i=1:size(tbl.TN)
    if(mod(i,3)==1)
        if(i==1)
            tn = [tn; tbl.TN(i)];
        else
            tn = [tn; tn(int16(i/3))+tbl.TN(i)];
        end
        tn_d = [tn_d; tbl.TN(i)];
    end
end

rj = [];
rj_d = [];
for i=1:size(tbl.RJ)
    if(mod(i,3)==1)
        if(i==1)
            rj = [rj; tbl.RJ(i)];
        else
            rj = [rj; rj(int16(i/3))+tbl.RJ(i)];
        end
        rj_d = [rj_d; tbl.RJ(i)];
    end
end

or = [];
or_d = [];
for i=1:size(tbl.OR)
    if(mod(i,3)==1)
        if(i==1)
            or = [or; tbl.OR(i)];
        else
            or = [or; or(int16(i/3))+tbl.OR(i)];
        end
        or_d = [or_d; tbl.OR(i)];
    end
end

wb = [];
wb_d = [];
for i=1:size(tbl.WB)
    if(mod(i,3)==1)
        if(i==1)
            wb = [wb; tbl.WB(i)];
        else
            wb = [wb; wb(int16(i/3))+tbl.WB(i)];
        end
        wb_d = [wb_d; tbl.WB(i)];
    end
end

M = 10;
M2 = [];
mh_plot = [0];
ka_plot = [0];
ap_plot = [0];
tn_plot = [0];
rj_plot = [0];
or_plot = [0];
wb_plot = [0];
mh_plot2 = [0];
ka_plot2 = [0];
ap_plot2 = [0];
tn_plot2 = [0];
rj_plot2 = [0];
or_plot2 = [0];
wb_plot2 = [0];
tot_plot = [0];
tot_plot2 = [0];
regret = [0];
for i=20:size(mh,1)
   
   f_mh = fit(x1(1:i).',mh(1:i),'exp1');
   f_ka = fit(x1(1:i).',ka(1:i),'exp1');
   f_ap = fit(x1(1:i).',ap(1:i),'exp1');
   f_tn = fit(x1(1:i).',tn(1:i),'exp1');
   f_rj = fit(x1(1:i).',rj(1:i),'exp1');
   f_or = fit(x1(1:i).',or(1:i),'exp1');
   f_wb = fit(x1(1:i).',wb(1:i),'exp1');
   
   alpha = [f_mh.b; f_ka.b; f_ap.b; f_tn.b; f_rj.b; f_or.b; f_wb.b];
   beta  = [f_mh.a; f_ka.a; f_ap.a; f_tn.a; f_rj.a; f_or.a; f_wb.a];
   n = [mh(i); ka(i); ap(i); tn(i); rj(i); or(i); wb(i)];
   
   %ldown = [0;0;0];
   gamma = 0.2; %Lockdown parameter
   mu = 10;
   b = 1000;
   %grad = subs(g,x,fin);
   %fin = fin - eta*grad;
   %Used cvx as this wasnt converging
   cvx_begin
        variable m(7,1)
        variable k(7,1)
        variable p
        variable del_mt
        
        minimize(p-mu*(sum(log(k))))
        subject to
            sum(m) < M+del_mt;
            k > gamma.*beta;
            m > k.*exp(alpha.*(i+1));
            p > b*(del_mt);
            del_mt > 0;
   cvx_end
   M = M+del_mt-sum(m);
   fin = m
   del_mt;
   mh_plot = [mh_plot m(1)];
   ka_plot = [ka_plot m(2)];
   ap_plot = [ap_plot m(3)];
   tn_plot = [tn_plot m(4)];
   rj_plot = [rj_plot m(5)];
   or_plot = [or_plot m(6)];
   wb_plot = [wb_plot m(7)];
   tot_plot = [tot_plot sum(m)];
   
   mh_plot2 = [mh_plot2 m(1)-mh_plot(length(mh_plot)-1)];
   ka_plot2 = [ka_plot2 m(2)-ka_plot(length(ka_plot)-1)];
   ap_plot2 = [ap_plot2 m(3)-ap_plot(length(ap_plot)-1)];
   tn_plot2 = [tn_plot2 m(4)-tn_plot(length(tn_plot)-1)];
   rj_plot2 = [rj_plot2 m(5)-rj_plot(length(rj_plot)-1)];
   or_plot2 = [or_plot2 m(6)-or_plot(length(or_plot)-1)];
   wb_plot2 = [wb_plot2 m(7)-wb_plot(length(wb_plot)-1)];
   tot_plot2 = [tot_plot2 sum(m)-mh_plot(length(mh_plot)-1)-ka_plot(length(ka_plot)-1)-ap_plot(length(ap_plot)-1)-tn_plot(length(tn_plot)-1)-rj_plot(length(rj_plot)-1)-or_plot(length(or_plot)-1)-wb_plot(length(wb_plot)-1)];
   M2 = [M2 M];
   
   tot_day = mh(i)+ka(i)+ap(i);
   %regret = [regret regret(length(mh_plot)-1)+(gamma*mh(i)/tot_day-m(1)/tot_day)^2+(gamma*ka(i)/tot_day-m(2)/tot_day)^2+(gamma*ap(i)/tot_day-m(3)/tot_day)^2];
end

fin

figure
plot(mh_plot);
hold on
plot(ka_plot);
hold on;
plot(ap_plot);
legend("MH","KA","AP")
xlabel("Day 20-120")
ylabel("Total medicine supplied to each state")

figure
plot(mh_plot2);
hold on;
plot(mh_d(20:end));
hold on;
legend("daily MH Med supplied","daily new MH cases")
xlabel("Day 20-120")
ylabel("Medicine supplied to MH")

figure
plot(ka_plot2);
hold on;
plot(ka_d(20:end));
legend("daily KA Med supplied","daily new KA cases")
xlabel("Day 20-120")
ylabel("Medicine supplied to KA")

figure
hold on;
plot(ap_plot2);
hold on;
plot(ap_d(20:end));
legend("daily AP Med supplied","daily new AP cases")
xlabel("Day 20-120")
ylabel("Medicine supplied to AP")

figure
plot(mh_d);
hold on;
plot(ka_d);
hold on;
plot(ap_d);
legend("Maharashtra", "Karnataka", "Andhra Pradesh")
xlabel("Days")
ylabel("Daily cases")
title("Daily cases based on covid19 data")

figure
plot(tot_plot);
legend("Total")
xlabel("Day 20-120")
ylabel("Total medicine supplied")

figure
plot(tot_plot2);
hold on;
plot(mh_plot2);
hold on;
plot(ka_plot2);
hold on;
plot(ap_plot2);
hold on;
plot(tn_plot2);
hold on;
plot(rj_plot2);
hold on;
plot(or_plot2);
hold on;
plot(wb_plot2);
legend("Total","MH","KA","AP","TN","RJ","OR","WB")
xlabel("Day 20-120")
ylabel("Daily medicine supplied")

figure
plot(regret);
xlabel("Day 20-120")
ylabel("Regret")

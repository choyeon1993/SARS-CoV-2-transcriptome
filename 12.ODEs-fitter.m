%% Script to solve equation 2 & equation 4 in manuscript, e1, e2, ed, k1 and k2 were estimated by randomly initalizing input of function 'solve' to find local optimazed estimations.

function [esol, e0best, ksol, k0best] = fit()
%% ---------------------- for non-drug treated samples ------------------

tspan = [0 24];
%y0 = [1e-6 0];
%yvalstrue = [1e-6 6; 0 6/16];
y0 = [1 0];
yvalstrue = [1 2.6e4; 0 1.6e3];
e = optimvar('e', 3, "LowerBound", 0, 'UpperBound',10);
myfcn = fcn2optimexpr(@EtoModel,e,tspan,y0);
obj = sum(sum((myfcn - yvalstrue).^2));
prob = optimproblem("Objective",obj);

%% Estimate with initial values 
% e0.e = [1 1 1];
% [esol,sumsq] = solve(prob, e0);
% display(esol.e)
% esol.e =[2.0246 0.0080 0.0237]

%% Estimate by an iteration approach
temp = 1e15;
for i = 1:1000
    e0.e =  [rand*10 rand*10 rand*10];
    [esol,sumsq] = solve(prob, e0);
    if sumsq < temp
        temp = sumsq;
        esolbest = esol;
        e0best = e0.e;
    end
end
disp(esolbest)
esol = esolbest %esol.e = [1.7531 0.0070 0.0073]

tspan2 = linspace(0, 24, 100);
fitfun = ode45(@(t,y)diffun(t,y, esol.e),tspan2,y0);
fity = deval(fitfun, tspan2);

plot(tspan2,fity(1,:))
hold on
plot(tspan2,fity(2,:))
hold on
obx=[0 0 24 24]
oby=[1 0 2.6e4 1.6e3]
scatter(obx,oby,'filled')
%     'MarkerSize',10,...
%     'MarkerEdgeColor','b',...
%     'MarkerFaceColor',[0.5,0.5,0.5])
hold off
title('Virus Treated Sample')
legend('Sense Strand', 'Anti-sense Strand', 'Observations', 'Location','northwest')
xlabel('Time')
ylabel('Coverage')



% ---------
e1= esol.e(1)
e2= esol.e(2)
ed= esol.e(3)
syms u(t) v(t)
ode1 = diff(u) == e2*v - ed*u;
ode2 = diff(v) == e1*u - ed*v;
odes = [ode1; ode2]

S = dsolve(odes)

uSol(t) = S.u
vSol(t) = S.v
[uSol(t), vSol(t)] = dsolve(odes)

cond1 = u(0) == 1;
cond2 = v(0) == 0;
conds = [cond1, cond2]
[uSol(t), vSol(t)] = dsolve(odes,conds)

%%  --------------------------- Drug treated samples ------
%y0 = [1 0];
tspan = [0 24];
%yvalstrue = [1 2; 0 2/20];
yvalstrue = [1 5.5e3; 0 1.65e2];
k = optimvar('k', 2, "LowerBound", 0, 'UpperBound',1);
myfcn = fcn2optimexpr(@KtoModel,esol.e,k,tspan,y0);
obj = sum(sum((myfcn - yvalstrue).^2));
prob = optimproblem("Objective",obj);

% k0.k = [0 0];
% [ksol,sumsq] = solve(prob,k0);
% disp(ksol.k)

temp = 1e15;
for i = 1:100
    k0.k =  [rand rand]
    [ksol,sumsq] = solve(prob, k0);
    if sumsq < temp
        temp = sumsq;
        ksolbest = ksol;
        k0best = k0.k
    end
end
ksol = ksolbest %ksol.k = [0.0343 0.0403]

tspan2 = linspace(0, 24, 100);
fitfun = ode45(@(t,y)diffun2(t,y, esol.e, ksol.k),tspan2,y0);
fity = deval(fitfun, tspan2);

plot(tspan2,fity(1,:))
hold on
plot(tspan2,fity(2,:))
hold on
obx=[0 0 24 24]
oby=[1 0 5.5e3 1.65e2]
scatter(obx,oby,'filled')
%     'MarkerSize',10,...
%     'MarkerEdgeColor','b',...
%     'MarkerFaceColor',[0.5,0.5,0.5])
hold off
title('RDV Treated Sample')
legend('Sense Strand', 'Anti-sense Strand', 'Observations', 'Location','northwest')
xlabel('Time')
ylabel('Coverage')

end


function dydt = diffun(~, y, e)
dydt = zeros(2,1);

dydt(1) = e(2)*y(2) - e(3)*y(1);
dydt(2) = e(1)*y(1) -e(3)*y(2);
end

function dydt = diffun2(t, y, e,k)
dydt = zeros(2,1);

dydt(1) = e(2)*(k(1)*t)*y(2) - e(3)*y(1);
dydt(2) = e(1)*(k(2)*t)*y(1) -e(3)*y(2);
end

function solpts = EtoModel(e, tspan, y0)
sol = ode45(@(t,y)diffun(t,y,e), tspan, y0);
solpts = deval(sol, tspan);
end


function solpts = KtoModel(e, k, tspan, y0)
sol = ode45(@(t,y)diffun2(t,y,e,k), tspan, y0);
solpts = deval(sol, tspan);
end


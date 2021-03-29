%This model evaluates the evolution of cooperation for strains with 1
%genotype in an increasing density-dependent manner.
%% Initiation
%This first part establishes the number of individuals (n on the paper) in each subpopulation 

format bank %
 
lambdavector=[2 3 4 5 6 7 8 9 10 50];

%% Creating subpopulations
%This second section creates one thousand subpopulations according to the
%Poisson distribution and assigns genotypes to each individual cells in all populations following
%the growth of the firsts subpopulations.
for i=1:length(lambdavector)
lambda=lambdavector(i);
xi= poissrnd(lambda,1,1000);
xi(xi == 0) = [ ]; 

for x = 1:numel(xi)
 gen{x}= randi([1,2],1,xi(x));
end

for z = 1:numel (gen)
gen1{z} = sum(gen{z}==1');
gen0{z} = sum(gen{z}==2');
end

gen1m = cell2mat(gen1); 
gen0m = cell2mat(gen0); 

nResult = length(gen1m);
Result  = cell(1, nResult);

for iResult = 1:nResult
    
  a = gen1m(iResult); 
  b = gen0m(iResult);;
   
tspan = [0 600];
y0=[a b];
[t,x] = ode45( @dstate ,tspan ,y0);
          iResult = iResult + 1;
          Result{iResult} = [x];
end


%% Structured population cycle
%This section evaluates n growth cycles for all subpopulations and finally
%saves the population frequency for each genotype as well as the population density.

n= 999

for  nt = 1:n
    
for w = 2:numel(Result)
popul{w}= Result{w}(end,:); 
end
for y = 2:numel(popul)
gen1pop{y}= popul{y}(:,1);
gen0pop{y}= popul{y}(:,2); 
end

popgen1 = cell2mat(gen1pop);
maximo1 = find(abs(popgen1)>1000000);
popgen1(maximo1) = [];

popgen0 = cell2mat(gen0pop);
maximo2 = find(abs(popgen0)>1000000);
popgen0(maximo2) = [];


gen1Popul= round(sum(popgen1,'omitnan')); 
gen0Popul= round(sum(popgen0,'omitnan'));


xi= poissrnd(lambda,1,1000);
xi(xi == 0) = [ ]; 

Total= sum([gen1Popul,gen0Popul]);
gen1f= gen1Popul / Total; 
gen0f= gen0Popul / Total; 

for x = 1:numel(xi)
newgen{x} = randsrc(1,xi(x),[1,2;gen1f,gen0f]);
end

for z2 = 1:numel (newgen)
newgen1{z2} = sum(newgen{z2}==1');
newgen0{z2} = sum(newgen{z2}==2');
end

newgen1m = cell2mat(newgen1);
newgen0m = cell2mat(newgen0);


nresult = length(newgen1m);
result  = cell(1, nresult);

for inresult = 1:nresult;
 a1 = newgen1m(inresult);
 b1 = newgen0m(inresult);
  
y0=[a1 b1];
[t,x1] = ode45( @dstate ,tspan ,y0);
          inresult = inresult + 1;
          result{inresult} = [x1];
end

 Total= sum([gen1Popul,gen0Popul]); 
    Gen1freq{nt}= gen1Popul /  Total;
    Gen0freq{nt}= gen0Popul / Total;
    
    Gen1Popul{nt}= gen1Popul;
    Gen0Popul{nt}= gen0Popul;
    Result= result;
end

 Gen1freq= cell2mat(Gen1freq);
 Gen0freq= cell2mat(Gen0freq);
 
filename = [ 'var ' num2str(lambda) '.mat' ];
save(filename,'Gen1freq','Gen0freq','Gen1Popul','Gen0Popul');
clear;
lambdavector=[2 3 4 5 6 7 8 9 10 50];
clc;
end

function dydt = dstate (t,x);
r = 1;
c= 0.05;
delta= 0.025;
kmax= 100000;
kbase= 1;
den= x(1) + x(2);
r1= r - c; 
kcoop= kmax*(x(1)/den);
kd= kmax/2;
g= (delta * den)/(kd + den); 

dydt = [(r1*x(1)*(1 - (den /((r1/r)*(kbase +((1 + g)*kcoop)))))); 
        (r*x(2)*(1 - (den /(kbase +(((1 - g)*kcoop))))))];
end

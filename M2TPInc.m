%This model evaluates the evolution of cooperation for strains with 2
%genotypes in an increasing density-dependent manner.
%% Initiation
%This first part establishes the number of individuals (n on the paper) in each subpopulation 

format bank  
lambdavector=[2 3 4 5 6 7 8 9 10 50];

for i=1:length(lambdavector)    
%% Creating subpopulations
%This second section creates one thousand subpopulations according to the
%Poisson distribution and assigns genotypes to each individual cells in all populations following
%the growth of the firsts subpopulations. 

lambda= lambdavector(i);
xi= poissrnd(lambda,1,1000);
xi(xi == 0) = [ ];

for x = 1:numel(xi)
gen{x}= randi([1,4],1,xi(x));
end

for z = 1:numel (gen)
gen11{z} = sum(gen{z}==1');
end
for z = 1:numel (gen)
gen10{z} = sum(gen{z}==2');
end
for z = 1:numel (gen)
gen01{z} = sum(gen{z}==3');
end
for z = 1:numel (gen)
gen00{z} = sum(gen{z}==4');
end

gen11m = cell2mat(gen11); 
gen10m = cell2mat(gen10); 
gen01m = cell2mat(gen01);
gen00m = cell2mat(gen00);


nResult = length(gen10m);
Result  = cell(1, nResult);

for iResult = 1:nResult

  a = gen11m(iResult); 
  b = gen10m(iResult);
  c = gen01m(iResult);
  d = gen00m(iResult);
  
tspan = [0 600];
y0=[a b c d];
[t,x] = ode45( @dstate ,tspan ,y0);
          iResult = iResult + 1;
          Result{iResult} = [x];
end


%% Structured deme cycle
%This section evaluates n growth cycles for all subpopulations and finally
%saves the population frequency for each genotype as well as the population density.

n= 999

for  nt = 1:n
 
for w = 2:numel(Result)
popul{w}= Result{w}(end,:); 
end

for y = 2:numel(popul)
gen11pop{y}= popul{y}(:,1); 
gen10pop{y}= popul{y}(:,2); 
gen01pop{y}= popul{y}(:,3); 
gen00pop{y}= popul{y}(:,4); 
end

popgen11 = cell2mat(gen11pop); 
maximo1 = find(abs(popgen11)>100000);
popgen11(maximo1) = 10000;

popgen10 = cell2mat(gen10pop);
maximo2 = find(abs(popgen10)>100000);
popgen10(maximo2) = 10000;

popgen01 = cell2mat(gen01pop);
maximo3 = find(abs(popgen01)>100000);
popgen01(maximo3) = 10000;

popgen00 = cell2mat(gen00pop);
maximo4 = find(abs(popgen00)>100000);
popgen00(maximo4) = 10000;

gen11Popul= round(sum(popgen11,'omitnan'));  
gen10Popul= round(sum(popgen10,'omitnan'));
gen01Popul= round(sum(popgen01,'omitnan'));
gen00Popul= round(sum(popgen00,'omitnan'));

xi= poissrnd(lambda,1,1000);
xi(xi == 0) = [ ];

Total= sum([gen11Popul,gen10Popul,gen01Popul,gen00Popul]);
gen11f= gen11Popul / Total; 
gen10f= gen10Popul / Total; 
gen01f= gen01Popul / Total; 
gen00f= gen00Popul / Total; 
for x = 1:numel(xi)
newgen{x} = randsrc(1,xi(x),[1,2,3,4;gen11f,gen10f,gen01f,gen00f]);
end

for z = 1:numel (newgen)
newgen11{z} = sum(newgen{z}==1');
newgen10{z} = sum(newgen{z}==2'); 
newgen01{z} = sum(newgen{z}==3');
newgen00{z} = sum(newgen{z}==4');
end

newgen11m = cell2mat(newgen11);
newgen10m = cell2mat(newgen10);
newgen01m = cell2mat(newgen01);
newgen00m = cell2mat(newgen00);


nresult = length(newgen11m);
result  = cell(1, nresult);
for inresult = 1:nresult

  a1 = newgen11m(inresult);
  b1 = newgen10m(inresult);
  c1 = newgen01m(inresult);
  d1 = newgen00m(inresult);
  
y0=[a1 b1 c1 d1];
[t,x1] = ode45( @dstate ,tspan ,y0);
          inresult = inresult + 1;
          result{inresult} = [x1];
end

Total= sum([gen11Popul,gen10Popul,gen01Popul,gen00Popul]);
    Gen11freq{nt}= gen11Popul / Total;
    Gen10freq{nt}= gen10Popul / Total;
    Gen01freq{nt}= gen01Popul / Total;
    Gen00freq{nt}= gen00Popul / Total;
    Result= result;
   generations{nt}= newgen{nt};
end

 Gen11freq= cell2mat(Gen11freq);
 Gen10freq= cell2mat(Gen10freq);
 Gen01freq= cell2mat(Gen01freq);
 Gen00freq= cell2mat(Gen00freq);

 
filename = [ 'var' num2str(lambda) '.mat' ];
save(filename,'Gen11freq','Gen10freq','Gen01freq','Gen00freq','gen11Popul','gen10Popul','gen01Popul','gen00Popul','generations');
clear;
lambdavector=[2 3 4 5 6 7 8 9 10 50];  
clc;
end

function dydt = dstate (t,x);
r = 1;
v = 1;
c= 0.1;
delta= 0.2;
kmax= 100000;
kbase= 1;
den= x(1) + x(2) +x(3) + x(4);
r11= r - (((1 + 1)^v)*c);
r10= r - (((1 + 0)^v)*c);
r01= r - (((0 + 1)^v)*c);
r00= r - (((0 + 0)^v)*c);
kcoop= kmax*((x(1) + x(2))/den)*((x(1) + x(3))/den);
kd= kmax/2;
g= (den * delta)/(kd + den); 

dydt = [(r11*x(1)*(1 - (den / ((r11/r)*(kbase +(((1 + (((2*1)-1)*g))*(1 + (((2*1)-1)*g)))*kcoop)))))); 
        (r10*x(2)*(1 - (den / ((r10/r)*(kbase +(((1 + (((2*1)-1)*g))*(1 + (((2*0)-1)*g)))*kcoop)))))); 
        (r01*x(3)*(1 - (den / ((r01/r)*(kbase +(((1 + (((2*0)-1)*g))*(1 + (((2*1)-1)*g)))*kcoop))))));
        (r00*x(4)*(1 - (den / ((r00/r)*(kbase +(((1 + (((2*0)-1)*g))*(1 + (((2*0)-1)*g)))*kcoop))))))];
end

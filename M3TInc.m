%This model evaluates the evolution of cooperation for strains with 3
%genotypes in an increasing density-dependent manner.
%% Initiation
%This first part establishes the number of individuals (n on the paper) in each subpopulation 

format bank 
  
lambdavector=[2 3 4 5 6 7 8 9 10 50];

%% Creating subpopulations
%This second section creates one thousand subpopulations according to the
%Poisson distribution and assigns genotypes to each individual cells in all populations following
%the growth of the firsts subpopulations. 

for i=1:length(lambdavector) 
    
lambda= lambdavector(i); 
xi= poissrnd(lambda,1,1000);
xi(xi == 0) = [ ]; 

for x = 1:numel(xi)
 gen{x}= randi([1,8],1,xi(x));
end

for z = 1:numel (gen)
gen111{z} = sum(gen{z}==1');
gen110{z} = sum(gen{z}==2');
gen101{z} = sum(gen{z}==3');
gen011{z} = sum(gen{z}==4');
gen100{z} = sum(gen{z}==5');
gen010{z} = sum(gen{z}==6');
gen001{z} = sum(gen{z}==7');
gen000{z} = sum(gen{z}==8');
end

gen111m = cell2mat(gen111); 
gen110m = cell2mat(gen110); 
gen101m = cell2mat(gen101);
gen011m = cell2mat(gen011);
gen100m = cell2mat(gen100);
gen010m = cell2mat(gen010);
gen001m = cell2mat(gen001);
gen000m = cell2mat(gen000);

nResult = length(gen111m);
Result  = cell(1, nResult);

for iResult = 1:nResult

  a = gen111m(iResult); 
  b = gen110m(iResult);
  c = gen101m(iResult);
  d = gen011m(iResult);
  e = gen100m(iResult);
  f = gen010m(iResult);
  g = gen001m(iResult);
  h = gen000m(iResult);
  
tspan = [0 600];
y0=[a b c d e f g h];
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
gen111pop{y}= popul{y}(:,1); 
gen110pop{y}= popul{y}(:,2); 
gen101pop{y}= popul{y}(:,3); 
gen011pop{y}= popul{y}(:,4); 
gen100pop{y}= popul{y}(:,5); 
gen010pop{y}= popul{y}(:,6); 
gen001pop{y}= popul{y}(:,7); 
gen000pop{y}= popul{y}(:,8); 
end

popgen111 = cell2mat(gen111pop);  
maximo1 = find(abs(popgen111)>200000);
minimo1 = find(abs(popgen111)<0.5);
popgen111(maximo1) = 10000;
popgen111(minimo1) = 0;

popgen110 = cell2mat(gen110pop);
maximo2 = find(abs(popgen110)>200000);
minimo2 = find(abs(popgen110)<0.5);
popgen110(maximo2) = 10000;
popgen110(minimo2) = 0;

popgen101 = cell2mat(gen101pop);
maximo3 = find(abs(popgen101)>200000);
minimo3 = find(abs(popgen101)<0.5);
popgen101(maximo3) = 10000;
popgen101(minimo3) = 0;

popgen011 = cell2mat(gen011pop);
maximo4 = find(abs(popgen011)>200000);
minimo4 = find(abs(popgen011)<0.5);
popgen011(maximo4) = 10000;
popgen011(minimo4) = 0;

popgen100 = cell2mat(gen100pop);
maximo5 = find(abs(popgen100)>200000);
minimo5 = find(abs(popgen100)<0.5);
popgen100(maximo5) = 10000;
popgen100(minimo5) = 0;

popgen010 = cell2mat(gen010pop);
maximo6 = find(abs(popgen010)>200000);
minimo6 = find(abs(popgen010)<0.5);
popgen010(maximo6) = 10000;
popgen010(minimo6) = 0;

popgen001 = cell2mat(gen001pop);
maximo7 = find(abs(popgen001)>200000);
minimo7 = find(abs(popgen001)<0.5);
popgen011(maximo7) = 10000;
popgen011(minimo7) = 0;

popgen000 = cell2mat(gen000pop);
maximo8 = find(abs(popgen000)>200000);
minimo8 = find(abs(popgen000)<0.5);
popgen000(maximo8) = 10000;
popgen000(minimo8) = 0;

gen111Popul= round(sum(popgen111,'omitnan'));  
gen110Popul= round(sum(popgen110,'omitnan'));
gen101Popul= round(sum(popgen101,'omitnan'));
gen011Popul= round(sum(popgen011,'omitnan'));
gen100Popul= round(sum(popgen100,'omitnan')); 
gen010Popul= round(sum(popgen010,'omitnan'));
gen001Popul= round(sum(popgen001,'omitnan'));
gen000Popul= round(sum(popgen000,'omitnan'));


xi= poissrnd(lambda,1,1000);
xi(xi == 0) = [ ]; 

Total= sum([gen111Popul,gen110Popul,gen101Popul,gen011Popul,gen100Popul,gen010Popul,gen001Popul,gen000Popul]);
gen111f= gen111Popul / Total; 
gen110f= gen110Popul / Total; 
gen101f= gen101Popul / Total; 
gen011f= gen011Popul / Total; 
gen100f= gen100Popul / Total;
gen010f= gen010Popul / Total;
gen001f= gen001Popul / Total;
gen000f= gen000Popul / Total;

for x = 1:numel(xi)
newgen{x} = randsrc(1,xi(x),[1,2,3,4,5,6,7,8;gen111f,gen110f,gen101f,gen011f,gen100f,gen010f,gen001f,gen000f]);
end

for z2 = 1:numel (newgen)
newgen111{z2} = sum(newgen{z2}==1');
newgen110{z2} = sum(newgen{z2}==2');
newgen101{z2} = sum(newgen{z2}==3');
newgen011{z2} = sum(newgen{z2}==4');
newgen100{z2} = sum(newgen{z2}==5');
newgen010{z2} = sum(newgen{z2}==6');
newgen001{z2} = sum(newgen{z2}==7');
newgen000{z2} = sum(newgen{z2}==8');
end

newgen111m = cell2mat(newgen111);
newgen110m = cell2mat(newgen110);
newgen101m = cell2mat(newgen101);
newgen011m = cell2mat(newgen011);
newgen100m = cell2mat(newgen100);
newgen010m = cell2mat(newgen010);
newgen001m = cell2mat(newgen001);
newgen000m = cell2mat(newgen000);

nresult = length(newgen111m);
result  = cell(1, nresult);
for inresult = 1:nresult

  a1 = newgen111m(inresult);
  b1 = newgen110m(inresult);
  c1 = newgen101m(inresult);
  d1 = newgen011m(inresult);
  e1 = newgen100m(inresult);
  f1 = newgen010m(inresult);
  g1 = newgen001m(inresult);
  h1 = newgen000m(inresult);

y0=[a1 b1 c1 d1 e1 f1 g1 h1];
[t,x1] = ode45( @dstate ,tspan ,y0);
          inresult = inresult + 1;
          result{inresult} = [x1];
end

 Total= sum([gen111Popul,gen110Popul,gen101Popul,gen011Popul,gen100Popul,gen010Popul, gen001Popul,gen000Popul]);
 
    Gen111freq{nt}= gen111Popul / Total;
    Gen101freq{nt}= gen101Popul / Total;
    Gen110freq{nt}= gen110Popul / Total;
    Gen011freq{nt}= gen011Popul / Total;
    Gen100freq{nt}= gen100Popul / Total;
    Gen010freq{nt}= gen010Popul / Total;
    Gen001freq{nt}= gen001Popul / Total;
    Gen000freq{nt}= gen000Popul / Total;
    Gen111Popul{nt}= gen111Popul;
    Gen101Popul{nt}= gen101Popul;
    Gen110Popul{nt}= gen110Popul;
    Gen011Popul{nt}= gen011Popul;
    Gen100Popul{nt}= gen100Popul;
    Gen010Popul{nt}= gen010Popul;
    Gen001Popul{nt}= gen001Popul;
    Gen000Popul{nt}= gen000Popul;
    Result= result;
end

 Gen111freq = cell2mat(Gen111freq);
 Gen101freq = cell2mat(Gen101freq);
 Gen110freq = cell2mat(Gen110freq);
 Gen011freq = cell2mat(Gen011freq);
 Gen100freq = cell2mat(Gen100freq);
 Gen010freq = cell2mat(Gen010freq);
 Gen001freq = cell2mat(Gen001freq);
 Gen000freq = cell2mat(Gen000freq);
 
filename= [ 'var' num2str(lambda) '.mat' ];
save(filename,'Gen111freq','Gen101freq','Gen110freq','Gen011freq','Gen100freq','Gen010freq','Gen001freq','Gen000freq','Gen111Popul','Gen101Popul','Gen110Popul','Gen011Popul','Gen100Popul','Gen010Popul','Gen001Popul','Gen000Popul');
clear;
lambdavector=[2 3 4 5 6 7 8 9 10 50];
clc;

end

function dydt = dstate (t,x);
r = 1;
c= 0.1; 
kmax= 100000;
v= 1;
kbase= 1;
delta= 0.2; 
den= x(1) + x(2) + x(3) + x(4) + x(5) + x(6) + x(7) + x(8);
freqgen1= x(1) + x(2) + x(3) + x(5);
freqgen2= x(1) + x(2) + x(4) + x(6);
freqgen3= x(1) + x(3) + x(4) + x(7);  
y=((freqgen1/den)*(freqgen2/den)*(freqgen3/den));
r111= r -(((1 + 1 + 1)^v)*c);
r110= r -(((1 + 1 + 0)^v)*c);
r101= r -(((1 + 0 + 1)^v)*c);
r011= r -(((0 + 1 + 1)^v)*c);
r100= r -(((1 + 0 + 0)^v)*c);
r010= r -(((0 + 1 + 0)^v)*c);
r001= r -(((0 + 0 + 1)^v)*c);
r000= r -(((0 + 0 + 0)^v)*c);
kcoop= kmax*(y);
kd= kmax/2;
g= (den * delta)/(kd + den);  

dydt = [(r111*x(1)*(1 - (den / ((r111/r)*(kbase + (((1 + (((2*1)-1)*g))*(1 + (((2*1)-1)*g))*(1 + (((2*1)-1)*g))) * kcoop))))));
        (r110*x(2)*(1 - (den / ((r110/r)*(kbase + (((1 + (((2*1)-1)*g))*(1 + (((2*1)-1)*g))*(1 + (((2*0)-1)*g))) * kcoop))))));
        (r101*x(3)*(1 - (den / ((r101/r)*(kbase + (((1 + (((2*1)-1)*g))*(1 + (((2*0)-1)*g))*(1 + (((2*1)-1)*g))) * kcoop))))));
        (r011*x(4)*(1 - (den / ((r011/r)*(kbase + (((1 + (((2*0)-1)*g))*(1 + (((2*1)-1)*g))*(1 + (((2*1)-1)*g))) * kcoop))))));
        (r100*x(5)*(1 - (den / ((r100/r)*(kbase + (((1 + (((2*1)-1)*g))*(1 + (((2*0)-1)*g))*(1 + (((2*0)-1)*g))) * kcoop))))));
        (r010*x(6)*(1 - (den / ((r010/r)*(kbase + (((1 + (((2*0)-1)*g))*(1 + (((2*1)-1)*g))*(1 + (((2*0)-1)*g))) * kcoop))))));
        (r001*x(7)*(1 - (den / ((r001/r)*(kbase + (((1 + (((2*0)-1)*g))*(1 + (((2*0)-1)*g))*(1 + (((2*1)-1)*g))) * kcoop))))));
        (r000*x(8)*(1 - (den / ((r000/r)*(kbase + (((1 + (((2*0)-1)*g))*(1 + (((2*0)-1)*g))*(1 + (((2*0)-1)*g))) * kcoop))))))];
end
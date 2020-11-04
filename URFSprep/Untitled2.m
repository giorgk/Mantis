%% N Base
N = linspace(10,20,15);
tmp = linspace(20,17,16);
N = [N tmp(2:end)];
tmp = linspace(17,25,16);
N = [N tmp(2:end)];
tmp = linspace(25,22,16);
N = [N tmp(2:end)];

%% N reduced
Nred = linspace(10,0.3*20,15);
tmp = linspace(0.3*20,0.5*17,16);
Nred = [Nred tmp(2:end)];
tmp = linspace(0.5*17,25*0.7,16);
Nred = [Nred tmp(2:end)];
tmp = linspace(25*0.7,22,16);
Nred = [Nred tmp(2:end)];
%%
pp = linspace(1,0.3,15);
tmp = linspace(0.3,0.5,16);
pp = [pp tmp(2:end)];
tmp = linspace(0.5,0.7,16);
pp = [pp tmp(2:end)];
tmp = linspace(0.7,1,16);
pp = [pp tmp(2:end)];
%%
clf
a = zeros(1,length(N));
a(25:40) = linspace(0,1,16);
a(41:end) = 1;
subplot(2,1,1);plot(a)
subplot(2,1,2);plot(N)
subplot(2,1,2);hold on
ylim([0 30])
subplot(2,1,2);plot(Nred,'r')
subplot(2,1,2);plot(N.*(1-a)+ Nred.*a,'g')
subplot(2,1,2);plot(N.*a.*pp,'c')
%% Read test client results
fid = fopen('CPP\TestClient\TestClient\testClientResults.dat','r');
CC = textscan(fid,'%f');
fclose(fid);
Conc = reshape(CC{1,1},155,1725)';
pp = prctile(Conc,[5 10:10:90 95],1);
plot(Conc','color',[0.5 0.5 0.5])
hold on
plot(pp')




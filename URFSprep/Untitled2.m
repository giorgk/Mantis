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
%%
color_order = colororder;
start_date = datetime(1945,1,1);
tm = start_date + calyears(0:254);
%% Read test client results
fid = fopen('f:\UCDAVIS\Mantis\CPP\TestClient\testClientResults.dat','r');
CC = textscan(fid,'%f');
fclose(fid);
nTimes = 255;
Nbtc = 1258;
SWAT1_Conc_Farm2_BUD1 = reshape(CC{1,1},nTimes,Nbtc)';
Conc = reshape(CC{1,1},nTimes,Nbtc)';
pp = prctile(Conc,[5 10:10:90 95],1);
%plot(Conc','color',[0.5 0.5 0.5])
%%
figure(2)
clf
prcnts = [50 75 90 95];
pp0 = prctile(SWAT1_Conc_Farm2_BUD0,prcnts,1);
plot(tm(1:156), pp0(:,1:156)','color',color_order(1,:), 'linewidth', 2)
hold on
pp1 = prctile(SWAT1_Conc_Farm2_BUD1,prcnts,1);
plot(tm(1:156), pp1(:,1:156)','color',color_order(2,:), 'linewidth', 2)
%legend('location','northwest')
grid on
ylabel('Concentration [mg/l]')
xlabel('Time')

for ii = 1:length(prcnts)
    text(tm(156),pp1(ii,156),[' ' num2str(prcnts(ii)) '%'])
    text(tm(156),pp0(ii,156),[' ' num2str(prcnts(ii)) '%'])
end

%patch([1945:2199 fliplr(1945:2199)], [pp0(11,:) fliplr(pp1(11,:))], 'g')



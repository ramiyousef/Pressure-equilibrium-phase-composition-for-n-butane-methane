%To  run this code simply click Run and the results should show after...
%aproximately 6.5 seconds
% Term Project Fall 2015
% copy right Rami Yousef Khalil 260558325
% CHEE390 Computational Methods in Chemical Engineeing
% Prof.Philip Servio
% Date: December 4 2015
T = 310.928;
xcomposition=zeros(1,1000);
ycomposition=zeros(1,1000);
Pressure=zeros(1,1000);
count = 1;
%% for loop to find x1 and y1 at the very begining 51 psi to 53 psi
for P = 51.99*6894.76: 100 : 53*6894.76
    z=[0.001 ; 0.999];
    [x , y,check] = flash(P,z,T);%using flash claculations to output the composition at each pressure
    
    if check == 1
        xcomposition(1,count)=x(1);
        ycomposition(1,count)=y(1);
        Pressure(1,count)=P;
        count = count + 1;
    end
end
%% for loop to find the x1 and y1 from 66.5 psi to 600 psi
for P = 66.5*6894.76: 10000 : 600*6894.76
    z = [0.2;0.8];
    [x , y,check] = flash(P,z,T);%using flash claculations to output the composition at each pressure
    
    if check == 1
        xcomposition(1,count)=x(1);
        ycomposition(1,count)=y(1);
        Pressure(1,count)=P;
        count = count + 1;
    end
end
%% for loop to find the x1 and y1 from 600 psi to 1194 psi
for P = 600*6894.76 : 10000 : 1194*6894.76
    z = [0.4;0.6];
    [x , y,check] = flash(P,z,T);%using flash claculations to output the composition at each pressure
    
    if check == 1
        xcomposition(1,count)=x(1);
        ycomposition(1,count)=y(1);
        Pressure(1,count)=P;
        count = count + 1;
    end
end
%% for loop to find the x1 and y1 from 1194 psi to 2110 psi
for P = 1194*6894.76 : 10000 : 2110*6894.76
    z=[0.75,0.25];
    [x , y,check ] = flash(P,z,T);%using flash claculations to output the composition at each pressure
    
    if check == 1
        xcomposition(1,count)=x(1);
        ycomposition(1,count)=y(1);
        Pressure(1,count)=P;
        count = count + 1;
        
    end
end
%%
Pressure(1,:)=Pressure(1,:)/6894.76;% convering from pa to psi

%% plotting
[z]=semilogy(xcomposition,Pressure,ycomposition,Pressure);
str = 'Two Phase Region';
str1='Liquid Phase';
str2='Vapor Phase';
str3='T=100^oF';
annotation('textbox',[.45 .3 .3 .3],'String',str,'FitBoxToText','on','fontsize',17,'color','k','linestyle','none','fontweight','bold');
annotation('textbox',[.25 .45 .3 .3],'String',str1,'FitBoxToText','on','fontsize',17,'color','b','linestyle','none','fontweight','bold');
annotation('textbox',[.65 .05 .3 .3],'String',str2,'FitBoxToText','on','fontsize',17,'color','r','linestyle','none','fontweight','bold');
annotation('textbox',[.15 .71 0 .2],'String',str3,'FitBoxToText','on','fontsize',17,'color','k','fontweight','bold');
grid on
set(gca,'box','on','TickDir','out','YTick',100:100:2000,'fontsize',13,'ylim',[49 2050],'xlim',[0 1])
set (z, 'LineWidth', 6);
title('Pressure Vs. x_1^m^e^t^h^a^n^e','fontsize',18,'fontangle','normal')
xlabel('x_1^m^e^t^h^a^n^e','fontsize',15,'fontangle','normal','fontweight','bold')
ylabel('Pressure (psia)','fontsize',15,'fontangle','normal','fontweight','bold')

hlegend=legend('Liquid','Vapor');

set(hlegend,'fontsize',13,'box','on','units','normalized','position',[0.27 0.8 0 0.2],'fontangle','normal','orientation','horizontal')

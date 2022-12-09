%% Le dados de irrad., Gera parametros alfa e beta, gera  curvas de irrad., calcula pot. gerada
clear;clc;

% Le dados do excel de 2017
RawData = csvread('109749_41.13_-8.62_2017_without_header.csv');
irradJan = RawData(1:2976, 2); tempJan = RawData(1:2976, 1);
irradJuly = RawData(17377:20352, 2); tempJuly = RawData(17377:20352, 1);

% Irradiancia em kW/m2
IrradJan = [irradJan(577:672) irradJan(1441:1536) irradJan(1537:1632) irradJan(1633:1728) irradJan(1729:1824) irradJan(1825:1920) irradJan(1921:2016)];
IrradJan = IrradJan./1000;
IrradJuly = [irradJuly(1:96) irradJuly(97:192) irradJuly(193:288) irradJuly(289:384) ...
        irradJuly(961:1056) irradJuly(1057:1152) irradJuly(1153:1248) irradJuly(1441:1536) ...
        irradJuly(1921:2016) irradJuly(2113:2208) irradJuly(2401:2496)  irradJuly(2497:2592) ];
IrradJuly = IrradJuly./1000;

% Calcula a irradiância media Jan e July
irrad_medJan = mean(IrradJan,2); 
irrad_medJuly = mean(IrradJuly,2);   

% Calcula o desvio padrão da irradiancia Jan e July
irrad_stdJan = std(IrradJan,1,2);  
irrad_stdJuly = std(IrradJuly,1,2);  

% Temperatura media em °C Jan e July 2017
TempJan = [tempJan(577:672) tempJan(1441:1536) tempJan(1537:1632) tempJan(1633:1728) tempJan(1729:1824) tempJan(1825:1920) tempJan(1921:2016)];
TempJuly = [tempJuly(1:96) tempJuly(97:192) tempJuly(193:288) tempJuly(289:384) ...
        tempJuly(961:1056) tempJuly(1057:1152) tempJuly(1153:1248) tempJuly(1441:1536) ...
        tempJuly(1921:2016) tempJuly(2113:2208) tempJuly(2401:2496)  tempJuly(2497:2592) ];
Temp_medJan(1,:) = mean(TempJan,2); Temp_medJuly(1,:) = mean(TempJuly,2); 

% Calcula alfa e beta da PDF para 1 dia de Jan
for i = 1:96
    if irrad_stdJan(i,1) == 0
        beta_Jan(i,1) = 0;
        alpha_Jan(i,1) = 0;
    else
        beta_Jan(i,1) = (1-irrad_medJan(i))*((irrad_medJan(i)*(1-irrad_medJan(i)))/(irrad_stdJan(i)^2) - 1);
        alpha_Jan(i,1) = (irrad_medJan(i)*beta_Jan(i))/(1-irrad_medJan(i));
    end
end

% Calcula alfa e beta da PDF para 1 dia de July
for i = 1:96
    if irrad_stdJuly(i,1) == 0
        beta_July(i,1) = 0;
        alpha_July(i,1) = 0;
    else
        beta_July(i,1) = (1-irrad_medJuly(i))*((irrad_medJuly(i)*(1-irrad_medJuly(i)))/(irrad_stdJuly(i)^2) - 1);
        alpha_July(i,1) = (irrad_medJuly(i)*beta_July(i))/(1-irrad_medJuly(i));
    end
end

% Gera  curvas de irrad. para 1 dia usando Beta PDF
N = 1000; 
for i=1:96
    Irrad_w(:,i) = betarnd(alpha_Jan(i),beta_Jan(i),N,1); 
    Irrad_s(:,i) = betarnd(alpha_July(i),beta_July(i),N,1); 
end

Irrad_w(isnan(Irrad_w))=0; Irrad_s(isnan(Irrad_s))=0;


% Dados do sistema fotovoltaico
Pmpref = 780;  % em (Wp)
NOCT   = 45;   % Nominal operating cell temperature (Datasheet)
n      = 220;   % Numero de m�dulos 
ninv   = 0.95; %Eficiencia do inversor

% Calcula Pot. gerada pela FV
Irrad_w = Irrad_w.*1000; Irrad_s = Irrad_s.*1000; % converte kW/m2 para W/m2
for mc = 1:N
    for hora = 1:96
        tc_w(mc,hora) = Temp_medJan(1,hora)+((NOCT-20)/800)*Irrad_w(mc,hora);
        tc_s(mc,hora) = Temp_medJuly(1,hora)+((NOCT-20)/800)*Irrad_s(mc,hora);
        Pmp_w(mc,hora) = n*Pmpref*(Irrad_w(mc,hora)/1000)*[1-0.0043*(tc_w(mc,hora)-25)]; 
        Pmp_s(mc,hora) = n*Pmpref*(Irrad_s(mc,hora)/1000)*[1-0.0043*(tc_s(mc,hora)-25)]; 
        Ppv_w(mc,hora) = (Pmp_w(mc,hora)*ninv)*10^-3; %Potencia de sa�da da usina Fotovoltaica (em kW)
        Ppv_s(mc,hora) = (Pmp_s(mc,hora)*ninv)*10^-3;
    end
end

% Cria 144 pontos a partir de 96 
FV_s = zeros(N,144); k = 1; l = 1; ll = 0;
for mc=1:N
    for i=2:2:96
       FV_s(mc,i+k) = Ppv_s(mc,i); 
       FV_w(mc,i+k) = Ppv_w(mc,i); k = k + 1;
       FV_s(mc,i+ll) = ( 0.5*Ppv_s(mc,i)+Ppv_s(mc,i-1) )/1.5; 
       FV_w(mc,i+ll) = ( 0.5*Ppv_w(mc,i)+Ppv_w(mc,i-1) )/1.5; ll = ll + 1;
    end
    for j=3:2:96
       FV_s(mc,j+l) = ( 2*Ppv_s(mc,j)+Ppv_s(mc,j-1) )/3 ; 
       FV_w(mc,j+l) = ( 2*Ppv_w(mc,j)+Ppv_w(mc,j-1) )/3; l = l + 1; 
    end
    k = 1; l = 1; ll = 0;
end

% Converte pot. gerada de kW para MW
FV_s = FV_s./1000 ; FV_w = FV_w./1000;

% Salva resultado da pot. gerada pela FV em MW
% save FVcurves1000pts.mat FV_s FV_w; % retorna FV_s FV_w

%% Parameters preparation
verbose = 1;
i = 10;
Cep6data = importdata('Copy_2_of_thp6_strip10.0001.txt');
Ce3data = importdata('Copy_2_of_th3_strip10.0001.txt');
Ep6 = Cep6data(:,1);Np6 = Cep6data(:,2);
E3 = Ce3data(:,1);N3 = Ce3data(:,2);
N_SAp6=[];
N_SA3=[];
X = [0,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.425,0.45];
consp6 = [0.4825,0.4830,0.4836,0.4843,0.4851,0.4859,0.4869,0.4879,0.4890,0.4902,0.4913,0.4926,0.4938,0.4951,0.4964,0.4977,0.4991,0.5004,0.5017];
cons3 = [0.6966,0.6967,0.6969,0.6972,0.6976,0.6981,0.6987,0.6994,0.7001,0.7009,0.7017,0.7026,0.7035,0.7045,0.7055,0.7064,0.7075,0.7085,0.7095];
A = 0.54*X(i)/(1-X(i));
B = (0.395*X(i)+0.518)/(1-X(i));
% sin(0.6) = 0.01047,sin(3) = 0.05234
Cp6 = (8.73*X(i)^2+8.79*X(i)+17.356)/0.01047*(1-X(i))*0.001/2.53;
C3 = (67.101*X(i)^2+67.1*X(i)+153.73)/0.05234*(1-X(i))*0.001/2.53;
%% 0.6 correction
for j = 1:284
if (verbose)
syms V
solx = vpasolve(V*(1+A+B*0.01047)/(V+A+B*0.01047)*(1-exp(-Cp6*(V+A)))/consp6(i) == Np6(j),V);
y = double(solx);
N_SAp6 = [N_SAp6; Ep6(j) Np6(j) y ];
j = j + 1;
end
end
%% 3 correction
for j = 1:248
if (verbose)
syms V2
solx = vpasolve(V2*(1+A+B*0.05234)/(V2+A+B*0.05234)*(1-exp(-C3*(V2+A)))/cons3(i) == N3(j),V2);
z = double(solx);
N_SA3 = [N_SA3; E3(j) N3(j) z ];
j = j + 1;
end
end
%% Plot the data
if(verbose)
    figure;
    plot(N_SAp6(:,1),N_SAp6(:,2), 'r-', ...
         N_SAp6(:,1),N_SAp6(:,3), 'b-',...
         Ce3data(:,1),Ce3data(:,2),'r--',...
         N_SA3(:,1),N_SA3(:,3),'b--');   
    xlabel('Energy (eV)');
    ylabel('\mu / absorption coefficient');
    legend('0.6', '0.6 correction','3', '3 correction');
    title('SA correction');
end





%% Input parameters
awCe=140.115,awZr=91.224,awO=16,dop=0.45;% molecular weight (g/cm^3)
zCe=58,zZr=40,zO=8;%atomic number
aCeO2=5.412,aZ125=5.381,aZ25=5.344,aZ45=5.292;%lattice costant (A)
Energy=5.733; %incident X-ray energy(Kev).
f1Ce=34.87;f2Ce=11.12;f1Zr=40.16,f2Zr=4.038;f1O=8.09;f2O=0.06884;%fp annd fpp
%% constants
r0=2.82*10^(-5);%classical electron radius,(A)
Na=0.6022;%Na:avogadro (A^3/cm^3/mol)
%% mass density
awCeO2= awCe*(1-dop)+awZr*dop+2*awO;
rhoCeO2=4*awCeO2/(Na*aZ125^3);%mass desnity (g/cm^3)
zCeO2=zCe*(1-dop)+zZr*dop+2*zO;
%% fp and fpp
fpCe=f1Ce-zCe;fppCe=f2Ce;fpZr=f1Zr-zZr;fppZr=f2Zr;fpO=f1O-zO;fppO=f2O;
fpCeO2=fpCe*(1-dop)+fpZr*dop+2*fpO
fppCeO2=fppCe*(1-dop)+fppZr*dop+2*fppO;
%% eletron density
ReNe=rhoCeO2*(Na/awCeO2).*(fpCeO2+zCeO2)
ImNe=rhoCeO2*(Na/awCeO2).*(fppCeO2);
%% critical angle
lamda=12.389./Energy;% wavelength of the X-ray(A).
delta=ReNe*r0.*lamda.^2/(2*pi);% n=1-delta+i*beta.
beta=ImNe*r0.*lamda.^2/(2*pi);
alphac=sqrt(delta*2);% critical angle of CeO2(radian).
alphac_deg=alphac/pi*180
%% penetration depth
x=0:0.01:9;
b=fppCeO2./(fpCeO2+zCeO2);
qc=4*pi.*sin(alphac)./lamda;
mue1=qc(1)./sqrt(2).*sqrt(sqrt(((x./alphac_deg(1)).^2-1).^2+b(1)^2)-((x./alphac_deg(1)).^2-1));
depth1=1./mue1/10; % penetration depth (nm),below E0.
p1=plot(x,depth1,'k');
xlabel('incident angle(degree)'),ylabel('penetration depth(nm)');
x=x';
depth1=depth1';
data=[x,depth1];
save('pd_E.txt','data','-ascii');

%% Varying permeate flowrate 
clear
clc
C14D=1/88.11; %mM  
CNOM=5/12; %mM as C
z_NOM=-3.296666667;

%------------------------------OH----------------
% k_SO4_NOM=23.8*10^7; %OH
% k_SO4_14D=2.8*10^9; %Emily Marron
%------------------------------OH----------------

%------------------------------SO4----------------
% k_SO4_NOM=1.66*10^7;
% k_SO4_14D=7.3*10^7; % Rate constants for hydrogen abstraction reactions of the sulfate radical, SO4 −: Experimental and theoretical results for cyclic ethers. The Journal of Physical Chemistry. 1991, 95 (23), 9336-9340.
%------------------------------SO4----------------

%------------------------------1O2----------------
k_SO4_NOM=2.65*10^5;   %ENVIRONMENTAL SCIENCE & TECHNOLOGY / VOL. 43, NO. 3, 2009
k_SO4_14D=10^3.08; % Rate constants for hydrogen abstraction reactions of the sulfate radical, SO4 −: Experimental and theoretical results for cyclic ethers. The Journal of Physical Chemistry. 1991, 95 (23), 9336-9340.
%------------------------------1O2----------------

Radic_yield = 0.50;

Vf=[10	20	30	40	50	60	80	100 125	150]/3.6*10^-6;
% Vf=[10	20	30	40	50	60	80	100 125	150]/3.6*10^-6;
%% Membrane
e=0.06; % effective porosity
Lm=120*10^-9; %150 membrane thickness (m)
X=-50; %fix charge (mM)
Pd=1; % pore size (nm)

Plate_num = 43;
del_x = Lm/(Plate_num-1);

%% Ions
d_14D=0.234*2; %nm
D_14D=9.1*10^-10; % diffusivity (m2/s) 
ld_14D=d_14D/Pd;  %lambda = ri/rp

d_NOM=0.963*2; %nm
D_NOM=2.55*10^-10; % diffusivity (m2/s) 
ld_NOM=d_NOM/Pd;  %lambda = ri/rp
% Nonelectrostatic hindrance coefficient
Kne_14D=(1+9/8*ld_14D*log(ld_14D)-1.56034*ld_14D+0.528155*ld_14D^2+1.91521*ld_14D^3-2.81903*ld_14D^4+0.270788*ld_14D^5+1.10115*ld_14D^6-0.435933*ld_14D^7)/(1-ld_14D)^2;%  (1+9/8 ld ln ld -1.56 ld +0.53 ld^2 + 1.95 ld^3 - 2.82 ld^4 +0.27 ld^5 +1.10 ld^6 -0.44 ld^7)/(1-ld)^2

%%Frictional factor 
Kf_14D=Kne_14D; %Ref Ruoyu 2024 EST Li/Mg separation

%% Partitioning coefficient
ps_14D=(1-ld_14D)^2; % Steric 

%% Hydrodynamics
% % Delta_f=1.29904E-05; % 0.1: 9.41559E-05  0.4: 5.95892E-05  0.7 2.85478E-05  1.0 2.18472E-05  2: 1.29904E-05  5:6.53384E-06
D_PMS=1.15*10^-9;
Delta_f=D_PMS*3.6*10^4;
k_mass_trans_14D=D_14D/Delta_f; %D/sdl; % mass transfer coefficient (m/s)  100 LMH=100/1000/3600=1/3.6*10^-4 m/s kf=D/delta_f
k_mass_trans_NOM=D_NOM/Delta_f;
%% Import XI values XI_record 0513 transport mechanism bottom table
XI_record=[1304.121022	2608.242045	3912.363067	4814.717187	5331.395937	5852.064538	6887.194258	7912.586487	9177.866822	10421.32196];
% XI_record=[5331.395937];

Surf_rate_record=[8.22131E-06	8.52E-06	8.82263E-06	9.26E-06	9.85E-06	1.05E-05	1.18E-05	1.32E-05	1.51E-05	1.72E-05];
% Surf_rate_record=[9.85E-06];

% In_membrane_rate_record=zeros(10,42);
In_membrane_rate_record=[1.668E-08	1.611E-08	1.556E-08	1.502E-08	1.449E-08	1.398E-08	1.348E-08	1.299E-08	1.251E-08	1.204E-08	1.158E-08	1.113E-08	1.069E-08	1.027E-08	9.845E-09	9.433E-09	9.030E-09	8.634E-09	8.246E-09	7.865E-09	7.490E-09	7.123E-09	6.761E-09	6.405E-09	6.055E-09	5.711E-09	5.371E-09	5.036E-09	4.705E-09	4.378E-09	4.055E-09	3.736E-09	3.420E-09	3.107E-09	2.797E-09	2.489E-09	2.184E-09	1.880E-09	1.578E-09	1.278E-09	9.784E-10	6.800E-10	0.000E+00
1.928E-08	1.862E-08	1.797E-08	1.734E-08	1.672E-08	1.612E-08	1.553E-08	1.496E-08	1.440E-08	1.385E-08	1.331E-08	1.279E-08	1.227E-08	1.177E-08	1.127E-08	1.079E-08	1.032E-08	9.849E-09	9.391E-09	8.941E-09	8.499E-09	8.065E-09	7.637E-09	7.216E-09	6.801E-09	6.392E-09	5.988E-09	5.590E-09	5.196E-09	4.807E-09	4.422E-09	4.041E-09	3.663E-09	3.289E-09	2.917E-09	2.548E-09	2.181E-09	1.816E-09	1.453E-09	1.090E-09	7.291E-10	3.686E-10	0.000E+00
2.228E-08	2.151E-08	2.076E-08	2.004E-08	1.932E-08	1.863E-08	1.795E-08	1.728E-08	1.663E-08	1.600E-08	1.538E-08	1.477E-08	1.417E-08	1.359E-08	1.302E-08	1.245E-08	1.190E-08	1.136E-08	1.083E-08	1.031E-08	9.797E-09	9.292E-09	8.795E-09	8.306E-09	7.824E-09	7.348E-09	6.879E-09	6.415E-09	5.957E-09	5.504E-09	5.056E-09	4.612E-09	4.172E-09	3.735E-09	3.302E-09	2.871E-09	2.443E-09	2.017E-09	1.592E-09	1.169E-09	7.464E-10	3.248E-10	0.000E+00
2.569E-08	2.481E-08	2.395E-08	2.311E-08	2.229E-08	2.149E-08	2.070E-08	1.994E-08	1.919E-08	1.846E-08	1.774E-08	1.704E-08	1.635E-08	1.568E-08	1.501E-08	1.437E-08	1.373E-08	1.311E-08	1.249E-08	1.189E-08	1.130E-08	1.072E-08	1.014E-08	9.576E-09	9.019E-09	8.469E-09	7.926E-09	7.390E-09	6.860E-09	6.336E-09	5.817E-09	5.302E-09	4.793E-09	4.287E-09	3.784E-09	3.285E-09	2.789E-09	2.294E-09	1.802E-09	1.311E-09	8.206E-10	3.311E-10	0.000E+00
2.955E-08	2.854E-08	2.755E-08	2.659E-08	2.565E-08	2.473E-08	2.383E-08	2.295E-08	2.209E-08	2.125E-08	2.042E-08	1.961E-08	1.882E-08	1.805E-08	1.729E-08	1.654E-08	1.581E-08	1.509E-08	1.439E-08	1.369E-08	1.301E-08	1.234E-08	1.168E-08	1.103E-08	1.039E-08	9.752E-09	9.126E-09	8.508E-09	7.897E-09	7.292E-09	6.693E-09	6.100E-09	5.511E-09	4.928E-09	4.348E-09	3.771E-09	3.198E-09	2.626E-09	2.057E-09	1.490E-09	9.231E-10	3.570E-10	0.000E+00
3.396E-08	3.280E-08	3.167E-08	3.056E-08	2.948E-08	2.843E-08	2.739E-08	2.639E-08	2.540E-08	2.443E-08	2.349E-08	2.256E-08	2.165E-08	2.076E-08	1.989E-08	1.903E-08	1.819E-08	1.737E-08	1.655E-08	1.576E-08	1.497E-08	1.420E-08	1.344E-08	1.269E-08	1.195E-08	1.122E-08	1.050E-08	9.791E-09	9.087E-09	8.391E-09	7.701E-09	7.018E-09	6.340E-09	5.667E-09	4.999E-09	4.334E-09	3.673E-09	3.014E-09	2.357E-09	1.702E-09	1.049E-09	3.952E-10	0.000E+00
4.467E-08	4.316E-08	4.167E-08	4.022E-08	3.881E-08	3.742E-08	3.607E-08	3.475E-08	3.345E-08	3.218E-08	3.094E-08	2.973E-08	2.853E-08	2.737E-08	2.622E-08	2.509E-08	2.399E-08	2.290E-08	2.184E-08	2.079E-08	1.976E-08	1.874E-08	1.774E-08	1.675E-08	1.578E-08	1.482E-08	1.387E-08	1.293E-08	1.200E-08	1.108E-08	1.017E-08	9.267E-09	8.371E-09	7.482E-09	6.598E-09	5.719E-09	4.843E-09	3.971E-09	3.102E-09	2.234E-09	1.368E-09	5.013E-10	0.000E+00
5.843E-08	5.645E-08	5.452E-08	5.263E-08	5.079E-08	4.899E-08	4.722E-08	4.550E-08	4.381E-08	4.216E-08	4.054E-08	3.895E-08	3.740E-08	3.587E-08	3.437E-08	3.290E-08	3.146E-08	3.004E-08	2.865E-08	2.728E-08	2.593E-08	2.460E-08	2.328E-08	2.199E-08	2.072E-08	1.946E-08	1.821E-08	1.698E-08	1.576E-08	1.456E-08	1.336E-08	1.218E-08	1.100E-08	9.833E-09	8.671E-09	7.515E-09	6.364E-09	5.217E-09	4.072E-09	2.930E-09	1.788E-09	6.469E-10	0.000E+00
8.095E-08	7.823E-08	7.557E-08	7.297E-08	7.043E-08	6.795E-08	6.551E-08	6.313E-08	6.080E-08	5.852E-08	5.629E-08	5.410E-08	5.195E-08	4.984E-08	4.777E-08	4.574E-08	4.374E-08	4.178E-08	3.985E-08	3.795E-08	3.607E-08	3.423E-08	3.241E-08	3.062E-08	2.885E-08	2.710E-08	2.537E-08	2.366E-08	2.197E-08	2.029E-08	1.863E-08	1.698E-08	1.534E-08	1.372E-08	1.210E-08	1.049E-08	8.880E-09	7.279E-09	5.681E-09	4.085E-09	2.490E-09	8.942E-10	0.000E+00
1.108E-07	1.071E-07	1.035E-07	9.997E-08	9.651E-08	9.313E-08	8.982E-08	8.657E-08	8.340E-08	8.028E-08	7.723E-08	7.424E-08	7.131E-08	6.843E-08	6.560E-08	6.283E-08	6.010E-08	5.741E-08	5.477E-08	5.217E-08	4.961E-08	4.708E-08	4.459E-08	4.214E-08	3.971E-08	3.731E-08	3.494E-08	3.259E-08	3.027E-08	2.796E-08	2.568E-08	2.341E-08	2.116E-08	1.891E-08	1.669E-08	1.446E-08	1.225E-08	1.004E-08	7.840E-09	5.638E-09	3.436E-09	1.231E-09	0.000E+00];

% In_membrane_rate_record=[2.955E-08	2.854E-08	2.755E-08	2.659E-08	2.565E-08	2.473E-08	2.383E-08	2.295E-08	2.209E-08	2.125E-08	2.042E-08	1.961E-08	1.882E-08	1.805E-08	1.729E-08	1.654E-08	1.581E-08	1.509E-08	1.439E-08	1.369E-08	1.301E-08	1.234E-08	1.168E-08	1.103E-08	1.039E-08	9.752E-09	9.126E-09	8.508E-09	7.897E-09	7.292E-09	6.693E-09	6.100E-09	5.511E-09	4.928E-09	4.348E-09	3.771E-09	3.198E-09	2.626E-09	2.057E-09	1.490E-09	9.231E-10	3.570E-10	0.000E+00];

Error_14D=zeros(1,length(Vf));
Error_NOM_14D_min=zeros(1,length(Vf));

C_m_NOM_record=zeros(1,length(Vf));
C_m_14D_ppm=zeros(1,length(Vf));
C_m_14D_min_ppm=zeros(1,length(Vf));
c_p_14D=zeros(1,length(Vf));
c_p_14D_min=zeros(1,length(Vf));

RXN_rate_14D_Surf=zeros(1,length(Vf));
RXN_rate_14D_In_membrane=zeros(1,length(Vf));

for i_vf =1:length(Vf)
vf=Vf(i_vf);
XI=XI_record(i_vf);
Surf_rate=Surf_rate_record(i_vf);
In_membrane_rate=In_membrane_rate_record(i_vf,:);

Iteration_limit=5000;
%% Solver loop
C_14D_guss=zeros(1,Iteration_limit);
C_14D_min_guss=zeros(1,Iteration_limit);
C_NOM_guss=zeros(1,Iteration_limit);

Branch_14D_min=zeros(1,Iteration_limit);
Branch_NOM=zeros(1,Iteration_limit);

C_NOM_guss(1)=CNOM*0.01;
% C_14D_min_guss(1)=C14D*0.01; % initial guess/same as feed
C_14D_min_guss(1)=C14D; % initial guess/same as feed

for i=1:length(C_14D_min_guss)
C_14D_min0=C_14D_min_guss(i)*ps_14D;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
a = sym('a', [1, Plate_num]);
    eq(1)=a(1)-C_14D_min0;
    for l=1:Plate_num-2
    eq(l+1)=-D_14D*Kf_14D*e*(a(l+1)-a(l))/del_x+Kf_14D*0.5*(a(l+1)+a(l))*vf-(-D_14D*Kf_14D*e*(a(l+2)-a(l+1))/del_x+Kf_14D*0.5*(a(l+2)+a(l+1))*vf)-In_membrane_rate(l)*Radic_yield/20;
    end
    eq(Plate_num)=-D_14D*Kf_14D*e*(a(Plate_num)-a(Plate_num-1))/del_x+Kf_14D*0.5*(a(Plate_num)+a(Plate_num-1))*vf-vf*a(Plate_num)/ps_14D-In_membrane_rate(Plate_num-1)*Radic_yield/20;

    if i==1
    sol=vpasolve(eq,a);
    Init_guss=struct2cell(sol);
    else
    sol=vpasolve(eq,a,Init_guss);
    Init_guss=struct2cell(sol);
    end
    disp(sol)
for m=1:Plate_num
    C_14D_min(m)= double([sol.(['a', num2str(m)])]);
end

cp_14D_min=C_14D_min(Plate_num)/ps_14D;
if cp_14D_min<0
   cp_14D_min=0;
end

Total_14D_min_removal=0;
Index_for_non_zero_14D_min=max(find(C_14D_min > 0));

    for i_14D_min =1:Index_for_non_zero_14D_min-1
    Total_14D_min_removal=Total_14D_min_removal+In_membrane_rate(i_14D_min)*Radic_yield/20;
    end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% CP solver for membrane surface conc
syms Cal_14D Cal_14D_min Cal_NOM % a-Na b-Cl, c-potential diff, d-donnan diff
Branch_14D_min(i)=C_14D_min_guss(i)*k_SO4_14D/(C_14D_min_guss(i)*k_SO4_14D+C_NOM_guss(i)*k_SO4_NOM);
Branch_NOM(i)=C_NOM_guss(i)*k_SO4_NOM/(C_14D_min_guss(i)*k_SO4_14D+C_NOM_guss(i)*k_SO4_NOM);

eqn1=((1+z_NOM*D_NOM*XI/vf)*Cal_NOM-Surf_rate*Radic_yield*Branch_NOM(i)/4/vf)/((1+z_NOM*D_NOM*XI/vf)*CNOM-Surf_rate*Radic_yield*Branch_NOM(i)/4/vf)-exp(vf/k_mass_trans_NOM*(1+z_NOM*D_NOM*XI/vf)); %NOM
eqn2=(Cal_14D_min-cp_14D_min-Total_14D_min_removal/vf-Surf_rate*Radic_yield*Branch_14D_min(i)/20/vf)/(C14D-cp_14D_min-Total_14D_min_removal/vf-Surf_rate*Radic_yield*Branch_14D_min(i)/20/vf)-exp(vf/k_mass_trans_14D);

    [aa,bb]=vpasolve(eqn1,eqn2,Cal_14D_min,Cal_NOM);
    
    if isempty([aa,bb])==1
        disp('NOM, and 14D_min wrong')
        Error_NOM_14D_min(i_vf)=1;
        break
    end
    C_14D_min_cal=aa;
    C_NOM_cal=bb;
    
        if (abs(C_NOM_guss(i)-C_NOM_cal)/C_NOM_cal<1E-3) && (abs(C_14D_min_guss(i)-C_14D_min_cal)/C_14D_min_cal<1E-3)
             break
        else
           C_14D_min_guss(i+1)=0.5*C_14D_min_cal+0.5*C_14D_min_guss(i);
           C_NOM_guss(i+1)=0.5*C_NOM_cal+0.5*C_NOM_guss(i);
        end  
disp(i)
end

C_m_NOM=C_NOM_cal;
C_m_14D_min=C_14D_min_cal;

% C_14D_guss(1)=0.000001*C_m_14D_min; % initial guess/same as feed
 C_14D_guss(1)=0.00001; % initial guess/same as feed
for i=1:length(C_14D_guss)
C_14D0=C_14D_guss(i)*ps_14D;
b = sym('b', [1, Plate_num]);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~
    eq(1)=b(1)-C_14D0;
    for l=1:Plate_num-2
    eq(l+1)=-D_14D*Kf_14D*e*(b(l+1)-b(l))/del_x+Kf_14D*0.5*(b(l+1)+b(l))*vf-(-D_14D*Kf_14D*e*(b(l+2)-b(l+1))/del_x+Kf_14D*0.5*(b(l+2)+b(l+1))*vf)-In_membrane_rate(l)*Radic_yield*0.5*(b(l+1)+b(l))/(0.5*(C_14D_min(l+1)+C_14D_min(l)));
    end
    eq(Plate_num)=-D_14D*Kf_14D*e*(b(Plate_num)-b(Plate_num-1))/del_x+Kf_14D*0.5*(b(Plate_num)+b(Plate_num-1))*vf-vf*b(Plate_num)/ps_14D-In_membrane_rate(Plate_num-1)*Radic_yield*0.5*(b(Plate_num)+b(Plate_num-1))/(0.5*(C_14D_min(Plate_num)+C_14D_min(Plate_num-1)));

    if i==1
    sol=vpasolve(eq,b);
    Init_guss=struct2cell(sol);
    else
    sol=vpasolve(eq,b,Init_guss);
    Init_guss=struct2cell(sol);
    end
    disp(sol)
for m=1:Plate_num
    C_14D(m)= double([sol.(['b', num2str(m)])]);
end

cp_14D=C_14D(Plate_num)/ps_14D;
if cp_14D<0
   cp_14D=0;
end

Total_14D_removal=0;
Index_for_non_zero_14D=max(find(C_14D > 0));

    for i_14D =1:Index_for_non_zero_14D-1
    Total_14D_removal= Total_14D_removal+In_membrane_rate(i_14D)*Radic_yield*0.5*(C_14D(i_14D)+C_14D(i_14D+1))/(0.5*(C_14D_min(i_14D)+C_14D_min(i_14D+1)));
    end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% CP solver for membrane surface conc
Branch_14D=Cal_14D*k_SO4_14D/(C_m_14D_min*k_SO4_14D+C_m_NOM*k_SO4_NOM);
eqn2=(Cal_14D-cp_14D-Surf_rate*Radic_yield*Branch_14D/vf-Total_14D_removal/vf)/(C14D-cp_14D-Surf_rate*Radic_yield*Branch_14D/vf-Total_14D_removal/vf)-exp(vf/k_mass_trans_14D); %14D

% [cc]=vpasolve(eqn2,Cal_14D,[0,100*C14D]);
[cc]=vpasolve(eqn2,Cal_14D,[0,inf]);
    if isempty(cc)==1
        disp('14D wrong')
        Error_14D(i_vf)=1;
        break
    end
    C_14D_cal=cc;
if abs(C_14D_guss(i)-C_14D_cal)/C_14D_cal<1E-3 
   break
else
C_14D_guss(i+1)=0.9*C_14D_guss(i)+0.1*C_14D_cal;
% C_14D_guss(i+1)=0.98*C_14D_guss(i)+0.02*C_14D_cal;
% C_14D_guss(i+1)=0.5*C_14D_guss(i)+0.5*C_14D_cal;
end
disp(['14D = ',num2str(i)]) 
end

C_m_14D=C_14D_cal;
c_p_14D=cp_14D;
c_p_14D_min=cp_14D_min;
% end
% end

% C_m_NOM0=C_m_NOM;
% C_m_14D_min0=C_m_14D_min;
% C_m_NOM0(C_m_NOM0<0)=0;
% C_m_14D_min0(C_m_14D_min0<0)=0;
C_m_NOM_record(i_vf)=C_m_NOM;
C_m_14D_ppm(i_vf)=C_m_14D*88.11;
C_m_14D_min_ppm(i_vf)=C_m_14D_min*88.11;

C_14D_ppm(i_vf,:)=C_14D*88.11;
C_14D_min_ppm(i_vf,:)=C_14D_min*88.11;



c_p_14D_ppm(i_vf)=c_p_14D*88.11;
c_p_14D_min_ppm(i_vf)=c_p_14D_min*88.11;
disp('NOM')
disp(C_m_NOM)
disp('14D')
disp(C_m_14D)

RXN_rate_14D_Surf(i_vf)=Surf_rate*Radic_yield*C_m_14D*k_SO4_14D/(C_m_14D_min*k_SO4_14D+C_m_NOM*k_SO4_NOM);
RXN_rate_14D_In_membrane(i_vf)=Total_14D_removal;


end
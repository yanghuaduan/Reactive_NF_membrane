%% Varying permeate flowrate 
clear
clc
C14D=1/88.11; %mM  
NOM_small_fraction=0.05;
CNOM_large=5*(1-NOM_small_fraction)/12; %mM as C  Split the NOM into two sections
CNOM_small=5*NOM_small_fraction/12; %mM as C
z_NOM=-3.296666667;


Initial_guess_A=ones(1,43)*5/12;
Initial_guess_B=[6.77E-09	6.61E-09	6.46E-09	6.30E-09	6.14E-09	5.98E-09	5.83E-09	5.67E-09	5.51E-09	5.35E-09	5.20E-09	5.04E-09	4.88E-09	4.72E-09	4.56E-09	4.40E-09	4.24E-09	4.08E-09	3.92E-09	3.76E-09	3.60E-09	3.44E-09	3.28E-09	3.11E-09	2.95E-09	2.79E-09	2.63E-09	2.47E-09	2.30E-09	2.14E-09	1.98E-09	1.81E-09	1.65E-09	1.49E-09	1.32E-09	1.16E-09	9.93E-10	8.28E-10	6.63E-10	4.98E-10	3.32E-10	1.66E-10	3.66E-18];
Init_guss=[Initial_guess_A,Initial_guess_B];


k_SO4_NOM=1.66*10^7;
k_SO4_14D=7.3*10^7; % Rate constants for hydrogen abstraction reactions of the sulfate radical, SO4 −: Experimental and theoretical results for cyclic ethers. The Journal of Physical Chemistry. 1991, 95 (23), 9336-9340.

Radic_yield = 0.9;

Vf=[20]/3.6*10^-6;
% Vf=[125	150]/3.6*10^-6;
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

d_NOM_large=0.963*2; %nm
d_NOM_small=0.25*2; %nm

D_NOM=2.55*10^-10; % diffusivity (m2/s) 
ld_NOM_small=d_NOM_small/Pd;  %lambda = ri/rp
% Nonelectrostatic hindrance coefficient
Kne_14D=(1+9/8*ld_14D*log(ld_14D)-1.56034*ld_14D+0.528155*ld_14D^2+1.91521*ld_14D^3-2.81903*ld_14D^4+0.270788*ld_14D^5+1.10115*ld_14D^6-0.435933*ld_14D^7)/(1-ld_14D)^2;%  (1+9/8 ld ln ld -1.56 ld +0.53 ld^2 + 1.95 ld^3 - 2.82 ld^4 +0.27 ld^5 +1.10 ld^6 -0.44 ld^7)/(1-ld)^2
Kne_NOM_small=(1+9/8*ld_NOM_small*log(ld_NOM_small)-1.56034*ld_NOM_small+0.528155*ld_NOM_small^2+1.91521*ld_NOM_small^3-2.81903*ld_NOM_small^4+0.270788*ld_NOM_small^5+1.10115*ld_NOM_small^6-0.435933*ld_NOM_small^7)/(1-ld_NOM_small)^2;
%%Frictional factor 
Kf_14D=Kne_14D; %Ref Ruoyu 2024 EST Li/Mg separation
Kf_NOM_small=Kne_NOM_small;

%% Partitioning coefficient
ps_14D=(1-ld_14D)^2; % Steric 
ps_NOM_small=(1-ld_NOM_small)^2;

%% Hydrodynamics
% % Delta_f=1.29904E-05; % 0.1: 9.41559E-05  0.4: 5.95892E-05  0.7 2.85478E-05  1.0 2.18472E-05  2: 1.29904E-05  5:6.53384E-06
D_PMS=1.15*10^-9;
Delta_f=D_PMS*3.6*10^4;
k_mass_trans_14D=D_14D/Delta_f; %D/sdl; % mass transfer coefficient (m/s)  100 LMH=100/1000/3600=1/3.6*10^-4 m/s kf=D/delta_f
k_mass_trans_NOM=D_NOM/Delta_f;
%% Import XI values XI_record 0513 transport mechanism bottom table
Surf_rate_record=[0];

XI_record=[1247.773487];

Phi=[-4.332307184	-4.332175245	-4.332043003	-4.331910472	-4.331777667	-4.3316446	-4.331511285	-4.331377734	-4.331243958	-4.331109968	-4.330975775	-4.330841389	-4.330706819	-4.330572076	-4.330437167	-4.330302102	-4.330166888	-4.330031534	-4.329896047	-4.329760433	-4.329624701	-4.329488856	-4.329352905	-4.329216855	-4.32908071	-4.328944478	-4.328808163	-4.328671772	-4.328535308	-4.328398777	-4.328262185	-4.328125535	-4.327988833	-4.327852082	-4.327715288	-4.327578455	-4.327441586	-4.327304687	-4.327167761	-4.327030812	-4.326893844	-4.326756861	-4.326619868];

% Phi=zeros(10,42);
d0=Phi(1);
d_permeate=-6.137692785;

% In_membrane_rate_record=zeros(10,42);
In_membrane_rate_record=[6.65E-08	6.33E-08	6.02E-08	5.73E-08	5.45E-08	5.18E-08	4.92E-08	4.68E-08	4.44E-08	4.22E-08	4.00E-08	3.80E-08	3.60E-08	3.41E-08	3.23E-08	3.06E-08	2.89E-08	2.73E-08	2.58E-08	2.43E-08	2.29E-08	2.15E-08	2.02E-08	1.90E-08	1.77E-08	1.65E-08	1.54E-08	1.43E-08	1.32E-08	1.21E-08	1.11E-08	1.01E-08	9.16E-09	8.21E-09	7.27E-09	6.36E-09	5.46E-09	4.57E-09	3.70E-09	2.83E-09	1.97E-09	1.12E-09];

Error_14D=zeros(1,length(Vf));
Error_NOM_14D_min=zeros(1,length(Vf));

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
C_NOM_large_guss=zeros(1,Iteration_limit);
C_NOM_small_guss=zeros(1,Iteration_limit);

Branch_14D_min=zeros(1,Iteration_limit);
Branch_NOM_large=zeros(1,Iteration_limit);
Branch_NOM_small=zeros(1,Iteration_limit);

C_NOM_large_guss(1)=CNOM_large;
C_NOM_small_guss(1)=CNOM_small;
C_14D_min_guss(1)=0.011503065; % initial guess/same as feed

for i=1:length(C_14D_min_guss)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C_14D_min0=C_14D_min_guss(i)*ps_14D;
C_NOM_small0=C_NOM_small_guss(i)*ps_NOM_small*exp(-z_NOM*d0);

idx = 1:Plate_num;
idx_nom = Plate_num+1:2*Plate_num;

% Build initial guess vector
Init_vec = Init_guss(1:2*Plate_num);  % Assuming Init_guss already has numeric guesses

% Anonymous function to represent the 2*Plate_num equations
eqfun = @(x) build_equations(x, Plate_num, del_x, D_14D, Kf_14D, e, vf, ...
    In_membrane_rate, Radic_yield, k_SO4_14D, k_SO4_NOM, Phi, ...
    D_NOM, Kf_NOM_small, z_NOM, ps_14D, ps_NOM_small, d0, d_permeate, ...
    C_14D_min0, C_NOM_small0);

% Solve using fsolve
options = optimoptions('fsolve','Display','none','FunctionTolerance',1e-10,'MaxIterations',1000);
[x_sol, fval, exitflag] = fsolve(eqfun, Init_vec, options);

if exitflag <= 0
    warning('fsolve did not converge');
    break
end

% Extract a and nom values
C_14D_min = x_sol(idx);
C_NOM_small = x_sol(idx_nom);

% a = sym('a', [1, Plate_num]);
% nom = sym('nom', [1, Plate_num]);
% 
%     eq(1)=a(1)-C_14D_min0;
%     eq(2)=nom(1)-C_NOM_small0;
% 
% for l=1:Plate_num-2
%     eq(l*2+1)=-D_14D*Kf_14D*e*(a(l+1)-a(l))/del_x+Kf_14D*0.5*(a(l+1)+a(l))*vf-(-D_14D*Kf_14D*e*(a(l+2)-a(l+1))/del_x+Kf_14D*0.5*(a(l+2)+a(l+1))*vf)-In_membrane_rate(l)*Radic_yield/20*k_SO4_14D*a(l+1)/(k_SO4_14D*a(l+1)+k_SO4_NOM*nom(l+1));
%     eq(l*2+2)=-D_NOM*Kf_NOM_small*e*(nom(l+1)-nom(l))/del_x+Kf_NOM_small*0.5*(nom(l+1)+nom(l))*vf-z_NOM*D_NOM*Kf_NOM_small*e*0.5*(nom(l+1)+nom(l))*(Phi(l+1)-Phi(l))/del_x-(-D_NOM*Kf_NOM_small*e*(nom(l+2)-nom(l+1))/del_x+Kf_NOM_small*0.5*(nom(l+2)+nom(l+1))*vf-z_NOM*D_NOM*Kf_NOM_small*e*0.5*(nom(l+2)+nom(l+1))*(Phi(l+2)-Phi(l+1))/del_x)-In_membrane_rate(l)*Radic_yield/4*k_SO4_NOM*nom(l+1)/(k_SO4_14D*a(l+1)+k_SO4_NOM*nom(l+1));
% end
%     eq(2*Plate_num-1)=-D_14D*Kf_14D*e*(a(Plate_num)-a(Plate_num-1))/del_x+Kf_14D*0.5*(a(Plate_num)+a(Plate_num-1))*vf-vf*a(Plate_num)/ps_14D-In_membrane_rate(Plate_num-1)*Radic_yield/20*k_SO4_14D*a(Plate_num)/(k_SO4_14D*a(Plate_num)+k_SO4_NOM*nom(Plate_num));
%     eq(2*Plate_num)=-D_NOM*Kf_NOM_small*e*(nom(Plate_num)-nom(Plate_num-1))/del_x+Kf_NOM_small*0.5*(nom(Plate_num)+nom(Plate_num-1))*vf-z_NOM*D_NOM*Kf_NOM_small*e*0.5*(nom(Plate_num)+nom(Plate_num-1))*(Phi(Plate_num)-Phi(Plate_num-1))/del_x-vf*nom(Plate_num)/ps_NOM_small/exp(-z_NOM*d_permeate)-In_membrane_rate(Plate_num-1)*Radic_yield/4*k_SO4_NOM*nom(Plate_num)/(k_SO4_14D*a(Plate_num)+k_SO4_NOM*nom(Plate_num));
% 
% 
%     % if i==1
%     % sol=vpasolve(eq,[a,nom]);
%     % Init_guss=struct2cell(sol);
%     % else
%     sol=vpasolve(eq,[a,nom],Init_guss);
%     Init_guss=struct2cell(sol);
%     % end
%     disp(sol)



% for m=1:Plate_num
%     C_14D_min(m)= double([sol.(['a', num2str(m)])]);
%     C_NOM_small(m)= double([sol.(['nom', num2str(m)])]);
% end

cp_14D_min=C_14D_min(Plate_num)/ps_14D;
cp_NOM_small=C_NOM_small(Plate_num)/ps_NOM_small/exp(-z_NOM*d_permeate);

if cp_14D_min<0
   cp_14D_min=0;
end

Total_14D_min_removal=0;
Total_NOM_In_removal=0;
Index_for_non_zero_14D_min=max(find(C_14D_min > 0));

    for i_14D_min =1:Index_for_non_zero_14D_min-1
    Total_14D_min_removal=Total_14D_min_removal+In_membrane_rate(i_14D_min)*Radic_yield/20*k_SO4_14D*C_14D_min(i_14D_min+1)/(k_SO4_14D*C_14D_min(i_14D_min+1)+k_SO4_NOM*C_NOM_small(i_14D_min+1));
    Total_NOM_In_removal=Total_NOM_In_removal+In_membrane_rate(i_14D_min)*Radic_yield/4*k_SO4_NOM*C_NOM_small(i_14D_min+1)/(k_SO4_14D*C_14D_min(i_14D_min+1)+k_SO4_NOM*C_NOM_small(i_14D_min+1));
    end


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% CP solver for membrane surface conc
syms Cal_14D Cal_14D_min Cal_NOM_large Cal_NOM_small % a-Na b-Cl, c-potential diff, d-donnan diff
Branch_14D_min(i)=C_14D_min_guss(i)*k_SO4_14D/(C_14D_min_guss(i)*k_SO4_14D+(C_NOM_large_guss(i)+C_NOM_small_guss(i))*k_SO4_NOM);
Branch_NOM_small(i)=C_NOM_small_guss(i)*k_SO4_NOM/(C_14D_min_guss(i)*k_SO4_14D+(C_NOM_large_guss(i)+C_NOM_small_guss(i))*k_SO4_NOM);
Branch_NOM_large(i)=C_NOM_large_guss(i)*k_SO4_NOM/(C_14D_min_guss(i)*k_SO4_14D+(C_NOM_large_guss(i)+C_NOM_small_guss(i))*k_SO4_NOM);

eqn1=((1+z_NOM*D_NOM*XI/vf)*Cal_NOM_large-Surf_rate*Radic_yield*Branch_NOM_large(i)/4/vf)/((1+z_NOM*D_NOM*XI/vf)*CNOM_large-Surf_rate*Radic_yield*Branch_NOM_large(i)/4/vf)-exp(vf/k_mass_trans_NOM*(1+z_NOM*D_NOM*XI/vf));
eqn2=((1+z_NOM*D_NOM*XI/vf)*Cal_NOM_small-cp_NOM_small-Surf_rate*Radic_yield*Branch_NOM_small(i)/4/vf-Total_NOM_In_removal/vf)/((1+z_NOM*D_NOM*XI/vf)*CNOM_small-cp_NOM_small-Surf_rate*Radic_yield*Branch_NOM_small(i)/4/vf-Total_NOM_In_removal/vf)-exp(vf/k_mass_trans_NOM*(1+z_NOM*D_NOM*XI/vf));
eqn3=(Cal_14D_min-cp_14D_min-Total_14D_min_removal/vf-Surf_rate*Radic_yield*Branch_14D_min(i)/20/vf)/(C14D-cp_14D_min-Total_14D_min_removal/vf-Surf_rate*Radic_yield*Branch_14D_min(i)/20/vf)-exp(vf/k_mass_trans_14D);

    [aa,bb,cc]=vpasolve(eqn1,eqn2,eqn3,Cal_14D_min,Cal_NOM_large,Cal_NOM_small);
    
    if isempty([aa,bb])==1
        disp('NOM, and 14D_min wrong')
        Error_NOM_14D_min(i_vf)=1;
        break
    end
    C_14D_min_cal=aa;
    C_NOM_large_cal=bb;
    C_NOM_small_cal=cc;

        if (abs(C_NOM_large_guss(i)-C_NOM_large_cal)/C_NOM_large_cal<1E-3) && (abs(C_14D_min_guss(i)-C_14D_min_cal)/C_14D_min_cal<1E-3)
             break
        else
           C_14D_min_guss(i+1)=0.5*C_14D_min_cal+0.5*C_14D_min_guss(i);
           C_NOM_large_guss(i+1)=0.5*C_NOM_large_cal+0.5*C_NOM_large_guss(i);
           C_NOM_small_guss(i+1)=0.5*C_NOM_small_cal+0.5*C_NOM_small_guss(i);
        end  
disp(i)
end

C_m_NOM_large=C_NOM_large_cal;
C_m_NOM_small=C_NOM_small_cal;
C_m_14D_min=C_14D_min_cal;


% C_14D_guss(1)=0.000001*C_m_14D_min; % initial guess/same as feed
C_14D_guss(1)=0.0001; % initial guess/same as feed

for i=1:length(C_14D_guss)
C_14D0=C_14D_guss(i)*ps_14D;
b = sym('b', [1, Plate_num]);

    eqn14D(1)=b(1)-C_14D0;

    for l=1:Plate_num-2
    eqn14D(l+1)=-D_14D*Kf_14D*e*(b(l+1)-b(l))/del_x+Kf_14D*0.5*(b(l+1)+b(l))*vf-(-D_14D*Kf_14D*e*(b(l+2)-b(l+1))/del_x+Kf_14D*0.5*(b(l+2)+b(l+1))*vf)-In_membrane_rate(l)*Radic_yield*k_SO4_14D*0.5*(b(l+1)+b(l))/(k_SO4_14D*0.5*(C_14D_min(l+1)+C_14D_min(l))+k_SO4_NOM*0.5*(C_NOM_small(l+1)+C_NOM_small(l)));
    end
    eqn14D(Plate_num)=-D_14D*Kf_14D*e*(b(Plate_num)-b(Plate_num-1))/del_x+Kf_14D*0.5*(b(Plate_num)+b(Plate_num-1))*vf-vf*b(Plate_num)/ps_14D-In_membrane_rate(Plate_num-1)*Radic_yield*k_SO4_14D*0.5*(b(Plate_num)+b(Plate_num-1))/(k_SO4_14D*0.5*(C_14D_min(Plate_num)+C_14D_min(Plate_num-1))+k_SO4_NOM*0.5*(C_NOM_small(Plate_num)+C_NOM_small(Plate_num-1)));

    if i==1
    sol=vpasolve(eqn14D,b);
    Init_guss=struct2cell(sol);
    else
    sol=vpasolve(eqn14D,b,Init_guss);
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
    Total_14D_removal= Total_14D_removal+In_membrane_rate(i_14D)*Radic_yield*0.5*k_SO4_14D*(C_14D(i_14D)+C_14D(i_14D+1))/(k_SO4_14D*0.5*(C_14D_min(i_14D)+C_14D_min(i_14D+1))+k_SO4_NOM*0.5*(C_NOM_small(i_14D)+C_NOM_small(i_14D+1)));
    end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% CP solver for membrane surface conc
Branch_14D=Cal_14D*k_SO4_14D/(C_m_14D_min*k_SO4_14D+(C_m_NOM_small+C_m_NOM_large)*k_SO4_NOM);
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

C_m_14D_ppm(i_vf)=C_m_14D*88.11;
C_m_14D_min_ppm(i_vf)=C_m_14D_min*88.11;

C_14D_ppm(i_vf,:)=C_14D*88.11;
C_14D_min_ppm(i_vf,:)=C_14D_min*88.11;



c_p_14D_ppm(i_vf)=c_p_14D*88.11;
c_p_14D_min_ppm(i_vf)=c_p_14D_min*88.11;
% disp('NOM')
% disp(C_m_NOM)
disp('14D')
disp(C_m_14D)

% RXN_rate_14D_Surf(i_vf)=Surf_rate*Radic_yield*C_m_14D*k_SO4_14D/(C_m_14D_min*k_SO4_14D+C_m_NOM*k_SO4_NOM);
% RXN_rate_14D_In_membrane(i_vf)=Total_14D_removal;


end
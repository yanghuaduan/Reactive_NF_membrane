clear
clc

Init_guss=[46.1687954519554	46.1706645183242	46.1725335846930	46.1744026510617	46.1762717174305	46.1781407837993	46.1800098501680	46.1818789165368	46.1825972073696	46.1833154982025	46.1840337890353	46.1847520798681	46.1854703707009	46.1861886615338	46.1869069523666	46.1873669223857	46.1878268924048	46.1882868624239	46.1887468324430	46.1892068024622	46.1896667724813	46.1901267425004	46.1905247912888	46.1909228400773	46.1913208888657	46.1917189376542	46.1921169864426	46.1925150352311	46.1929130840195	46.1932924453950	46.1936718067706	46.1940511681461	46.1944305295217	46.1948098908972	46.1951892522728	46.1955686136483	46.1959388063625	46.1963089990767	46.1966791917909	46.1970493845051	46.1974195772193	46.1977897699335	46.1981599626477	3.84127762196004	3.84046716264976	3.83965670333948	3.83884624402920	3.83803578471891	3.83722532540863	3.83641486609835	3.83560440678807	3.83466467087668	3.83372493496528	3.83278519905389	3.83184546314250	3.83090572723110	3.82996599131971	3.82902625540832	3.82805903776776	3.82709182012721	3.82612460248665	3.82515738484610	3.82419016720554	3.82322294956499	3.82225573192443	3.82128340786323	3.82031108380203	3.81933875974084	3.81836643567964	3.81739411161844	3.81642178755724	3.81544946349604	3.81447695272699	3.81350444195795	3.81253193118890	3.81155942041986	3.81058690965081	3.80961439888177	3.80864188811272	3.80767027325553	3.80669865839835	3.80572704354117	3.80475542868398	3.80378381382680	3.80281219896961	3.80184058411243	0.00986672010844444	0.00876750795231766	0.00766829579619088	0.00656908364006409	0.00546987148393731	0.00437065932781053	0.00327144717168375	0.00217223501555697	0.00193024471851420	0.00168825442147143	0.00144626412442866	0.00120427382738590	0.000962283530343131	0.000720293233300364	0.000478302936257598	0.000425010166649112	0.000371717397040626	0.000318424627432140	0.000265131857823655	0.000211839088215169	0.000158546318606683	0.000105253548998197	9.35037819607559e-05	8.17540149233145e-05	7.00042478858732e-05	5.82544808484319e-05	4.65047138109906e-05	3.47549467735493e-05	2.30051797361079e-05	2.03857548033135e-05	1.77663298705192e-05	1.51469049377248e-05	1.25274800049304e-05	9.90805507213597e-06	7.28863013934158e-06	4.66920520654719e-06	4.01929637480727e-06	3.36938754306736e-06	2.71947871132744e-06	2.06956987958753e-06	1.41966104784761e-06	7.69752216107698e-07	1.19843384367783e-07	0.000103176903501706	0.00118208651080882	0.00226099611811593	0.00333990572542304	0.00441881533273015	0.00549772494003727	0.00657663454734438	0.00765554415465149	0.00766581676389073	0.00767608937312997	0.00768636198236920	0.00769663459160844	0.00770690720084768	0.00771717981008692	0.00772745241932616	0.00750047499340858	0.00727349756749101	0.00704652014157344	0.00681954271565587	0.00659256528973829	0.00636558786382072	0.00613861043790315	0.00585734768504523	0.00557608493218731	0.00529482217932939	0.00501355942647147	0.00473229667361355	0.00445103392075563	0.00416977116789771	0.00387450618361337	0.00357924119932904	0.00328397621504470	0.00298871123076037	0.00269344624647603	0.00239818126219170	0.00210291627790736	0.00180253016083128	0.00150214404375521	0.00120175792667913	0.000901371809603052	0.000600985692526975	0.000300599575450898	2.13458374820466e-07	7.01535078333620e-06	6.23360045060884e-06	5.45185011788147e-06	4.67009978515411e-06	3.88834945242674e-06	3.10659911969938e-06	2.32484878697202e-06	1.54309845424465e-06	1.37098560023072e-06	1.19887274621679e-06	1.02675989220285e-06	8.54647038188921e-07	6.82534184174988e-07	5.10421330161055e-07	3.38308476147123e-07	3.00415197883759e-07	2.62521919620396e-07	2.24628641357033e-07	1.86735363093669e-07	1.48842084830306e-07	1.10948806566943e-07	7.30555283035796e-08	6.47270394247287e-08	5.63985505458777e-08	4.80700616670267e-08	3.97415727881757e-08	3.14130839093247e-08	2.30845950304738e-08	1.47561061516228e-08	1.29590681683335e-08	1.11620301850442e-08	9.36499220175489e-09	7.56795421846559e-09	5.77091623517630e-09	3.97387825188701e-09	2.17684026859771e-09	1.86586308736947e-09	1.55488590614122e-09	1.24390872491298e-09	9.32931543684734e-10	6.21954362456489e-10	3.10977181228245e-10	0	-3.65814535965028	-3.65785171680517	-3.65755807396007	-3.65726443111496	-3.65697078826985	-3.65667714542474	-3.65638350257964	-3.65608985973453	-3.65577129132371	-3.65545272291289	-3.65513415450207	-3.65481558609125	-3.65449701768043	-3.65417844926962	-3.65385988085880	-3.65353571801365	-3.65321155516850	-3.65288739232335	-3.65256322947820	-3.65223906663305	-3.65191490378790	-3.65159074094275	-3.65126523670220	-3.65093973246165	-3.65061422822109	-3.65028872398054	-3.64996321973999	-3.64963771549943	-3.64931221125888	-3.64898630177560	-3.64866039229233	-3.64833448280905	-3.64800857332577	-3.64768266384250	-3.64735675435922	-3.64703084487594	-3.64670473626872	-3.64637862766150	-3.64605251905427	-3.64572641044705	-3.64540030183983	-3.64507419323260	-3.64474808462538	-5.96409287925331];

% % Init_guss=[0.102831650361419	0.100118373623197	0.0974050968849745	0.0946918201467520	0.0919785434085295	0.0892652666703070	0.0865519899320845	0.0838387131938621	0.0685311087455085	0.0532235042971549	0.0379158998488012	0.0226082954004476	0.00730069095209403	-0.00800691349625958	-0.0233145179446132	-0.0351464612512344	-0.0469784045578557	-0.0588103478644769	-0.0706422911710982	-0.0824742344777194	-0.0943061777843407	-0.106138121090962	-0.0841095253419820	-0.0620809295930021	-0.0400523338440222	-0.0180237380950423	0.00400485765393767	0.0260334534029176	0.0480620491518975	0.0726647422178179	0.0972674352837384	0.121870128349659	0.146472821415579	0.171075514481500	0.195678207547420	0.220280900613341	0.191431941781425	0.162582982949510	0.133734024117594	0.104885065285679	0.0760361064537636	0.0471871476218482	0.0183381887899328	0.00196634358090327	0.00193300533778266	0.00189966709466205	0.00186632885154145	0.00183299060842084	0.00179965236530023	0.00176631412217962	0.00173297587905901	0.00142000249280716	0.00110702910655532	0.000794055720303478	0.000481082334051635	0.000168108947799793	-0.000144864438452050	-0.000457837824703892	-0.000636094306954722	-0.000814350789205552	-0.000992607271456381	-0.00117086375370721	-0.00134912023595804	-0.00152737671820887	-0.00170563320045970	-0.00134055649225611	-0.000975479784052515	-0.000610403075848922	-0.000245326367645330	0.000119750340558262	0.000484827048761855	0.000849903756965447	0.00136393501041310	0.00187796626386075	0.00239199751730840	0.00290602877075605	0.00342006002420370	0.00393409127765135	0.00444812253109900	0.00387018258931437	0.00329224264752975	0.00271430270574512	0.00213636276396049	0.00155842282217586	0.000980482880391232	0.000402542938606604	0.104422810051322	0.0982817147684688	0.0921406194856154	0.0859995242027620	0.0798584289199086	0.0737173336370552	0.0675762383542019	0.0614351430713484	0.0208876962924624	-0.0196597504864236	-0.0602071972653097	-0.100754644044196	-0.141302090823082	-0.181849537601968	-0.222396984380854	-0.203231267204334	-0.184065550027815	-0.164899832851295	-0.145734115674775	-0.126568398498256	-0.107402681321736	-0.0882369641452162	-0.0356460153462682	0.0169449334526798	0.0695358822516278	0.122126831050576	0.174717779849524	0.227308728648472	0.279899677447420	0.240364725967685	0.200829774487950	0.161294823008215	0.121759871528480	0.0822249200487454	0.0426899685690106	0.00315501708927568	0.00281564107059173	0.00247626505190779	0.00213688903322385	0.00179751301453990	0.00145813699585596	0.00111876097717201	0.000779384958488071	0.000187591945500255	0.00188483209625541	0.00358207224701056	0.00527931239776571	0.00697655254852086	0.00867379269927601	0.0103710328500312	0.0120682730007863	0.0245317074729266	0.0369951419450669	0.0494585764172072	0.0619220108893475	0.0743854453614878	0.0868488798336281	0.0993123143057684	0.0837243558230725	0.0681363973403766	0.0525484388576808	0.0369604803749849	0.0213725218922890	0.00578456340959314	-0.00980339507310272	-0.0249020332439850	-0.0400006714148672	-0.0550993095857494	-0.0701979477566317	-0.0852965859275139	-0.100395224098396	-0.115493862269278	-0.0831680243697269	-0.0508421864701755	-0.0185163485706240	0.0138094893289275	0.0461353272284789	0.0784611651280304	0.110787003027582	0.0962432416500739	0.0816994802725659	0.0671557188950578	0.0526119575175498	0.0380681961400418	0.0235244347625337	0.00898067338502568	1.72913973173758e-06	1.49219979870519e-06	1.25525986567281e-06	1.01831993264042e-06	7.81379999608034e-07	5.44440066575648e-07	3.07500133543261e-07	7.05602005108741e-08	3.00505688095882e-07	5.30451175680890e-07	7.60396663265897e-07	9.90342150850905e-07	1.22028763843591e-06	1.45023312602092e-06	1.68017861360593e-06	2.12394139721460e-06	2.56770418082327e-06	3.01146696443194e-06	3.45522974804061e-06	3.89899253164928e-06	4.34275531525796e-06	4.78651809886663e-06	4.51271422272062e-06	4.23891034657462e-06	3.96510647042861e-06	3.69130259428261e-06	3.41749871813660e-06	3.14369484199060e-06	2.86989096584459e-06	2.46552711650646e-06	2.06116326716832e-06	1.65679941783018e-06	1.25243556849205e-06	8.48071719153910e-07	4.43707869815774e-07	3.93440204776375e-08	3.37234461236893e-08	2.81028717697411e-08	2.24822974157928e-08	1.68617230618446e-08	1.12411487078964e-08	5.62057435394821e-09	0	1.22389907168412	1.24854233054002	1.27318558939591	1.29782884825180	1.32247210710769	1.34711536596358	1.37175862481947	1.39640188367536	1.88787897078481	2.37935605789427	2.87083314500373	3.36231023211318	3.85378731922264	4.34526440633210	4.83674149344155	4.66114964813953	4.48555780283751	4.30996595753548	4.13437411223346	3.95878226693144	3.78319042162941	3.60759857632739	4.38186177593748	5.15612497554757	5.93038817515766	6.70465137476775	7.47891457437784	8.25317777398793	9.02744097359802	8.84109654291511	8.65475211223220	8.46840768154929	8.28206325086639	8.09571882018348	7.90937438950057	7.72302995881766	7.96143946670286	8.19984897458806	8.43825848247325	8.67666799035845	8.91507749824365	9.15348700612885	9.39189651401404	1.22662815061331];

pH=7;
CPMS=1; % mM
CSO4=0; % mM
CH=10^-pH*1000; %mM
CK=CPMS+2*CSO4-CH; % mM
%% Membrane
% e=0.06; % effective porosity
% Lm=200*10^-9; %150 membrane thickness (m)
% X=-50; %fix charge (mM)
% Pd=1; % pore size (nm)
Lm=120*10^-9; %150 membrane thickness (m)
e=0.06; % effective porosity

X=-50; %fix charge (mM)
Pd=1; % pore size (nm)

VF_LMH=[50];
% VF_LMH=[10	20	30 40	50	60	80	100	125	150];
VF=VF_LMH/3.6/10^6;

% Plate_num = 21;
Plate_num = 43;
del_x = Lm/(Plate_num-1);
% New_plate_number=85; % if initial  5 plate, need to add X*4 points
% New_Init_guess=zeros(1,5*New_plate_number+1);
% Pts_middle=(New_plate_number-Plate_num)/(Plate_num-1);  %=19
%% Ions
d_K=0.125*2; % K+ size (nm)
d_H=0.0264*2; % K+ size (nm)
d_PMS=0.213*2; % PMS- size (nm)
d_SO4=0.230*2; % PMS- size (nm)

D_K=1.96*10^-9; % K+ diffusivity (m2/s) 1.96   
D_H=9.312*10^-9;
D_PMS=1.15*10^-9; % PMS- diffusivity (m2/s) 1.15   
D_SO4=1.065*10^-9; % SO4- diffusivity (m2/s)

%% Reaction
Total_eq_in_rxn=2000;
Ratio_Surf_in=1/20;  % percentage of surf/percentage of Inside
 % Rxn_const_surf=240*10^-6;
Rxn_const_surf=Total_eq_in_rxn*Ratio_Surf_in/(Ratio_Surf_in+1)/1000*120*10^-6;  %Typical value 0.00001 0.000031299 OK
%For 1000 in membrane it would be equivalent to 1.2*10^-4 on surface. So
%120X10-6. So if surf-10 X10-6,r+ 1000/12
% r=0;
r=Total_eq_in_rxn/(Ratio_Surf_in+1);
% r=Rxn_const_surf/(1-0.24/2)*2/(0.5*Pd*10^-9)*0.24;  % per total membrane volume
Rxn_const=r*del_x; %Rate constant adjsuted for length of each slice


%% Define empty matrix
P_SF=zeros(length(Rxn_const_surf),length(VF));
C_m_PMS=zeros(length(Rxn_const_surf),length(VF));
C_m_SO4=zeros(length(Rxn_const_surf),length(VF));
C_m_H=zeros(length(Rxn_const_surf),length(VF));
C_m_K=zeros(length(Rxn_const_surf),length(VF));
Surf_RXN=zeros(length(Rxn_const_surf),length(VF));
Find_soln=zeros(length(Rxn_const_surf),length(VF));
XI=zeros(length(Rxn_const_surf),length(VF));
Numb_of_Iter=zeros(length(Rxn_const_surf),length(VF));
J=zeros(length(Rxn_const_surf),length(VF));
J_PMS=zeros(length(Rxn_const_surf),length(VF));
J_ratio=zeros(length(Rxn_const_surf),length(VF));
J_PMS_Conv=zeros(length(Rxn_const_surf),length(VF));
J_PMS_Diff=zeros(length(Rxn_const_surf),length(VF));
J_PMS_Mig=zeros(length(Rxn_const_surf),length(VF));
J_PMS_Diff_check=zeros(length(Rxn_const_surf),length(VF));
% J_PMS_sum=zeros(length(Rxn_const_surf),length(VF));
% J_PMS_ratio_linear_analytic=zeros(length(Rxn_const_surf),length(VF));
Donnan_surf=zeros(length(Rxn_const_surf),length(VF));
Donnan_permeate=zeros(length(Rxn_const_surf),length(VF));

C_p_PMS=zeros(length(Rxn_const_surf),length(VF));
C_p_SO4=zeros(length(Rxn_const_surf),length(VF));
C_p_H=zeros(length(Rxn_const_surf),length(VF));
C_p_K=zeros(length(Rxn_const_surf),length(VF));

C_K=zeros(1,Plate_num);
C_H=zeros(1,Plate_num);
C_PMS=zeros(1,Plate_num); % 
C_SO4=zeros(1,Plate_num); % 
SUM_RXN=zeros(1,Plate_num);
phi=zeros(1,Plate_num);

C_K_record=zeros(1,Plate_num);
C_H_record=zeros(1,Plate_num);
C_PMS_record=zeros(1,Plate_num); 
C_SO4_record=zeros(1,Plate_num); 
phi_record=zeros(1,Plate_num);

Boundary=zeros(Plate_num*6+1,2);
Boundary(:,2)=Inf;
Boundary(5*Plate_num+1:end,1)=-Inf;
Boundary=Boundary-1;

%% Loop for rxn
for r_surf=1:length(Rxn_const_surf)
%% Hindrance coefficient
ld_K=d_K/Pd;  %lambda = ri/rp
ld_H=d_H/Pd;  %lambda = ri/rp
ld_PMS=d_PMS/Pd;
ld_SO4=d_SO4/Pd;

% Nonelectrostatic hindrance coefficient
Kne_K=(1+9/8*ld_K*log(ld_K)-1.56034*ld_K+0.528155*ld_K^2+1.91521*ld_K^3-2.81903*ld_K^4+0.270788*ld_K^5+1.10115*ld_K^6-0.435933*ld_K^7)/(1-ld_K)^2;%  (1+9/8 ld ln ld -1.56 ld +0.53 ld^2 + 1.95 ld^3 - 2.82 ld^4 +0.27 ld^5 +1.10 ld^6 -0.44 ld^7)/(1-ld)^2
Kne_H=(1+9/8*ld_H*log(ld_H)-1.56034*ld_H+0.528155*ld_H^2+1.91521*ld_H^3-2.81903*ld_H^4+0.270788*ld_H^5+1.10115*ld_H^6-0.435933*ld_H^7)/(1-ld_H)^2;
Kne_PMS=(1+9/8*ld_PMS*log(ld_PMS)-1.56034*ld_PMS+0.528155*ld_PMS^2+1.91521*ld_PMS^3-2.81903*ld_PMS^4+0.270788*ld_PMS^5+1.10115*ld_PMS^6-0.435933*ld_PMS^7)/(1-ld_PMS)^2;
Kne_SO4=(1+9/8*ld_SO4*log(ld_SO4)-1.56034*ld_SO4+0.528155*ld_SO4^2+1.91521*ld_SO4^3-2.81903*ld_SO4^4+0.270788*ld_SO4^5+1.10115*ld_SO4^6-0.435933*ld_SO4^7)/(1-ld_SO4)^2;

% Electrostatic hindrance coefficient
A_Ke = 0.003*abs(X)^(2/3);  % X use mM in the equation
if X > 0 %positively charged
Ke_K=1; 
Ke_H=1; 
Ke_PMS=exp(-A_Ke);  %exp(-Zi^2 A)
Ke_SO4=exp(-A_Ke*2^2); 
else %negatively charged
Ke_K=exp(-A_Ke); 
Ke_H=exp(-A_Ke); 
Ke_PMS=1;   
Ke_SO4=1;  
end

%%Frictional factor 
Kf_K=Kne_K*Ke_K; %Ref Ruoyu 2024 EST Li/Mg separation
Kf_H=Kne_H*Ke_H;
Kf_PMS=Kne_PMS*Ke_PMS;
Kf_SO4=Kne_SO4*Ke_SO4;

%% Partitioning coefficient
ps_K=(1-ld_K)^2; % Steric 
ps_H=(1-ld_H)^2; % Steric 
ps_PMS=(1-ld_PMS)^2;
ps_SO4=(1-ld_SO4)^2;

%% Hydrodynamics
% Delta_f=9.41559E-05; % 0.1: 9.41559E-05  0.4: 5.95892E-05  0.7 2.85478E-05  1.0 2.18472E-05  2: 1.29904E-05  5:6.53384E-06
Delta_f=D_PMS*3.6*10^4;
k_mass_trans_PMS=D_PMS/Delta_f; %D/sdl; % mass transfer coefficient (m/s)  100 LMH=100/1000/3600=1/3.6*10^-4 m/s kf=D/delta_f
k_mass_trans_K=D_K/Delta_f;
k_mass_trans_H=D_H/Delta_f;
k_mass_trans_SO4=D_SO4/Delta_f;


%% Loop solver
for j=1:length(VF)
vf=VF(j);
C_PMS_gss=zeros(1,50);
C_SO4_gss=zeros(1,50);
C_H_gss=zeros(1,50);

% if r_surf==1
C_PMS_gss(1)=0.8; % initial guess/same as feed
C_SO4_gss(1)=2.3;
C_H_gss(1)=0.05;
% else
% C_PMS_gss(1)=C_m_PMS(r_surf-1,j); % initial guess/same as feed
% C_SO4_gss(1)=C_m_SO4(r_surf-1,j);
% C_H_gss(1)=C_m_H(r_surf-1,j);
% end

for i=1:length(C_PMS_gss)
%% Donnan and BC at feed side
syms d0
eq=(2*C_SO4_gss(i)+C_PMS_gss(i)-C_H_gss(i))*ps_K*exp(-d0)+C_H_gss(i)*ps_H*exp(-d0)-C_PMS_gss(i)*ps_PMS*exp(d0)-2*C_SO4_gss(i)*ps_SO4*exp(2*d0)+X;
[x]=vpasolve(eq,d0);
d0=x; % dimensionless

C_K(1)=(2*C_SO4_gss(i)+C_PMS_gss(i)-C_H_gss(i))*ps_K*exp(-d0); % K at x=0
C_H(1)=C_H_gss(i)*ps_H*exp(-d0);
C_PMS(1)=C_PMS_gss(i)*ps_PMS*exp(d0); % PMS at x=0
C_SO4(1)=C_SO4_gss(i)*ps_SO4*exp(2*d0); % PMS at x=0

%% Solve for ENP
% syms a b1 b2 c h d % a-K b-Cl, c-potential diff, d-donnan diff
a = sym('a', [1, Plate_num]);
b1 = sym('b1', [1, Plate_num]);
b2 = sym('b2', [1, Plate_num]);
c = sym('c', [1, Plate_num]);
h = sym('h', [1, Plate_num]);
syms d % d-donnan diff
sum_rxn=sym('sum_rxn', [1, Plate_num]);

eq(1)=a(1)-C_K(1);
eq(2)=b1(1)-C_PMS(1);
eq(3)=b2(1)-C_SO4(1);
eq(4)=h(1)-C_H(1);
eq(5)=c(1)-d0;
eq(6)=sum_rxn(Plate_num);

for n=Plate_num-1:-1:1
eq(6+n)=-sum_rxn(n)+0.5*(b1(n+1)+b1(n))*Rxn_const+sum_rxn(n+1);
end

for l=1:Plate_num-1
eq(5*(l)+1+Plate_num)=-D_K*Kf_K*e*(a(l+1)-a(l))/del_x-D_K*Kf_K*e*0.5*(a(l+1)+a(l))*(c(l+1)-c(l))/del_x+Kf_K*0.5*(a(l+1)+a(l))*vf-vf*a(Plate_num)/(ps_K*exp(-d));
eq(5*(l)+2+Plate_num)=-D_H*Kf_H*e*(h(l+1)-h(l))/del_x-D_H*Kf_H*e*0.5*(h(l+1)+h(l))*(c(l+1)-c(l))/del_x+Kf_H*0.5*(h(l+1)+h(l))*vf-vf*h(Plate_num)/(ps_H*exp(-d))+sum_rxn(l);
eq(5*(l)+3+Plate_num)=D_PMS*Kf_PMS*e*(b1(l+1)-b1(l))/del_x+(-1)*D_PMS*Kf_PMS*e*0.5*(b1(l+1)+b1(l))*(c(l+1)-c(l))/del_x-Kf_PMS*0.5*(b1(l+1)+b1(l))*vf+vf*b1(Plate_num)/(ps_PMS*exp(d))+sum_rxn(l); %flux negative so consumption positive
eq(5*(l)+4+Plate_num)=D_SO4*Kf_SO4*e*(b2(l+1)-b2(l))/del_x+(-2)*D_SO4*Kf_SO4*e*0.5*(b2(l+1)+b2(l))*(c(l+1)-c(l))/del_x-Kf_SO4*0.5*(b2(l+1)+b2(l))*vf+vf*b2(Plate_num)/(ps_SO4*exp(2*d))-sum_rxn(l); %flux negative so production negative
eq(5*(l)+5+Plate_num)=a(l+1)+h(l+1)-b1(l+1)-2*b2(l+1)+X;
end

eq(6*Plate_num+1)=a(Plate_num)/(ps_K*exp(-d))+h(Plate_num)/(ps_H*exp(-d))-b1(Plate_num)/(ps_PMS*exp(d))-2*b2(Plate_num)/(ps_SO4*exp(2*d));

% if i==1
% sol=vpasolve(eq,a,h,b1,b2,sum_rxn,c,d,Boundary);
% Init_guss=struct2cell(sol);
% else
sol=vpasolve(eq,a,h,b1,b2,sum_rxn,c,d,Init_guss);
Init_guss=struct2cell(sol);
% end
% disp(sol)

for m=1:Plate_num
    C_K(m)= double([sol.(['a', num2str(m)])]);
    C_H(m)= double([sol.(['h', num2str(m)])]);
    C_PMS(m)= double([sol.(['b1', num2str(m)])]);
    C_SO4(m)= double([sol.(['b2', num2str(m)])]);
    phi(m)=double([sol.(['c', num2str(m)])]);
    SUM_RXN(m)=double([sol.(['sum_rxn', num2str(m)])]);
end

phi_diff=phi(Plate_num)-phi(1); % potential diff in pore dimensionless
d1=double(sol.d);
disp("b1=PMS, b2=SO4")
cp_H=C_H(Plate_num)/(ps_H*exp(-d1));
cp_K=C_K(Plate_num)/(ps_K*exp(-d1));
cp_PMS=C_PMS(Plate_num)/(ps_PMS*exp(d1));
cp_SO4=C_SO4(Plate_num)/(ps_SO4*exp(2*d1));

%% CP solver for membrane surface conc
% Initial guess 4
% if i==1
%     if r_surf==1
%      Cmx0=[3,3,13,0.3,0.1];  
%     else
%      Cmx0=[C_m_PMS(r_surf-1,j),C_m_SO4(r_surf-1,j),C_m_K(r_surf-1,j),C_m_H(r_surf-1,j),XI(r_surf-1,j)];  %Boundry 5 or 8*CPMS
%     end
% else
%     if r_surf==1 
%      Cmx0=[Cm_Sol(1),Cm_Sol(2),Cm_Sol(3),Cm_Sol(4),0.1];  %Boundry 5 or 8*CPMS
%     else
%      Cmx0=[Cm_Sol(1),Cm_Sol(2),Cm_Sol(3),Cm_Sol(4),XI(r_surf-1,j)];  %Boundry 5 or 8*CPMS
%     end
% end

if i==1
    if r_surf==1
     Cmx0=[CPMS,CSO4,CK,CH,0.1];  
    else
     Cmx0=[C_m_PMS(r_surf-1,j),C_m_SO4(r_surf-1,j),C_m_K(r_surf-1,j),C_m_H(r_surf-1,j),XI(r_surf-1,j)]; 
    end
else
    if r_surf==1 
     % Cmx0=[Cm_Sol(1),Cm_Sol(2),Cm_Sol(3),Cm_Sol(4),Cm_Sol(5)];  
     % Cmx0=[Cm_Sol(1),Cm_Sol(2),Cm_Sol(3),Cm_Sol(4),0.1]; %Used for most cases.
     Cmx0=[CPMS,CSO4,CK,CH,0.1];  
    else
     % Cmx0=[Cm_Sol(1),C_m_SO4(r_surf-1,j),Cm_Sol(3),Cm_Sol(4),XI(r_surf-1,j)];  %Boundry 5 or 8*CPMS
     % Cmx0=[Cm_Sol(1),C_m_SO4(r_surf-1,j),Cm_Sol(3),Cm_Sol(4),Cm_Sol(5)];  %Boundry 5 or 8*CPMS
     Cmx0=[Cm_Sol(1),C_m_SO4(r_surf-1,j),Cm_Sol(3),Cm_Sol(4),XI(r_surf-1,j)];  %Boundry 5 or 8*CPMS
    end
end


% Cmx0=[CPMS,CSO4,CK,CH,500];
% Cmx0=[C_PMS_gss(i),C_SO4_gss(i),CK,CH,0]; 
% options = optimset('Display','off');
Eqn_name=@(Cmx)myEquations(Cmx,D_K,D_PMS,D_H,D_SO4,vf,k_mass_trans_K,k_mass_trans_PMS,k_mass_trans_SO4,k_mass_trans_H,Rxn_const_surf,r_surf,SUM_RXN,cp_PMS,cp_SO4,cp_K,cp_H,CK,CPMS,CH,CSO4);
% Cm_Sol=fsolve(Eqn_name, Cmx0,options);
options=optimset('Display','off','TolFun',1e-20); %10^-20
[Cm_Sol,favl,exitflag]=fsolve(Eqn_name, Cmx0,options);

% Cal_PMS,Cal_SO4,Cal_K,Cal_H,XI --> Cmx(1),Cmx(2),Cmx(3),Cmx(4),Cmx(5)
C_PMS_Cal=abs(Cm_Sol(1)); % membran surface concentration (mol/m3)
C_SO4_Cal=abs(Cm_Sol(2)); % membran surface concentration (mol/m3)
C_K_Cal=abs(Cm_Sol(3));
C_H_Cal=abs(Cm_Sol(4));
disp("PMS, SO4, K, H")
disp(Cm_Sol)

if (abs(C_PMS_Cal-C_PMS_gss(i))/C_PMS_Cal<1E-4) && (abs(C_SO4_Cal-C_SO4_gss(i))/C_SO4_Cal<1E-4)  
    break
else
    
    if i<20 
    C_PMS_gss(i+1)=(C_PMS_Cal*C_PMS_gss(i))^0.5;
    else
    C_PMS_gss(i+1)=0.5*(C_PMS_Cal+C_PMS_gss(i));
    end
    C_SO4_gss(i+1)=0.5*(C_SO4_Cal+C_SO4_gss(i));
    C_H_gss(i+1)=0.5*(C_H_Cal+C_H_gss(i));

%OPT 2

% C_PMS_gss(i+1)=(C_PMS_Cal*C_PMS_gss(i))^0.5;
% C_SO4_gss(i+1)=0.5*(C_SO4_Cal+C_SO4_gss(i));
% C_H_gss(i+1)=0.5*(C_H_Cal+C_H_gss(i));

% C_PMS_gss(i+1)=0.5*(C_PMS_Cal+C_PMS_gss(i));
% C_SO4_gss(i+1)=0.5*(C_SO4_Cal+C_SO4_gss(i));
% C_H_gss(i+1)=0.5*(C_H_Cal+C_H_gss(i));

% C_PMS_gss(i+1)=0.2*C_PMS_Cal+C_PMS_gss(i)*0.8;
% C_SO4_gss(i+1)=0.2*C_SO4_Cal+C_SO4_gss(i)*0.8;
% C_H_gss(i+1)=0.2*C_H_Cal+C_H_gss(i)*0.8;

% C_PMS_gss(i+1)=0.1*C_PMS_Cal+C_PMS_gss(i)*0.9;
% C_SO4_gss(i+1)=0.1*C_SO4_Cal+C_SO4_gss(i)*0.9;
% C_H_gss(i+1)=0.1*C_H_Cal+C_H_gss(i)*0.9;

% % 
% C_PMS_gss(i+1)=0.05*C_PMS_Cal+C_PMS_gss(i)*0.95;
% C_SO4_gss(i+1)=0.05*C_SO4_Cal+C_SO4_gss(i)*0.95;
% C_H_gss(i+1)=0.05*C_H_Cal+C_H_gss(i)*0.95;


end 
disp(i)
end


%% Initial guess
%--------------------------------------------------------------------%
% Delta_K=(C_K(2:Plate_num)-C_K(1:Plate_num-1))/(Pts_middle+1);
% Delta_H=(C_H(2:Plate_num)-C_H(1:Plate_num-1))/(Pts_middle+1);
% Delta_PMS=(C_PMS(2:Plate_num)-C_PMS(1:Plate_num-1))/(Pts_middle+1);
% Delta_SO4=(C_SO4(2:Plate_num)-C_SO4(1:Plate_num-1))/(Pts_middle+1);
% Delta_RXN=(SUM_RXN(2:Plate_num)-SUM_RXN(1:Plate_num-1))/(Pts_middle+1);
% Delta_phi=(phi(2:Plate_num)-phi(1:Plate_num-1))/(Pts_middle+1);
% 
% for p=1:Pts_middle+1 %upper limit 20 points includng 1st one
% new_K(p:Pts_middle+1:New_plate_number-1)=C_K(1:end-1)+(p-1)*Delta_K;
% end
% new_K(New_plate_number)=C_K(end);
% 
% for p=1:Pts_middle+1 %upper limit 20 points includng 1st one
% new_H(p:Pts_middle+1:New_plate_number-1)=C_H(1:end-1)+(p-1)*Delta_H;
% end
% new_H(New_plate_number)=C_H(end);
% 
% for p=1:Pts_middle+1 %upper limit 20 points includng 1st one
% new_PMS(p:Pts_middle+1:New_plate_number-1)=C_PMS(1:end-1)+(p-1)*Delta_PMS;
% end
% new_PMS(New_plate_number)=C_PMS(end);
% 
% for p=1:Pts_middle+1 %upper limit 20 points includng 1st one
% new_SO4(p:Pts_middle+1:New_plate_number-1)=C_SO4(1:end-1)+(p-1)*Delta_SO4;
% end
% new_SO4(New_plate_number)=C_SO4(end);
% 
% for p=1:Pts_middle+1 %upper limit 20 points includng 1st one
% new_RXN(p:Pts_middle+1:New_plate_number-1)=SUM_RXN(1:end-1)+(p-1)*Delta_RXN;
% end
% new_RXN(New_plate_number)=SUM_RXN(end);
% 
% for p=1:Pts_middle+1 %upper limit 20 points includng 1st one
% new_phi(p:Pts_middle+1:New_plate_number-1)=phi(1:end-1)+(p-1)*Delta_phi;
% end
% new_phi(New_plate_number)=phi(end);
% 
% New_Init_guess(1:New_plate_number)=new_K;
% New_Init_guess(New_plate_number+1:2*New_plate_number)=new_H;
% New_Init_guess(New_plate_number*2+1:3*New_plate_number)=new_PMS;
% New_Init_guess(New_plate_number*3+1:4*New_plate_number)=new_SO4;
% New_Init_guess(New_plate_number*4+1:5*New_plate_number)=new_RXN;
% New_Init_guess(New_plate_number*5+1:6*New_plate_number)=new_phi;
% New_Init_guess(6*New_plate_number+1)=d1;
%--------------------------------------------------------------------%

%% Pressure
A_lmhb = 5; %L/m2/hr/bar  1 lmH = (1/3.6)*10^-11 m s-1 pa-1
c1=mean(C_K);
c2=mean(C_PMS);
c3=mean(C_SO4);
c4=mean(C_H);
f1=1/D_K; % friction with fluid (s/m2)
f2=1/D_PMS; % friction with fluid (s/m2)
f3=1/D_SO4; % friction with fluid (s/m2)
f4=1/D_H; % friction with fluid (s/m2)
K_f_avg=(Kf_K*c1*(+1)+Kf_PMS*c2*(-1)+Kf_SO4*c3*(-2)+Kf_H*c4*(+1))/(c1*(+1)+c2*(-1)+c3*(-2)+c4*(+1));
f=(A_lmhb/3.6*10^-11*8.314*298*Lm)^-1; % fluid membrane friction (mol*s/m5) Ruoyu EST eq.6 %% f=(A_W RTL/e)^-1
p=-f*vf*Lm-c1*f1*(1-Kf_K)*vf/e*Lm-c2*f2*(1-Kf_PMS)*vf/e*Lm-c3*f3*(1-Kf_SO4)*vf/e*Lm-c4*f4*(1-Kf_H)*vf/e*Lm-Kf_K*(C_K(Plate_num)-C_K(1))-Kf_PMS*(C_PMS(Plate_num)-C_PMS(1))-Kf_SO4*(C_SO4(Plate_num)-C_SO4(1))-Kf_H*(C_H(Plate_num)-C_H(1))+K_f_avg*X*phi_diff; % -f*vf*Lm+X*phi_p mole/m3
ph=-p+(C_PMS_Cal+C_SO4_Cal+C_H_Cal+C_K_Cal)-(cp_PMS+cp_SO4+cp_H+cp_K); %+(Cna1+Ccl1)-(Cna0+Ccl0) mole/m3
Ph=ph*1*(8.3144*298)/1000*0.01; % Pa-->bar
P_SF(r_surf,j)=Ph; %Pressure from solution-friction

Find_soln(r_surf,j)=exitflag;
C_m_PMS(r_surf,j)=C_PMS_Cal;
C_m_SO4(r_surf,j)=C_SO4_Cal;
C_m_H(r_surf,j)=C_H_Cal;
C_m_K(r_surf,j)=C_K_Cal;
Surf_RXN(r_surf,j)=C_m_PMS(r_surf,j)*Rxn_const_surf(r_surf);
XI(r_surf,j)=Cm_Sol(5);
Numb_of_Iter(r_surf,j)=i;

C_p_PMS(r_surf,j)=cp_PMS;
C_p_SO4(r_surf,j)=cp_SO4;
C_p_H(r_surf,j)=cp_H;
C_p_K(r_surf,j)=cp_K;

J_K=vf*cp_K;
J_PMS(r_surf,j)=vf*cp_PMS+SUM_RXN(1)+Rxn_const_surf(r_surf)*C_PMS_Cal;
J_H=vf*cp_H-SUM_RXN(1)-Rxn_const_surf(r_surf)*C_PMS_Cal;
J_SO4=vf*cp_SO4-SUM_RXN(1)-Rxn_const_surf(r_surf)*C_PMS_Cal;
J(r_surf,j)=J_K+J_H-J_PMS(r_surf,j)-2*J_SO4;
J_ratio(r_surf,j)=J(r_surf,j)/J_PMS(r_surf,j);

J_PMS_Conv(r_surf,j)=vf*0.5*(CPMS+C_PMS_Cal);
J_PMS_Mig(r_surf,j)=-1*0.5*(CPMS+C_PMS_Cal)*D_PMS*XI(r_surf,j);
J_PMS_Diff(r_surf,j)=J_PMS(r_surf,j)-J_PMS_Conv(r_surf,j)-J_PMS_Mig(r_surf,j);
J_PMS_Diff_check(r_surf,j)=D_PMS*(CPMS-C_PMS_Cal)/Delta_f;

Donnan_surf(r_surf,j)=d0;
Donnan_permeate(r_surf,j)=d1;
lmh=vf*3.6*10^6; % L/m2/hr

In_rxn(j,:)=SUM_RXN(1:end-1)-SUM_RXN(2:end);
Sum_RXN_In_record(j,:)=SUM_RXN;
C_K_record(j,:)=C_K;
C_H_record=C_H;
C_PMS_record=C_PMS;
C_SO4_record=C_SO4;
phi_record=phi;

for i_influx=1:Plate_num-1
J_PMS_IN_conv(j,i_influx)=Kf_PMS*0.5*(C_PMS(i+1)+C_PMS(i))*vf;
J_PMS_IN_diff(j,i_influx)=(-1)*D_PMS*Kf_PMS*e*(C_PMS(i+1)-C_PMS(i))/del_x;
J_PMS_IN_mig(j,i_influx)=-1*(-1)*D_PMS*Kf_PMS*e*0.5*(C_PMS(i+1)+C_PMS(i))*(phi(i+1)-phi(i))/del_x;
end

end
end



%------------------------------------------------------------------
% Define the system of equations  % Cmx1-5 Cal_PMS,Cal_SO4,Cal_K,Cal_H,XI); % Cmx1-5 Cmx(1),Cmx(2),Cmx(3),Cmx(4),Cmx(5)); 
function F = myEquations(Cmx,D_K,D_PMS,D_H,D_SO4,vf,k_mass_trans_K,k_mass_trans_PMS,k_mass_trans_SO4,k_mass_trans_H,Rxn_const_surf,r_surf,SUM_RXN,cp_PMS,cp_SO4,cp_K,cp_H,CK,CPMS,CH,CSO4)
F(1)=((1+D_K*Cmx(5)/vf)*abs(Cmx(3))-cp_K)/((1+D_K*Cmx(5)/vf)*CK-cp_K)-exp(vf/k_mass_trans_K*(1+D_K*Cmx(5)/vf)); %K
F(2)=((1-D_PMS*Cmx(5)/vf)*abs(Cmx(1))-cp_PMS-SUM_RXN(1)/vf-Rxn_const_surf(r_surf)*abs(Cmx(1))/vf)/((1-D_PMS*Cmx(5)/vf)*CPMS-cp_PMS-SUM_RXN(1)/vf-Rxn_const_surf(r_surf)*abs(Cmx(1))/vf)-exp(vf/k_mass_trans_PMS*(1-D_PMS*Cmx(5)/vf)); %PMS
F(3)=((1+D_H*Cmx(5)/vf)*abs(Cmx(4))-cp_H+SUM_RXN(1)/vf+Rxn_const_surf(r_surf)*abs(Cmx(1))/vf)/((1+D_H*Cmx(5)/vf)*CH-cp_H+SUM_RXN(1)/vf+Rxn_const_surf(r_surf)*abs(Cmx(1))/vf)-exp(vf/k_mass_trans_H*(1+D_H*Cmx(5)/vf)); %H
F(4)=((1-2*D_SO4*Cmx(5)/vf)*abs(Cmx(2))-cp_SO4+SUM_RXN(1)/vf+Rxn_const_surf(r_surf)*abs(Cmx(1))/vf)/((1-2*D_SO4*Cmx(5)/vf)*CSO4-cp_SO4+SUM_RXN(1)/vf+Rxn_const_surf(r_surf)*abs(Cmx(1))/vf)-exp(vf/k_mass_trans_SO4*(1-2*D_SO4*Cmx(5)/vf)); %SO4
F(5)=abs(Cmx(1))+2*abs(Cmx(2))-abs(Cmx(3))-abs(Cmx(4));
% F(6)=(vf*cp_K)+(vf*cp_H-SUM_RXN(1)-Rxn_const_surf(r_surf)*abs(Cmx(1)))-(vf*cp_PMS+SUM_RXN(1)+Rxn_const_surf(r_surf)*abs(Cmx(1)))-2*(vf*cp_SO4-SUM_RXN(1)-Rxn_const_surf(r_surf)*abs(Cmx(1)));
F(6)=1000*((vf*cp_K)+(vf*cp_H-SUM_RXN(1)-Rxn_const_surf(r_surf)*abs(Cmx(1)))-(vf*cp_PMS+SUM_RXN(1)+Rxn_const_surf(r_surf)*abs(Cmx(1)))-2*(vf*cp_SO4-SUM_RXN(1)-Rxn_const_surf(r_surf)*abs(Cmx(1))));
% J_K=vf*cp_K;
% J_PMS=vf*cp_PMS+SUM_RXN(1)+Rxn_const_surf(r_surf)*Cal_PMS;
% J_H=vf*cp_H-SUM_RXN(1)-Rxn_const_surf(r_surf)*Cal_PMS;
% J_SO4=vf*cp_SO4-SUM_RXN(1)-Rxn_const_surf(r_surf)*Cal_PMS;
% J=J_K+J_H-J_PMS-2*J_SO4;
end
%------------------------------------------------------------------------------

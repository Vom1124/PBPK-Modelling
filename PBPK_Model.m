%%%%%%%%%%%%%%%%%%%  PBPK Modelling  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  Vomsheendhur (Vom) Raju %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This PBPK (Physiological Based Pharmacokinetic) model solves the PK
% (Pharmacokiectics) of the drug administered, in this case it will be the
% Contrast Medium (CM). Only iodinated CMs are considered, and therefore,
% the intracellular compartment are ignored as iodinated CMs do not
% penetrate the cells. Exchange between the intravascular and interstitial
% spaces are considered and defined as such using Fick's Law of Diffusion.

%--------------------------------------------------------
clear all; clc; close all; format long;
%--------------------------------------------------------
% simpleApp
% function simpleApp
% % SIMPLEAPP Interactively explore plotting functions
% %   Choose the function used to plot the sample data to see the
% %   differences between surface plots, mesh plots, and waterfall plots
% 
% % Create figure window
% fig = uifigure;
% fig.Name = "My App";
% 
% % Manage app layout
% gl = uigridlayout(fig,[2 3]);
% gl.RowHeight = {30,'1x'};
% gl.ColumnWidth = {'fit','1x'};
% 
% % Create UI components
% lbl = uilabel(gl);
% dd = uidropdown(gl);
% ax = uiaxes(gl);
% ta=uieditfield(gl);
% 
% % Lay out UI components
% % Position label
% lbl.Layout.Row = 1;
% lbl.Layout.Column = 1;
% % Position drop-down
% dd.Layout.Row = 1;
% dd.Layout.Column = 2;
% % Position axes
% ax.Layout.Row = 2;
% ax.Layout.Column = [1 3];
% %Position Selected Text
% ta.Layout.Row=1;
% ta.Layout.Column=3;
% % Configure UI component appearance
% lbl.Text = "Choose Plot Type:";
% dd.Items = ["" "Surf" "Mesh" "Waterfall"];
% % Assign callback function to drop-down
% dd.ValueChangedFcn = {@changePlotType,ax,ta};
% end
% 
% % Program app behavior
% function changePlotType(src,event,ax,ta)
% type = event.Value;
% switch type
%     case "Surf"
%         surf(ax,peaks);
%         
%     case "Mesh"
%         mesh(ax,peaks);
%     case "Waterfall"
%         waterfall(ax,peaks);
% end
% ta.Value=sprintf('%s',type);
% end

%------------------Concentration profiles 
%-- For bolus injection

%% Mass Balance for one compartment modelling
%--------Blood vessel, well-mixed condition Volume is constant and the flow
%in is equal to flow out. 
% Ci = input conc.
% Co = output conc.
% Cc = compartment conc.
% For a well-mixed case, Co=Cc.
% Fick's law: V*(dCo/dt)=Q*(Ci-Co), 
% where 'V' is the volume of the compartment,
% 'Q' is the blood (plasma) flow rate, and
% 't' is the time.

%% Interactive UI

Q=input('Enter the problem to solve=   ');
if Q==1
syms Co1(t) Ci Q V x
ode=diff(Co1,t)==(Q/V)*(Ci-Co1); % Ordinary Differential Equation
Co1=dsolve(ode,Co1(0)==0); %ODE solution
Co1=subs(Co1,[V,Q],[20,10]); % Substituting volume and flow rate
tCi1=[0,0,1,2,3,4,5,6,7,8,9,10,10]; % time range for plot
Ci1=5*ones(1,13); % Input concentration for the above defined time
Ci1(1)=0;Ci1(end)=0;%Initial and final conditions
eqCo1=subs(Co1,{Ci,t},{Ci1(1:end-1),tCi1(1:end-1)});% Finding the output conc.
plot(tCi1,Ci1,'-r','linewidth',2);% plot input conc. vs. time
hold on
%-- Redistributed concentration in the system
syms A Ke tCo2
% "A" is the known concentration at "0" time which is F*D/Vd, where F is the
% bioavailability, D is the Dose, and Vd is the volume of distribution.
% "Ke" is the elimination constant.
Co2=A*exp(-Ke*t); % Plasma conc. after the bolus is done
Ke=0.1386;
ACo2=subs(solve(eval(Co2(1))==eqCo1(end-1),A),t,10);% max. at t=10 sec.
Co2_sub=subs(eval(Co2),A,ACo2);% Plasma conc. after distribution
eqtCo2=solve(Co2_sub==1e-2,t);
tCo2=[10:eqtCo2];
eqCo2=subs(Co2_sub,t,tCo2(1:end));
Co=[eqCo1,eqCo2];
tCo=[tCi1(1:end-1),tCo2];
plot(tCo,Co,'--b')% plot plasma conc. after distribution vs. time
l=legend('\color{red}Ci: Input Conc.\bf',...
  '\color{blue} Co: Output Conc. (bolus geometry + distribution phase)\bf');
l.FontSize=20;
hold off
grid minor
tend=vpa(tCo2(end)+2);
xlim([0 eval(tend)]);
ylim([0 eval(max(Co)+0.2)]);
%%

elseif Q==2
%------------ Organ ODE including Diffusion within organ (intravascular and
%extracellular)
%---- Intravascular space
%Paramters for each organ: Ca, Civ(t) Viv Qiv, Cec(t) Vec Qec P S
% "C": is the conc. of the contrast (in mg/ml).
% "V": is the volume of the space (in ml).
% "Q": is the flowrate via the space (in ml/s).
% "P": is the permeability defined as D/del(x), where "D" is the
% diffusivity and del(x) is the thickness of the membrane separating the
% two spaces (in cm^2/s).
%"S": is the surface area normal to the direction of mass transfer b/w
%spaces.
%"PS": is the permeability surface area product (in ml/min)
% Subscripts refers to the specific organ / blood vessel.
%"i": Coming in to the organ / blood vessel.
%"o": going out from the oran / blood vessel.
% Convert all values in correct units before solving.

%% Testing
% %----- Testing expotential decay using half-life
% ode1=(108.33*(C(2)-C(1))+eval(g))/180;
% ode2=(108.33*((C(1)-C(2))/71))-C(1)*exp(-l*t)/71;
% ode3=(C(3))*exp(-l*t)/524;

%----------Using sym variables
% syms C(i)  g(t) 
% I_conc=320; % Iodine conc. in the contrast medium chosen
% flowrate_i=5; % Initial phase flow rate injected in ml/s. 
% flowrate_f=2; % Second-phase flow rate injected in ml/s. 
% g(t)=flowrate_i*I_conc*(t>0 & t<=25);
% V_RH=180;
% Q_B=6500/60;
% V_iv=12493;
% V_ecf=7987;
% PS=Q_B;
% ode(1)=OCM(Q_B,C(2),C(1),V_RH)+eval(g/V_RH);
% [ode(2),ode(3)]=TCM(Q_B,C(1),C(2),V_iv,C(3),V_ecf,PS);
% % % % % 
% % % % % %---- Using parametric form
% % % % % % ode(1)=(108.33*(C(2)-C(1))+g)/180;
% % % % % % ode(2)=(108.33*(C(1)-C(2))-110*(C(2)-C(3)))/71;
% % % % % % ode(3)=(110/524)*(C(2)-C(3));
% % % % % 
% tic
% odefun=@(t,C) eval(ode)';
% [t,CC]=ode45(odefun,(0:15:300), [0,0,0]');
% toc
% Cf=((V_iv*CC(:,2))+(V_ecf*CC(:,3)))/(V_iv+V_ecf);
% plot(t,25*Cf)
% hold on 
% plot(t,25*CC(:,1),'-or')
% legend('Rest of the body','Right Heart')
% xlabel('Time after injection')
% ylabel('Enhancement in H.U.')
% title('PS/Q=1')

%% Global Ciculation Model
%Refer to the Circulation schematic image for the numbering scheme.
%---Initializing the ode variable to create and store all ODEs.
% ode=sym.empty(44,0);

syms C(i) g(t) t
I_conc=320; % Iodine conc. in the contrast medium chosen
% flowrate_i=input("Enter the flow rate here ="); % Initial phase flow rate injected in ml/s. 
flowrate_i=5;
flowrate_f=1; % Second-phase flow rate injected in ml/s. 
g(t)=flowrate_i*I_conc*(t>0 & t<=25);

% Bi-phasic injection
% g(t)=I_conc*((flowrate_i*(t>0 & t<=20)) + (flowrate_f*(t>20 & t<=75)));

%% Paramters for each organ and blood vessel
%---- Model parameters defined for 70 kg adult human
CO_std=6500; 
BV_std=4760;

%----- Patient's Information
G_p="M"; % User input gender of the patient ("M"-Male and "F"-Female)
H_p=68; % User input height in inches
W_p=177; % User input weight in lbs.

%------ Model size adjustments
%-- Blood Volume
if G_p=="M"
    %Calculated BV for male with weight range 100-310 lbs and height range
    %60-74 inches.
    BV_cal=33.164*(H_p^0.725)*(W_p^0.425)-1229;
elseif G_p=="F"
    %Calculated BV for female with weight range 80-290 lbs and height range
    %60-74 inches.
    BV_cal=34.85*(H_p^0.725)*(W_p^0.425)-1954;
end
BV_r=(BV_std/BV_cal)^(-1); % Ratio to find the regional blood volumes for patient info.
% BV_r=1;

%-- Cardian Output (CO)
CO_cal=36.36*(H_p^0.725)*(W_p^0.425);%Calculated CO for the entered info.
CO_r=(CO_std/CO_cal)^(1);% Ratio to find the CO with user input patient info.
% CO_r=1;

%----Converting the flow rates to ml/s
%% ODEs Starting with the right heart (from left side to right side of the model)
%-- Flowrates and volumes
Q_b1=162/60/CO_r;   V_b1=40/BV_r;   % B1 (Antecubital Vein)
Q_b2=6500/60/CO_r;  V_b2=180/BV_r;  % B2 (Right Heart)
Q_b3=6500/60/CO_r;  V_b3=130/BV_r;  % B3 (Pulmonary Artery)

%--B4 (Lungs)
Q_b4=6500/60/CO_r;  V_b4=150/BV_r;  V_b5=144/BV_r;  PS_b4=0.015*Q_b4;   
Q_b6=6500/60/CO_r;  V_b6=160/BV_r;  % B6 (Pulmonary Vein)
Q_b7=6500/60/CO_r;  V_b7=180/BV_r;  % B7 (Left Heart)
Q_b8=6500/60/CO_r;  V_b8=100/BV_r;  % B8 (Ascending Aorta)

%--B9 (Coronary)
Q_b9=260/60/CO_r;   V_b9=10/BV_r;   V_b10=103/BV_r; PS_b9=0.015*Q_b9;

Q_b11=4940/60/CO_r; V_b11=100/BV_r; %B11 (Descending Aorta)

%--B12 (Bronchial)
Q_b12=130/60/CO_r;  V_b12=5/BV_r; V_b13=144/BV_r;   PS_b12=0.05*Q_b12;  

Q_b14=4810/60/CO_r; V_b14=100/BV_r; %B14 (Abdominal Aorta)
Q_b15=455/60/CO_r;  V_b15=20/BV_r;  %B15

%--B16 (Liver)
Q_b16=1885/60/CO_r; V_b16=71/BV_r;  V_b17=524/BV_r;     PS_b16=0.125*Q_b16;

Q_b18=495/60/CO_r;  V_b18=20/BV_r;  %B18
%--B19 (St/Sp/Pan)
Q_b19=495/60/CO_r;  V_b19=19/BV_r;  V_b20=119/BV_r;     PS_b19=0.9*Q_b19; 

Q_b21=935/60/BV_r;  V_b21=20/BV_r;  %B21;
%--B22 (Intestine)
Q_b22=935/60/CO_r;  V_b22=35/BV_r;  V_b23=547/BV_r;     PS_b22=0.9*Q_b22;

Q_b24=325/60/CO_r;  V_b24=20/BV_r;  %B24
%--B25 (Upper Extremities)
Q_b25=325/60/CO_r;  V_b25=12/BV_r;  V_b26=2751/BV_r;    PS_b25=0.05*Q_b25;

Q_b27=2925/60/CO_r; V_b27=80/BV_r;  %B27
Q_b28=1430/60/CO_r; V_b28=20/BV_r;  %B28
%--B29 (Kidney)
Q_b29=1430/60/CO_r; V_b29=54/BV_r;  V_b30=89/BV_r;      PS_b29=0.0125*Q_b29;
Q_b31=0.19*Q_b29;   V_b31=7;      %B31 (Bladder)

Q_b32=975/60/CO_r;  V_b32=20/BV_r;  %B32;
%--B33 (Head)
Q_b33=975/60/CO_r;  V_b33=37/BV_r;  V_b34=484/BV_r;     PS_b33=0.05*Q_b33;

Q_b35=1430/60/CO_r; V_b35=100/BV_r; %B35 (Portal Vein)

Q_b36=1495/60/CO_r; V_b36=200/BV_r; %B26
%--Trunk and Lower Extremities
Q_b37=1495/60/CO_r; V_b37=57/BV_r;  V_b38=11002/BV_r;   PS_b37=0.05*Q_b37;

Q_b39=1885/60/CO_r; V_b39=100/BV_r; %B39
Q_b40=163/60/CO_r;  V_b40=40/BV_r;  %B40
Q_b41=1430/60/CO_r; V_b41=100/BV_r; %B41
Q_b42=1495/60/CO_r; V_b42=1000/BV_r;%B42
Q_b43=975/60/CO_r;  V_b43=80/BV_r;  %B43
Q_b44=2925/60/CO_r; V_b44=700/BV_r; %B44
Q_b45=4810/60/CO_r; V_b45=800/BV_r; %B45

%-----------Populating the odes into single ODE function
tic
% odefun=@(t,C) eval(ode)'; %making a column vector using transpose.
odefun=@(t,C) [ 
    
    %B1 (Antecubital Vein)... ode1
    (Q_b1*(C(25)-C(1))/V_b1)+(eval(g)/V_b1);
    
    %B2 (Right Heart)... ode2
    (  (Q_b43*C(43)) + (Q_b1*C(1)) + (Q_b40*C(40)) + (Q_b9*C(9)) + ...
    (Q_b45*C(45)) + (Q_b12*C(12)) - (Q_b2*C(2))  )/V_b2;
    
    %B3... ode3
    Q_b3*(C(2)-C(3))/V_b3;
    
    %B4 (Lungs)... ode4 and ode5
    ( (Q_b4*(C(3)-C(4))) - (PS_b4*(C(4)-C(5))) )/V_b4;
    PS_b4*(C(4)-C(5))/V_b5; 
    
    %B6... ode6
    Q_b6*(C(4)-C(6))/V_b6;
    
    %B7 (Left Heart... ode7)
    Q_b7*(C(6)-C(7))/V_b7;
    
    %B8 (Aorta)... ode8
    Q_b8*(C(7)-C(8))/V_b8;
    
    %B9 (Coronary)... ode9 and ode10
    (  (Q_b9*(C(8)-C(9)))  -  (PS_b9*(C(9)-C(10)))  )/V_b9;
    PS_b9*(C(9)-C(10))/V_b10;
    
    %B11... ode11
    Q_b11*(C(8)-C(11))/V_b11;
    
    %B12 (Bronchial)... ode12 and ode13
    (  (Q_b12*(C(11)-C(12)))  -  (PS_b12*(C(12)-C(13)))  )/V_b12;
    PS_b12*(C(12)-C(13))/V_b13;
    
    %B14... ode14
    Q_b14*(C(11)-C(14))/V_b14;
    
    %B15... ode15
    Q_b15*(C(14)-C(15))/V_b15;
    
    %B16 (Liver)... ode16 and ode17
    (  (Q_b15*C(15)) + (Q_b35*C(35)) - (Q_b16*C(16)) - ...
    (PS_b16*(C(16)-C(17)))  )/V_b16;
    PS_b16*(C(16)-C(17))/V_b17;

    %B18... ode18
    Q_b18*(C(14)-C(18))/V_b18;
    
    %B19 (St/Sp/Pan)... ode19 and ode20
    (  (Q_b19*(C(18)-C(19)))  -  (PS_b19*(C(19)-C(20)))  )/V_b19;
    PS_b19*(C(19)-C(20))/V_b20;
    
    %B21... ode21
    Q_b21*(C(14)-C(21))/V_b21;
    
    %B22 (Intestine)... ode22 and ode23
    (  (Q_b22*(C(21)-C(22)))  -  (PS_b22*(C(22)-C(23)))  )/V_b22;
    PS_b22*(C(22)-C(23))/V_b23;
    
    %B24... ode24
    Q_b24*(C(8)-C(24))/V_b24;
    
    %B25 (Upper Extremities)... ode25 and ode26
    (  (Q_b25*(C(24)-C(25)))  -(PS_b25*(C(25)-C(26)))  )/V_b25;
    PS_b25*(C(25)-C(26))/V_b26;
    
    %B27... ode27
    Q_b27*(C(14)-C(27))/V_b27;
    
    %B28... ode28
    Q_b28*(C(27)-C(28))/V_b28;
    
    %B29 (Kidney)... ode29 and ode30
    (   (Q_b29*(C(28)-C(29)))  -  (PS_b29*(C(29)-C(30)))  - ...
    (Q_b31*(C(29)-C(31)))  )/V_b29;
    PS_b29*(C(29)-C(30))/V_b30;
    
    %B31... ode31 (Filtration Into Bladder)
    0.19*Q_b29*(C(29)-C(31))/7;
    
    %B32... ode32
    Q_b32*(C(8)-C(32))/V_b32;
    
    %B33 (Head)... ode33 and ode34
    (  (Q_b33*(C(32)-C(33))) - (PS_b33*(C(33)-C(34)))  )/V_b33;
    PS_b33*(C(33)-C(34))/V_b34;
    
    %B35 (Portal Vein)... ode35
    (  (Q_b22*C(22)) + (Q_b19*C(19)) - (Q_b35*C(35))  )/V_b35;
    
    %B36... ode36
    Q_b36*(C(27)-C(36))/V_b36;
   
    %B37 (Trunk and Lower Extremities)... ode37 and ode38
    (  (Q_b37*(C(36)-C(37))) - (PS_b37*(C(37)-C(38)))  )/V_b37;
    PS_b37*(C(37)-C(38))/V_b38;  

    %B39... ode39
    Q_b39*(C(16)-C(39))/V_b39;
    
    %B40... ode40
    Q_b40*(C(25)-C(40))/V_b40;
    
    %B41... ode41
    Q_b41*(C(29)-C(31)-C(41))/V_b41;
    
    %B42... ode42
    Q_b42*(C(37)-C(42))/V_b42;
    
    %B43... ode43
    Q_b43*(C(33)-C(43))/V_b43;
    
    %B44... ode44
    ( (Q_b41*C(41)) + (Q_b42*C(42)) - (Q_b44*C(44)) )/V_b44;
    
    %B45... ode45
    ( (Q_b39*C(39)) + (Q_b43*C(43)) - (Q_b45*C(45)) )/V_b45;
    
    ];

C0=zeros(1,45);% Initial conditions
C0(1)=eval(g(1))/V_b1;
tspan=(0:7.5:300);%Total time after the start of the injection.
[tt,cc] = ode45(odefun, tspan, C0);%Using Runge-Kutta(4,5)
toc

%% Enhancement 
% Each mg of contrast roughly translates to 25 H.U. of enhancement
%------------Liver
% C_liver=((V_liver_iv*cc(:,31))+(V_liver_ecf*cc(:,32)))...
%     /(V_liver_iv+V_liver_ecf); %Total conc. in the liver
C_liver=(  (V_b16*cc(:,16)) + (V_b17*cc(:,17))  )/(V_b16+V_b17);
E_liver=25*C_liver; %Total enhancement obtainable from the liver

%------------Right Heart 
C_rh=cc(:,2);%Total conc. in the right heart
E_rh=25*C_rh;%Total enhancement obtainable from the right heart

%------------Aorta
C_a=cc(:,8);%Total conc. in the aorta
E_a=25*C_a; %Total enhancement obtainable from the aorta

%------------Coronary
C_co=( (V_b35*cc(:,35)) + (V_b36*cc(:,36))  )/(V_b35+V_b36);%Total conc. in the coronary
E_co=25*C_co; %Total enhancement obtainable from the coronary

%------------Kidney
C_kid=(  (V_b28*cc(:,28)) + (V_b29*cc(:,29))  )/(V_b28+V_b29);

plot(tt,E_a,'-or')
hold on 
plot(tt,E_liver,'-sb')
% plot(t,E_co,'-ok')
legend('Aorta','Liver')
grid on
xlim([0 300]);ylim([0 350]);xticks(0:30:300);yticks(0:50:700);
end 


% %-------------- Right Heart
% Q_rh_out=6500/60/CO_r;Q_rh_c18=975/60/CO_r;Q_rh_c14=162/60/CO_r;
% Q_rh_c13=163/60/CO_r;Q_rh_c8=260/60/CO_r;Q_rh_c20=130/60/CO_r;
% Q_rh_c45=4810/60/CO_r;
% V_rh=180/BV_r;
% % ode(1)=OCM([Q_rh_c18,Q_rh_c14,Q_rh_c13,Q_rh_c8,Q_rh_c20,Q_rh_c45,Q_rh_out],...
% %     [C(18),C(14),C(13),C(8),C(20),C(45)],C(1),V_rh);
% % odef1=(  (Q_rh_c18*C(18)) + (Q_rh_c14*C(14)) + (Q_rh_c13*C(13)) + ...
% %     (Q_rh_c8*C(8)) +  (Q_rh_c20*C(20)) + (Q_rh_c45*C(45)) - ...
% %     (Q_rh_out*C(1))  ) /V_rh;
% 
% %-------------- B12
% Q_B12=6500/60/CO_r;V_B12=130/BV_r;
% % ode(2)=OCM(Q_B12,C(1),C(2),V_B12);
% % odef2=  Q_B12*(C(1)-C(2))/V_B12;
% 
% %-------------- Lungs
% Q_L=6500/60/CO_r;
% V_L_iv=150/BV_r;
% V_L_ecf=144/BV_r;
% PS_L=1.5*Q_L;
% % [ode(3),ode(4)]=TCM(Q_L,C(2),C(3),V_L_iv,C(4),V_L_ecf,PS_L);
% % odef3=( (Q_L*(C(2)-C(3)))-(PS_L*(C(3)-C(4))) )/V_L_iv;
% % odef4= PS_L*(C(3)-C(4))/V_L_ecf;
% 
% %-------------- B35
% Q_B35=6500/60/CO_r;
% V_B35=160/BV_r;
% % ode(5)=OCM(Q_B35,C(3),C(5),V_B35);
% % odef5=Q_B35*(C(3)-C(5))/V_B35;
% 
% %-------------- Left Heart
% Q_lh=6500/60/CO_r;
% V_lh=180/BV_r;
% % ode(6)=OCM(Q_lh,C(5),C(6),V_lh);
% % odef6=Q_lh*(C(5)-C(6))/V_lh;
% 
% %-------------- B67 ("Aorta")
% Q_B67=6500/60/CO_r;
% V_B67=100/BV_r;
% % ode(7)=OCM(Q_B67,C(6),C(7),V_B67);
% % odef7= Q_B67*(C(6)-C(7))/V_B67;
% 
% %----------------------Upper Circulation System
% %-------------- Coronary
% Q_C=260/60/CO_r;
% V_C_iv=10/BV_r;V_C_ecf=103/BV_r;
% PS_C=1.5*Q_C;
% % [ode(8),ode(9)]=TCM(Q_C,C(7),C(8),V_C_iv,C(9),V_C_ecf,PS_C);
% % odef8= (  (Q_C*(C(7)-C(8))) - (PS_C*(C(8)-C(9)))  )/V_C_iv;
% % odef9=PS_C*(C(8)-C(9))/V_C_ecf;
% 
% % -------------- B710
% Q_B710=325/60/CO_r;
% V_B710=20/BV_r;
% % ode(10)= OCM(Q_B710,C(7),C(10),V_B710);
% % odef10= Q_B710*(C(7)-C(10))/V_B710;
% 
% %-------------- Up Ext
% Q_ue=325/60/CO_r;
% V_ue_iv=12/BV_r;V_ue_ecf=2751/BV_r;
% PS_ue=1.5*Q_ue;
% % [ode(11),ode(12)]=TCM(Q_ue,C(10),C(11),V_ue_iv,C(12),V_ue_ecf,PS_ue);
% % odef11=(  (Q_ue*(C(10)-C(11))) - (PS_ue*(C(11)-C(12)))  )/V_ue_iv;
% % odef12= PS_ue*(C(11)-C(12))/V_ue_ecf;
% 
% %-------------- B1113
% Q_B1113=163/60/CO_r;
% V_B1113=40/BV_r;
% % ode(13)=OCM(Q_B1113,C(11),C(13),V_B1113);
% % odef13=Q_B1113*(C(11)-C(13))/V_B1113;
% 
% %-------------- B1114 (Antecubital) Source of Injection
% Q_B1114=162/60/CO_r;
% V_B1114=40/BV_r;
% % ode(14)=OCM(Q_B1114,C(11),C(14),V_B1114)+(eval(g)/V_B1114);
% % odef14=(Q_B1114*(C(11)-C(14))/V_B1114) + (eval(g)/V_B1114);
% 
% %-------------- B715
% Q_B715=975/60/CO_r;
% V_B715=20/BV_r;
% % ode(15)=OCM(Q_B715,C(7),C(15),V_B715);
% % odef15= Q_B715*(C(7)-C(15))/V_B715;
% 
% %-------------- Head
% Q_H=975/60/CO_r;
% V_H_iv=37/BV_r;V_H_ecf=484/BV_r;
% PS_H=1.55*Q_H;
% % [ode(16),ode(17)]=TCM(Q_H,C(15),C(16),...
% %     V_H_iv,C(17),V_H_ecf,PS_H);
% % odef16= (  (Q_H*(C(15)-C(16))) - (PS_H*(C(16)-C(17)))  )/V_H_iv;
% % odef17= PS_H*(C(16)-C(17))/V_H_ecf;
% 
% %-------------- B1618
% Q_B1618=975/60/CO_r;
% V_B1618=80/BV_r;
% % ode(18)=OCM(Q_B1618,C(16),C(18),V_B1618);
% % odef18= Q_B1618*(C(16)-C(18))/V_B1618;
% 
% %------------------------Lower Circulation System
% %-------------- B719
% Q_B719=4940/60/CO_r;
% V_B719=100/BV_r;
% % ode(19)=OCM(Q_B719,C(7),C(19),V_B719);
% % odef19= Q_B719*(C(7)-C(19))/V_B719;
% 
% %-------------- Bronchial
% Q_B=130/60/CO_r;
% V_B_iv=5/BV_r;V_B_ecf=144/BV_r;
% PS_B=1.5*Q_B;
% % [ode(20),ode(21)]=TCM(Q_B,C(7),C(20),V_B_iv,C(21),V_B_ecf,PS_B);
% % odef20= (  (Q_B*(C(19)-C(20))) - (PS_B*(C(20)-C(21)))  )/V_B_iv;
% % odef21= PS_B*(C(20)-C(21))/V_B_ecf;
% 
% 
% %-------------- B1922
% Q_B1922=4810/60/CO_r;
% V_B1922=100/BV_r;
% % ode(22)=OCM(Q_B1922,C(19),C(22),V_B1922);
% % odef22= Q_B1922*(C(19)-C(22))/V_B1922;
% 
% %-------------- B2223
% Q_B2223=455/60/CO_r;
% V_B2223=20/BV_r;
% % ode(23)=OCM(Q_B2223,C(22),C(23),V_B2223);
% % odef23= Q_B2223*(C(22)-C(23))/V_B2223;
% 
% %-------------- B2224
% Q_B2224=495/60/CO_r;
% V_B2224=20/BV_r;
% % ode(24)=OCM(Q_B2224,C(22),C(24),V_B2224);
% % odef24= Q_B2224*(C(22)-C(24))/V_B2224;
% 
% %-------------- St/Sp/Pan
% Q_st=495/60/CO_r;
% V_st_iv=19/BV_r;V_st_ecf=112/BV_r;
% PS_st=1.125*Q_st;
% % [ode(25),ode(26)]=TCM(Q_st,C(24),C(25),V_st_iv,C(26),V_st_ecf,PS_st);
% % odef25= (  (Q_st*(C(24)-C(25))) - (PS_st*(C(25)-C(26)))  )/V_st_iv;
% % odef26= PS_st*(C(25)-C(26))/V_st_ecf;
% 
% %-------------- B2227
% Q_B2227=935/60/CO_r;
% V_B2227=20/BV_r;
% % ode(27)=OCM(Q_B2227,C(22),C(27),V_B2227);
% % odef27= Q_B2227*(C(22)-C(27))/V_B2227;
% 
% %-------------- Intestine
% Q_int=935/60/CO_r;
% V_int_iv=35/BV_r;V_int_ecf=547/BV_r;
% PS_int=1.125*Q_int;
% % [ode(28),ode(29)]=TCM(Q_int,C(27),C(28),V_int_iv,C(29),V_int_ecf,PS_int);
% % odef28=(  (Q_int*(C(27)-C(28))) - (PS_int*(C(28)-C(29)))  )/V_int_iv;
% % odef29= PS_int*(C(28)-C(29))/V_int_ecf;
% 
% %-------------- Portal Vein
% Q_BP_out=1430/60/CO_r; Q_BP_int=935/60/CO_r;Q_BP_st=495/60/CO_r;
% V_BP=100/BV_r;
% % ode(30)=OCM([Q_BP_st,Q_BP_int,Q_BP_out],[C(25),C(28)],C(30),V_BP);
% % odef30= (  (Q_BP_st*C(25)) + (Q_BP_int*C(28)) - (Q_BP_out*C(30))  )/V_BP;
% 
% %-------------- Liver
% Q_liver_out=1885/60/CO_r;Q_liver_B2223=455/60/CO_r;Q_liver_BP=1430/60/CO_r;
% V_liver_iv=71/BV_r;V_liver_ecf=524/BV_r;
% PS_liver=0.65*Q_liver_out;
% % [ode(31),ode(32)]=TCM([Q_liver_B2223,Q_liver_BP,Q_liver_out],...
% %     [C(23),C(30)],C(31),V_liver_iv,C(32),V_liver_ecf,PS_liver);
% % odef31= (  (Q_liver_B2223*C(23)) + (Q_liver_BP*C(30)) - ...
% %     (Q_liver_out*C(31)) - (PS_liver*(C(31)-C(32)))  )/V_liver_iv;
% % odef32= PS_liver*(C(31)-C(32))/V_liver_ecf;
% 
% %-------------- B3133
% Q_B3133=1885/60/CO_r;
% V_B3133=100/BV_r;
% % ode(33)=OCM(Q_B3133,C(31),C(33),V_B3133);
% % odef33= Q_B3133*(C(31)-C(33))/V_B3133;
% 
% %-------------- B2234
% Q_B2234=2925/60/CO_r;
% V_B2234=80/BV_r;
% % ode(34)=OCM(Q_B2234,C(22),C(34),V_B2234);
% % odef34= Q_B2234*(C(22)-C(34))/V_B2234;
% 
% %-------------- B3435
% Q_B3435=1430/60/CO_r;
% V_B3435=20/BV_r;
% % ode(35)=OCM(Q_B3435,C(34),C(35),V_B3435);
% % odef35= Q_B3435*(C(34)-C(35))/V_B3435;
% 
% %-------------- Kidney
% Q_kid=1430/60/CO_r;
% V_kid_iv=54/BV_r;V_kid_ecf=89/BV_r;
% PS_kid=0.19*Q_kid;
% % [ode(36),ode(37)]=TCM(Q_kid,C(35),C(36),...
% %     V_kid_iv,C(37),V_kid_ecf,PS_kid);
% % % ode(36)=(  (Q_kid*(C(35)-C(36))) - (PS_kid*(C(36)-C(37)))  )/V_kid_iv ;
% % ode(37)=PS_kid*(C(36)-C(37))/V_kid_ecf;
% % odef36= (  (Q_kid*(C(35)-C(36))) - (PS_kid*(C(36)-C(37)))  )/V_kid_iv;
% % odef37= PS_kid*(C(36)-C(37))/V_kid_ecf;
% 
% % % Plasma half life for bladder excretion
% % lambda=-log(0.35)/(3600*1); % 35% excreted in 60 min.(3600 sec) for Iopamidol
% % C_kid_decay=C(35)*exp(-lambda/t);% Half-life decay
% % C_kid_decay=0.1*Q_kid*(C(38));
% 
% % --Bladder excretion
% % ode(38)=0.19*Q_kid*C(36)-Q_kid*C(38);
% % odef38= (0.19*Q_kid*C(36)) - (Q_kid*C(38));
% 
% %-------------- B3639
% Q_B3639=1430/60/CO_r;
% V_B3639=100/BV_r;
% % ode(39)=Q_B3639*(C(38)-C(39))/V_B3639;
% % odef39=Q_B3639*(C(38)-C(39))/V_B3639;
% 
% %-------------- B3440
% Q_B3440=1495/60/CO_r;
% V_B3440=200/BV_r;
% % ode(40)=OCM(Q_B3440,C(34),C(40),V_B3440);
% % odef40=Q_B3440*(C(34)-C(40))/V_B3440;
% 
% %-------------- Trunk and Lower Extremities
% Q_le=1495/60/CO_r;
% V_le_iv=57/BV_r;V_le_ecf=11002/BV_r;
% PS_le=1.5*Q_le;
% % [ode(41),ode(42)]=TCM(Q_le,C(40),C(41),V_le_iv,C(42),V_le_ecf,PS_le);
% % odef41=(  (Q_le*(C(40)-C(41))) - (PS_le*(C(41)-C(42)))  )/V_le_iv;
% % odef42= PS_le*(C(41)-C(42))/V_le_ecf;
% 
% %-------------- B4143
% Q_B4143=1495/60/CO_r;
% V_B4143=1000/BV_r;
% % ode(43)=OCM(Q_B4143,C(41),C(43),V_B4143);
% % odef43=Q_B4143*(C(41)-C(43))/V_B4143;
% 
% %-------------- B4344
% Q_B4344_out=2925/60/CO_r;Q_B4344_B39=1430/60/CO_r;Q_B4344_B43=1495/60/CO_r;
% V_B4344=700/BV_r;
% % ode(44)=OCM([Q_B4344_B39,Q_B4344_B43,Q_B4344_out],[C(39),C(43)],C(44),V_B4344);
% % odef44= (  (Q_B4344_B43*C(43)) + (Q_B4344_B39*C(39)) - ...
% %     (Q_B4344_out*C(44))  )/V_B4344;
% 
% %-------------- B4445
% Q_B4445_out=4810/60/CO_r;Q_B4445_B44=2925/60/CO_r;Q_B4445_B33=1885/60/CO_r;
% V_B4445=800/BV_r;
% % ode(45)=OCM([Q_B4445_B33,Q_B4445_B44,Q_B4445_out],[C(33),C(44)],C(45),V_B4445);
% % odef45= (  (Q_B4445_B44*C(44)) + (Q_B4445_B33*C(33)) - ...
% %     (Q_B4445_out*C(45))  )/V_B4445;

% %------ODEs
%   [
        %B2 (Right Heart)
%     (  (Q_rh_c18*C(18)) + (Q_rh_c14*C(14)) + (Q_rh_c13*C(13)) + ...
%     (Q_rh_c8*C(8)) +  (Q_rh_c20*C(20)) + (Q_rh_c45*C(46)) - ...
%     (Q_rh_out*C(1))  ) /V_rh;

      %Coronary
%     (  (Q_C*(C(7)-C(8))) - (PS_C*(C(8)-C(9)))  )/V_C_iv;
%     PS_C*(C(8)-C(9))/V_C_ecf;
%     
%     %B710
%     Q_B710*(C(7)-C(10))/V_B710;
%     
%     %Up Ext
%     (  (Q_ue*(C(10)-C(11))) - (PS_ue*(C(11)-C(12)))  )/V_ue_iv;    
%      PS_ue*(C(11)-C(12))/V_ue_ecf;
%      
%     %B1113
%     Q_B1113*(C(11)-C(13))/V_B1113;
%     
% 
%     
%     %B715
%     Q_B715*(C(7)-C(15))/V_B715;
%     
%     %Head
%     (  (Q_H*(C(15)-C(16))) - (PS_H*(C(16)-C(17)))  )/V_H_iv;
%     PS_H*(C(16)-C(17))/V_H_ecf;
%     
%     %B1618
%     Q_B1618*(C(16)-C(18))/V_B1618;
%     
%     %B719
%     Q_B719*(C(7)-C(19))/V_B719;
%     
%     %Bronchial
%     (  (Q_B*(C(19)-C(20))) - (PS_B*(C(20)-C(21)))  )/V_B_iv;
%     PS_B*(C(20)-C(21))/V_B_ecf;
%     
%     %B1922
%     Q_B1922*(C(19)-C(22))/V_B1922;
%     
%     %B2223
%     Q_B2223*(C(22)-C(23))/V_B2223;
%     
%     %B2224
%     Q_B2224*(C(22)-C(24))/V_B2224;
%     
%     %St/SP/Pan
%     (  (Q_st*(C(24)-C(25))) - (PS_st*(C(25)-C(26)))  )/V_st_iv;
%     PS_st*(C(25)-C(26))/V_st_ecf;
%     
%     %B2227
%     Q_B2227*(C(22)-C(27))/V_B2227;
%     
%     %Intestine
%     (  (Q_int*(C(27)-C(28))) - (PS_int*(C(28)-C(29)))  )/V_int_iv;
%     PS_int*(C(28)-C(29))/V_int_ecf;
%     
%     % Portal Vein
%     (  (Q_BP_st*C(25)) + (Q_BP_int*C(28)) - (Q_BP_out*C(30))  )/V_BP;
%     
%     %Liver
%     (  (Q_liver_B2223*C(23)) + (Q_liver_BP*C(30)) - ...
%     (Q_liver_out*C(31)) - (PS_liver*(C(31)-C(32)))  )/V_liver_iv;
%     PS_liver*(C(31)-C(32))/V_liver_ecf;
%     
%     %B3133
%     Q_B3133*(C(31)-C(33))/V_B3133;
%     
%     %B2234
%     Q_B2234*(C(22)-C(34))/V_B2234;
%     
%     %B3435
%     Q_B3435*(C(34)-C(35))/V_B3435;
%     
%     %Kidney
%     (  (Q_kid*(C(35)-C(36))) - (PS_kid*(C(36)-C(37)))  )/V_kid_iv;
%     PS_kid*(C(36)-C(37))/V_kid_ecf;
%     
%     %Bladder Excretion
%     (  (Q_kid*(C(36)-C(38))) - (0.19*Q_kid*(C(38)-C(39)))  )/300;
%     0.19*Q_kid*(C(38)-C(39))/300;
%     
%     %B3639
%     Q_B3639*(C(38)-C(40))/V_B3639;
%     
%     %B3440
%     Q_B3440*(C(34)-C(41))/V_B3440;
%     
%     %Trunk and Lower Extremities
%     (  (Q_le*(C(41)-C(42))) - (PS_le*(C(42)-C(43)))  )/V_le_iv;
%     PS_le*(C(42)-C(43))/V_le_ecf;
%     
%     %B4143
%     Q_B4143*(C(42)-C(44))/V_B4143;
%     
%     %B4344
%     (  (Q_B4344_B43*C(44)) + (Q_B4344_B39*C(40)) - ...
%     (Q_B4344_out*C(45))  )/V_B4344;
%     
%     %B4445
%     (  (Q_B4445_B44*C(45)) + (Q_B4445_B33*C(33)) - ...
%     (Q_B4445_out*C(46))  )/V_B4445;
%     
%     ];

%% Defining the ODEs for the compartment models
%-----One Compartment Model
function odefun=OCM(Qin,Cin,Civ,Viv)
% OCM is the One Compartment Transport Model.
% Qin is the inflow blood flowate.
% Cin is the inflow CM conc.
% Civ is the intravascular CM conc.
% Viv is the volume of the iv space.
if nargin<4
    error(['Not enough input arguments!' newline...
        'Check if all the input arguments are entered.'])
else
    if numel(Cin)==1
        odefun=Qin*(Cin-Civ)/Viv;
    elseif numel(Cin)==2
        odefun=( (Qin(1)*Cin(1))+(Qin(2)*Cin(2))-(Qin(3)*Civ) )./Viv;
    elseif numel(Cin)==6
        odefun=( (Qin(1)*Cin(1))+(Qin(2)*Cin(2))+(Qin(3)*Cin(3))+...
        (Qin(4)*Cin(4))+(Qin(5)*Cin(5))+(Qin(6)*Cin(6))-(Qin(7)*Civ) )./Viv;
    end
end
end
%-----Two Compartment Model
function [odefun_iv,odefun_ecf]=TCM(Qin,Cin,Civ,Viv,Cecf,Vecf,PS_r)
if nargin<7
    error(['Not enough input arguments!' newline...
        'Check if all the input arguments are entered.'])
else
% TCM is the Two Compartment Transport Model.
% Qin is the inflow blood flowate.
% Cin is the inflow CM conc.
% Civ is the intravascular CM conc.
% Cecf is the extracellular CM conc.
% Viv is the volume of the iv space.
% Vecf is the colume of the ecf space.
% PS_r is the ratio of product of permeability and surface area to the
% blood fow rate.
if numel(Cin)==1
odefun_iv= (  (Qin.*(Cin-Civ)) - (PS_r.*(Civ-Cecf))  )./Viv;    
odefun_ecf=PS_r.*(Civ-Cecf)./Vecf;

elseif numel(Cin)==2
odefun_iv=(  (Qin(1)*Cin(1))+(Qin(2)*Cin(2))-(Qin(3)*Civ)...
    -(PS_r.*(Civ-Cecf))  )./Viv;    
odefun_ecf=PS_r.*(Civ-Cecf)/Vecf;
end
end
end
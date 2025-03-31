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
g(t)=I_conc*(t>0 & t<=25);

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
BV_r=(BV_std/BV_cal)^(1); % Ratio to find the regional blood volumes for patient info.
% BV_r=1;

%-- Cardian Output (CO)
CO_cal=36.36*(H_p^0.725)*(W_p^0.425);%Calculated CO for the entered info.
CO_r=(CO_std/CO_cal)^(1);% Ratio to find the CO with user input patient info.
% CO_r=1;

%----Converting the flow rates to ml/s
%% ODEs Starting with the right heart (from left side to right side of the model)
%-- Flowrates and volumes
Q1=260/60/CO_r;     V1=10/BV_r;    V14=103/BV_r;    ps14=0.005*Q1;
Q2=130/60/CO_r;     V2=5/BV_r;     V16=144/BV_r;    ps16=0.005*Q2;
Q3=162/60/CO_r;     V3=40/BV_r;
Q4=163/60/CO_r;     V4=40/BV_r;
Q5=975/60/CO_r;     V5=80/BV_r;
Q6=4810/60/CO_r;    V6=800/BV_r;
Q7=6500/60/CO_r;    V7=180/BV_r;
Q8=Q7;              V8=130/BV_r;
Q9=6500/60/CO_r;    V9=150/BV_r;    V10=144/BV_r;   ps10=0.005*Q9;
Q11=Q7;             V11=160/BV_r;   
Q12=Q7;             V12=180/BV_r;
Q13=Q7;             V13=100/BV_r;
Q15=4940/60/CO_r;   V15=100/BV_r;    
Q17=4810/60/CO_r;   V17=100/BV_r;
Q18=455/60/CO_r;    V18=20/BV_r;
Q19=495/60/CO_r;    V19=20/BV_r;
Q20=935/60/CO_r;    V20=20/BV_r;
Q21=Q19;            V21=19/BV_r;    V22=112/BV_r;   ps22=0.00250*Q21;
Q23=Q20;            V23=35/BV_r;    V24=547/BV_r;   ps24=0.0025*Q23;
Q25=1430/60/CO_r;   V25=100/BV_r;
Q26=1885/60/CO_r;   V26=71/BV_r;    V27=524/BV_r;   ps27=0.05*Q26;
Q28=1885/60/CO_r;   V28=100/BV_r;
Q29=325/60/CO_r;    V29=20/BV_r;
Q30=Q29;            V30=12/BV_r;    V31=2751/BV_r;  ps31=0.0025*Q30;
Q32=2925/60/CO_r;   V32=80/BV_r;
Q33=1430/60/CO_r;   V33=20/BV_r;
Q34=Q33;            V34=54/BV_r;    V35=89/BV_r;    ps35=0.0019*Q34;
Q36=Q33;            V36=100/BV_r;
Q37=975/60/CO_r;    V37=20/BV_r;
Q38=Q37;            V38=37/BV_r;    V39=484/BV_r;   ps39=0.005*Q38;
Q40=1495/60/CO_r;   V40=200/BV_r;
Q41=Q40;            V41=57/BV_r;    V42=11002/BV_r; ps42=0.0005*Q41;
Q43=Q41;            V43=1000/BV_r;
Q44=Q32;            V44=700/BV_r;

%-----------Populating the odes into single ODE function
tic
% odefun=@(t,C) eval(ode)'; %making a column vector using transpose.
odefun=@(t,C) [ 

    % ode1
    (  (Q1*(C(13)-C(1))) - (ps14*(C(1)-C(14)))  )/V1;
    
    % ode2
    (  (Q2*(C(15)-C(2))) - (ps16*(C(2)-C(16)))  )/V2;
    
    % ode3 
    (  (Q3*C(30)) -  (Q3*C(3))  )/V3  ;
    
    
    % ode4
    Q4*(C(30)-C(4))/V4;
    
    % ode5
    Q5*(C(38)-C(5))/V5;
    
    % ode6
    (  (Q28*C(28)) + (Q44*C(44)) - (Q6*C(6))  )/V6;
    
    % ode7
    ( (Q1*C(1)) + (Q2*(C(2) +(5*C(45)/Q2))) + (Q3*C(3))  + (Q4*C(4)) + (Q5*C(5)) + ...
    (Q6*C(6)) - (Q7*C(7))  )/V7;
    
    % ode8
    Q8*(C(7)-C(8))/V8;
    
    % ode9
    (  (Q9*(C(8)-C(9))) - (ps10*(C(9)-C(10)))  )/V9;
    % ode10
    ps10*(C(9)-C(10))/V10;
    
    % ode11
    Q11*(C(9)-C(11))/V11;
    
    % ode12
    Q12*(C(11)-C(12))/V12;
    
    % ode13
    Q13*(C(12)-C(13))/V13;
    
    % ode14
    ps14*(C(1)-C(14))/V14;
    
    % ode15
    Q15*(C(13)-C(15))/V15;
    
    % ode16
    ps16*(C(2)-C(16))/V16;
    
    % ode17
    Q17*(C(15)-C(17))/V17;
    
    % ode18
    Q18*(C(17)-C(18))/V18;
    
    % ode19
    Q19*(C(17)-C(19))/V19;
    
    % ode20
    Q20*(C(17)-C(20))/V20;
    
    % ode21
    (  (Q21*(C(19)-C(21))) - (ps22*(C(21)-C(22)))  )/V21;
    % ode22
    ps22*(C(21)-C(22))/V22;
    
    % ode23
    (  (Q23*(C(20)-C(23))) - (ps24*(C(23)-C(24)))  )/V23;
    % ode24
    ps24*(C(23)-C(24))/V24;
    
    % ode25 (oBladder)
    0.19*Q33*(C(32)-C(25));
    
    % ode26
    (  (Q25*C(46)) + (Q18*C(18)) - (Q26*C(26)) - (ps27*(C(26)-C(27)))  )/V26;
    % ode27
    ps27*(C(26)-C(27))/V27;
    
    % ode 28
    ((Q28*(C(26)-C(28))/V28)) ;
    
    % ode29
    Q29*(C(13)-C(29))/V29;
    
    % ode30
    (  (Q30*(C(29)-C(30))) - (ps31*(C(30)-C(31)))  )/V30  ;
    % ode31
    ps31*(C(30)-C(31))/V31;
    
    % ode32
    Q32*(C(17)-C(32))/V32;
    
    % ode33
    Q33*(C(32)-C(33))/V33;
        
    % ode34
    (  (Q34*( C(33)- C(34) ))  -  (ps35*(C(34)-C(35)))  )/V34   -  C(25);
    % ode35    /V34
    ps35*(C(34)-C(35))/V35;
    
    % ode36
    (  (Q33*(C(34)))  - (Q36*C(36))  )/V36  ;
    
    % ode37
    Q37*(C(13)-C(37))/V37;
    
    % ode38
    (  (Q38*(C(37)-C(38))) - (ps39*(C(38)-C(39)))  )/V38;
    % ode39
    ps39*(C(38)-C(39))/V39;
    
    % ode40
    Q40*(C(32)-C(40))/V40;
    
    % ode41
    (  (Q41*(C(40)-C(41))) - (ps42*(C(41)-C(42)))  )/V41;
    % ode42
    ps42*(C(41)-C(42))/V42;
    
    % ode43
    Q43*(C(41)-C(43))/V43;
    
    % ode44
    (  (Q36*C(36)) + (Q43*C(43)) - (Q44*C(44))  )/V44;
    
    % ode45
    5*(eval(g)-C(45))/V3;
    
    %ode46 (Portal Vein)
    (  (Q23*C(23)) + (Q21*C(21)) - (Q25*C(46)) )/V25;
    ];

C0=zeros(1,46);% Initial conditions
tspan=(0:7.5:300);%Total time after the start of the injection.
[tt,cc] = ode45(odefun, tspan, C0);%Using Runge-Kutta(4,5)
toc

%% Enhancement 
% Each mg of contrast roughly translates to 25 H.U. of enhancement
%------------Liver
% C_liver=((V_liver_iv*cc(:,31))+(V_liver_ecf*cc(:,32)))...
%     /(V_liver_iv+V_liver_ecf); %Total conc. in the liver
C_liver=(  (V26*cc(:,26)) + (V27*cc(:,27))  )/(V26+V27);
E_liver=25*C_liver; %Total enhancement obtainable from the liver

%------------Right Heart 
C_rh=cc(:,2);%Total conc. in the right heart
E_rh=25*C_rh;%Total enhancement obtainable from the right heart

%------------Aorta
C_a=cc(:,13);%Total conc. in the aorta
E_a=25*C_a; %Total enhancement obtainable from the aorta

%------------Coronary
C_co=( (V1*cc(:,1)) + (V14*cc(:,14))  )/(V1+V14);%Total conc. in the coronary
E_co=25*C_co; %Total enhancement obtainable from the coronary

%------------Kidney
C_kid=(  (V34*cc(:,34)) + (V35*cc(:,35))  )/(V34+V35);

plot(tt,E_a,'-or')
hold on 
plot(tt,E_liver,'-sb')
% plot(t,E_co,'-ok')
legend('Aorta','Liver')
grid on
xlim([0 300]);ylim([0 400]);xticks(0:30:300);yticks(0:50:400);
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
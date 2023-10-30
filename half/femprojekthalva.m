%% setup
clear
% close all
p =  matfile('p.mat').p;
e =  matfile('e.mat').e;
t =  matfile('t.mat').t;





%%  för temperatur.
% a)  bestäm temperaturen

% constanterna
% 
% Tinf  = 40;
% Tinfd = 40;
% ac = 100;  %ac = αc
% Q = 0;    % är Q = 0
% k_g = 0.8;
% k_ti  =  17;
% ko  =  [k_g k_ti];
% L = [];
% D = eye(2);
% Tc = 20;
% rog = 3860;
% roti = 4620;
% crog = 670;
% croti = 523;
% nen = 3;
% coord = p';
% ndof = length(coord);
% enod = t(1:3,:)'; % nodes of elements
% nelm = size(enod,1); % number of elements
% nnod = size(coord,1);
% dof = (1:nnod)'; % dof number is node number
% Konstanter
Tinf =-96;
Tc = 20;

ac = 100;  %ac = αc
Q = 0;    % är Q = 0
k_g = 0.8;
k_ti = 17;
ko = [k_g k_ti];
L = [];
D = eye(2);

rog = 3860;
roti = 4620;
crog = 670;
croti = 523;
nen = 3;
coord = p';
ndof = length(coord);
enod = t(1:3,:)'; % nodes of elements
nelm = size(enod,1); % number of elements
nnod = size(coord,1);
dof = (1:nnod)'; % dof number is node number




%segmenter
conv_segments = [1 8 19 28];
tc_segment =  7;
tiarea = 1; %vilket subdomaint  är Ti


Finit = [];

for ie = 1:nelm
    edof(ie,:) = [ie,enod(ie,:)];
end



[ex,ey] = coordxtr(edof,coord,dof,nen);
K = zeros(ndof);
F =  zeros(ndof,1);
C = zeros(ndof);
thickness = 0.01;



for elnr = 1:nelm      % här assembleras allting
    
    if t(4,elnr) == tiarea
        Ke = flw2te(ex(elnr,:),ey(elnr,:),thickness,ko(2).*D);
        Ce = plantml(ex(elnr,:),ey(elnr,:),croti*roti*thickness);
    else
        Ke= flw2te(ex(elnr,:),ey(elnr,:),thickness,ko(1).*D);
        Ce = plantml(ex(elnr,:),ey(elnr,:),crog*rog*thickness);
    end
    
    indx = edof(elnr,2:end);
    K(indx,indx) = K(indx,indx)+Ke;
    C(indx,indx) =  C(indx,indx)+Ce;
end

convind = [];
% här syns koden för konvektions randvilkor
Finit = F;


for kant = 1:length(e)
    x1 = coord(e(1,kant),1);
    y1 = coord(e(1,kant),2);
    x2 = coord(e(2,kant),1);
    y2 = coord(e(2,kant),2);
    L =  hypot((x1-x2),y1-y2);
    
    if ismember(e(5,kant),conv_segments)
        F(e(1,kant)) = F(e(1,kant)) + thickness*(L/2)*ac*Tinf;
        F(e(2,kant)) = F(e(2,kant)) +thickness* (L/2)*ac*Tinf;
        
        % om det är dag så vill vi beräkna f för natt och viseversa
        if Tinf == -96
            Tinit = 40;
        else
            Tinit = -96;
        end
        Finit(e(1,kant)) =  Finit(e(1,kant)) + thickness*(L/2)*ac*Tinit;
        Finit(e(2,kant)) =  Finit(e(2,kant)) + thickness*(L/2)*ac*Tinit;
        
        
        Kc = thickness*((ac*L)/6)* [2 1;1 2];
        te =  [e(1,kant) e(2,kant)];
        K(te,te) = K(te,te)+Kc;
        
    elseif e(5,kant) ==  tc_segment
        F(e(1,kant)) = F(e(1,kant)) + thickness*(L/2)*ac*Tc;
        F(e(2,kant)) = F(e(2,kant)) + thickness*(L/2)*ac*Tc;
        
        
        Finit(e(1,kant)) = Finit(e(1,kant)) + thickness*(L/2)*ac*Tc;
        Finit(e(2,kant)) = Finit(e(2,kant)) + thickness*(L/2)*ac*Tc;
        
        Kc =thickness*((ac*L)/6)* [2 1;1 2];
        te =  [e(1,kant) e(2,kant)];
        K(te,te) = K(te,te)+Kc;
        
    end
    
end

%lösning till stattionär värmeledning.

a_sol = solveq(K,F);

ed = extract(edof,a_sol);
figure(1)
h=patch(ex',ey',ed');
hold on
q=patch(ex',-ey',ed')
if Tinf==-96
    title('Stationary heat  distribution during night');
else
    title('Stationary heat distribution during day');
end
c = colorbar;
xlabel('x (m)');
ylabel('y (m)');
ylabel(c,'^{\circ}C','FontSize',18);

set(gca,'fontsize');
% set(h,'EdgeColor','none')
% set(q,'EdgeColor','none')

colormap(hot);
hold off;

%% transient temperatur

close all
a = a_sol;
time  = 1440;
dt = 1;

avt = zeros(length(a),length(time));
avt(:,1) = a;
for i = 2:time
    a =  (C+K*dt)\(Finit*dt+C*a);
    avt(:,i) = a;
    
    
end

times=[10,12,15,20];  % time step is in minutes, so  plots corresponds...
% transient solutions after 600s,720s,900s,1200s
timenbr = 1;
for i = 1:length(times)
    
    ed = extract(edof,avt(:,times(i)));
    %     timenbr = times(i); %timenbr+10*i^2
    figure
    h = patch(ex',ey',ed');
    hold on
    z = patch(ex',-ey',ed');
    
    title('Transient heat flow','FontSize', 20);
    c = colorbar;
    xlabel('x (m)','FontSize',18);
    ylabel('y (m)','FontSize',18);
    ylabel(c,'^{\circ}C','FontSize',18);
    set(gca,'fontsize',14)
    colormap(hot);
    if Tinf==40
        caxis([-10, 32]);
    else
        caxis([-46,-21]);
    end
    set(h,'EdgeColor','none');
    set(z,'EdgeColor','none');
    
end
%
% figure
% patch(ex',ey',ed');
% hold on
% patch(ex',-ey',ed')
%
% title('Transient heat flow','FontSize', 20);
% xlabel('x (m)','FontSize',18);
% ylabel('y (m)','FontSize',18);
% ylabel(c,'^{\circ}C','FontSize',18);
% set(gca,'fontsize',14)
% colorbar;
% colormap(jet);
% set(h,'EdgeColor','none')

%% von miesse
%konstanter
ptype = 2;
tcond =20;
E_ti = 110E9;
E_gl = 67E9 ;
v_ti = 0.34;
v_gl = 0.2;
alpha_ti = 9.4E-6;
alpha_gl = 7E-6;
thickness = 0.01;

uyo = [9 10 ];
uxo=7;


D_gl = E_gl / ((1+v_gl)*(1-2*v_gl)) .* [1-v_gl, v_gl, 0;
    v_gl, 1-v_gl, 0;
    0    0      0.5*(1-2*v_gl)];

D_ti = E_ti / ((1+v_ti)*(1-2*v_ti)) .* [1-v_ti, v_ti, 0;
    v_ti, 1-v_ti, 0;
    0    0      0.5*(1-2*v_ti)];

ep = [ ptype , thickness];

dof_S = [(1:nnod)',(nnod+1:2*nnod)']; % give each dof a number
for ie = 1:nelm
    edof_S(ie,:) = [ie dof_S(enod(ie,1),:), dof_S(enod(ie,2),:),dof_S(enod(ie,3),:)];
    
end

%kod

Ks = zeros(2*nnod);
Fs = zeros(2*nnod,1);
temperature = a_sol  ;

delta_t_enod = zeros(nelm,3);

for i = 1:nelm
    delta_t_enod(i,:) = [temperature(t(1,i)) temperature(t(2,i)) temperature(t(3,i))]; % delta_t för varje element
    delta_t = max(delta_t_enod(i,:))-tcond;
    if t(4,i) == tiarea
        De = D_ti;
        D_eps_delta_t_ti = ((E_ti*alpha_ti*delta_t)/(1-2*v_ti))*[1;1;0];
        fue = plantf(ex(i,:),ey(i,:),ep,D_eps_delta_t_ti');
    else
         De = D_gl;
        D_eps_delta_t_gl = (E_gl/(1-2*v_gl))*alpha_gl*delta_t.*[1;1;0];
        fue = plantf(ex(i,:),ey(i,:),ep,D_eps_delta_t_gl');
    end
    
    Keu = plante(ex(i,:),ey(i,:),ep,De);
 
    indx = edof_S(i,2:end);
    Ks(indx,indx) = Ks(indx,indx)+Keu;
    Fs(indx) = Fs(indx) + fue;
%   [Ks, Fs] = assem(edof_S(i,:),Ks,Keu, Fs, fue);
end




%nu fixar vi begynnelsevilkoret
bc = [];


for i = 1:length(e)
    
    if ismember(e(5,i),uyo)    
        
        bc = [bc;e(2,i)+nnod , 0];
        bc = [bc;e(1,i)+nnod , 0 ];
    end
    if(e(5,i)  ==  uxo)
        
        bc  = [bc; e(1,i),0];
        bc =  [bc; e(2,i),0];
        
    end
end



%beräkning av vonmises
F = zeros(nelm,1);
u  = solveq(Ks,Fs,bc);
ed = extract(edof_S,u);
vonmises_e = zeros(nelm,1); %vonmises för varje element

sigma = [];
for i = 1:nelm
     
     delta_t=max(delta_t_enod(i,:))-tcond;
    if t(4,i) ==  tiarea
       [sig,eps] = plants(ex(i,:),ey(i,:),ep,D_ti,ed(i,:)) ;
        eps_delta_t_ti = (1+v_ti)*alpha_ti*delta_t*[1 1 0];
        
        sigma=D_ti*(eps-eps_delta_t_ti)';
        
        sigma_zz = v_ti*(sigma(1)+sigma(2))-alpha_ti*E_ti*delta_t;
  
    else
        [sig,eps] = plants(ex(i,:),ey(i,:),ep,D_gl,ed(i,:));
         eps_delta_t_gl =(1+v_gl)*alpha_gl*delta_t*[1 1 0];
         sigma = D_gl*(eps-eps_delta_t_gl)';
        sigma_zz = v_gl*(sigma(1)+sigma(2))-alpha_gl*E_gl*delta_t;
        
    end
    es = sigma;
    vonmises_e(i) = sqrt((es(1)^2 + es(2)^2+ sigma_zz^2 -(es(1)*es(2)) ...
        -(es(1)*sigma_zz) - (es(2)*sigma_zz) + 3*es(3)^2));
end


Vonmises_n = zeros(ndof,1);

for i = 1:nnod
    [r,~] =  find(edof(:,2:4) == i);
    
    
    vonmises_n(i) = sum(vonmises_e(r)/length(r));
    
end


ed = extract(edof,vonmises_n);
figure
% h = fill(ex',ey',ed');
h = patch(ex',ey',ed');
hold on
% z = patch(ex',-ey',ed');

set(h,'EdgeColor','none')
% set(z,'EdgeColor','none')
% set(gca,'fontsize',14)
title('Von Mises stress field','FontSize', 18);
c = colorbar;
xlabel('x (m)','FontSize',18);
ylabel('y (m)','FontSize',18);
ylabel(c,'Pa','FontSize',18);


%%

% Calculate displaced coordinates
ed = extract(edof_S,u);
mag = 100; % Magnification (due to small deformations)
exd = ex + mag*ed(:,1:2:end);
eyd = ey + mag*ed(:,2:2:end);
meyd = -ey + mag*ed(:,2:2:end);
figure()
patch(ex',ey',[0 0 0],'EdgeColor','none','FaceAlpha',0.3)
hold on
patch(ex',-ey',[0 0 0],'EdgeColor','none','FaceAlpha',0.3)
hold on
patch(exd',eyd',[0 0 0],'FaceAlpha',0.3)
hold on
patch(exd',meyd',[0 0 0],'FaceAlpha',0.3)
axis equal
title('Displacement field [Magnitude enhancement 100]')



%% frontlens displacement
ed = extract(edof_S,u);
squaresum = 0;
for i = 1:nelm
    
    if t(4,i) == 3
        intNTN = plantmlmod(ex(i,:),ey(i,:));
        squaresum = squaresum+ed(i,:)*intNTN*ed(i,:)';
        
    end
end







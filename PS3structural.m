%PS 3 Structural

clear;clc;
Data = load('-mat', 'cleandata');
vars = fieldnames(Data);
for i = 1:length(vars)
    assignin('base', vars{i}, Data.(vars{i}));
end
%load probabilities estimated using NFXP
p1nfx=load('p1nfx.txt');%('NFXP1.txt')';
% first part identical to PS 2
% (2) Compute the transition probabilities
N = size(xt,2);
T = size(xt,1);

%create dummies for stay in same state (0), move one state, two states and
%replace (drep)
dummy = zeros(size(xt));
d0 = zeros(size(xt));
d1 = zeros(size(xt));
d2 = zeros(size(xt));
drep = zeros(size(xt));

%create a discrete state vector with 90 sattes corrsponding to each 5000
%miles 
gridx = 0:5000:89*5000;
disEExt = discretize(EExt,gridx);
n=size(gridx,2);
no_obs=double(isnan(EExt));
dEExt = diff(disEExt);

TT=ones(T,N);
for i=1:T
     for j=1:N
         if isnan(EExt(i,j));
             TT(i,j)=0;
         end
     end
end       
TTT=sum(sum(TT));


for i = 1:T-1
    for j = 1:N
        if dEExt(i,j) == 0 
            dummy(i,j) = 0; 
        elseif dEExt(i,j) == 1
            dummy(i,j) = 1;
        elseif dEExt(i,j) == 2
            dummy(i,j) = 2;
        elseif dEExt(i,j) < 0
            dummy(i,j) = 3;  %replacement
        else
            dummy(i,j) = 10;
        end
    end
end

for i = 1:T
    for j = 1:N
        if dummy(i,j) == 0
            d0(i,j) = 1;
        elseif dummy(i,j) == 1
            d1(i,j) = 1;
        elseif dummy(i,j) == 2
            d2(i,j) = 1;
        elseif dummy(i,j) == 3
            drep(i,j) = 1; 
        end
    end
end

repl=sum(sum(drep));

%get transition probs by summing all 3 possible transitions and dividing by
%totl amount of observations
phi0 = sum(sum(d0))/(TTT-repl);
phi1 = sum(sum(d1))/(TTT-repl);
phi2 = sum(sum(d2))/(TTT-repl);
phi = [phi0 phi1 phi2];

% (3) create transition matrix
F0 = zeros(90,90);
for i = 1:90
    if i < 89
        F0(i,i:i+2) = phi;
    elseif i == 89
        F0(i,i:i+1) = [phi0 1-phi0];
    else
        F0(i,i) = 1;
    end
end

%% calculate realtive frequencies
beta=0.99;
w=1;
for i = 1:T
    for j = 1:N
       if no_obs(i,j) == 0
           XT(w)=disEExt(i,j);
           y(w)=drep(i,j);
           w=w+1;
       end        
    end
end
NN=w-1;
repmil=XT.*y;

p1=zeros(1,90);
for i=1:90
    
    ppp=find(repmil==i);
    ttt=find(XT==i);
    pp=size(ppp,2);
    total=size(ttt,2);
    if total==0
        total=1;
    end
    p1(i)=pp/total;
    %for working, but that cant be rigth 
    if p1(i)==0
        p1(i)=0.0000000001;
    end
end

% get expression for the last part of v diff

pdiff=zeros(NN,1);
pdiffnfx=zeros(NN,1);
for xc = 1:NN
    
           lnp1=log(p1);
           lnp1nfx=log(p1nfx)';
           sump=lnp1.*F0(XT(xc),:);
           sumpnfx=lnp1nfx.*F0(XT(xc),:);
           pdiff(xc)=beta*(log(p1(1))-sum(sump));
           pdiffnfx(xc)=beta*(log(p1nfx(1))-sum(sumpnfx));
                     
end

pdiff=pdiff';
pdiffnfx=pdiffnfx';
%% Estimate the loge likelihood
loglik=@(theta) -sum((1-y).*log((exp(theta(1)-theta(2)*XT+pdiff))./(1+(exp(theta(1)-theta(2)*XT+pdiff))))+y.*log(1-(exp(theta(1)-theta(2)*XT+pdiff))./(1+(exp(theta(1)-theta(2)*XT+pdiff)))));
%loglik=@(theta) -sum(y.*log(1./(1+(exp(theta(1)-theta(2)*XT+pdiff))))+(1-y).*log(1-(1./(1+(exp(theta(1)-theta(2)*XT+pdiff))))));
logliknfx=@(thetanfx) -sum((1-y).*log((exp(thetanfx(1)-thetanfx(2)*XT+pdiffnfx))./(1+(exp(thetanfx(1)-thetanfx(2)*XT+pdiffnfx))))+y.*log(1-(exp(thetanfx(1)-thetanfx(2)*XT+pdiffnfx))./(1+(exp(thetanfx(1)-thetanfx(2)*XT+pdiffnfx)))));

%theta(1) is R 
theta0=[0,0];
thetanfx0=[0,0];

theta=fminunc(loglik,theta0);
thetanfx=fminunc(logliknfx,thetanfx0);


%% Nested Pseude likelihood

%use thatas obtained from first part as initial guess 
nn=90;
maxiter=1000;
tol=0.001;
gridx=gridx./5000;
theta=thetanfx;   %turn of if you want to use the unstable thetas calculated with the frequencies
for nx=1:6
    R=theta(1);
    M=theta(2);
    
    V0=zeros(nn,1);
    V1=zeros(nn,1);
    v0=zeros(n,1);
    v1=zeros(n,1);
    for s=1:maxiter
        
        vv=log(exp(v0)+exp(v1));
        
        for xc=1:nn
            v=vv.*F0(xc,:)';
            
            V1(xc)=-R+beta*(vv(1));
            V0(xc)=-0.001*M.*gridx(xc)+beta*sum(v);
            
        end
        %convergence criterion
        if double(all(abs(V0-v0)<tol)) > 0  && double(all(abs(V1-v1)<tol)) > 0
            display('mambo no. 5')
            break
        end
        
        if s==990;
            display('pause')
            break
        end
        %update value functions
        v0=V0;
        v1=V1;
        
    end
    
    
    
    
    P1=zeros(90,1);
    for j=1:90
        P1(j)=1-exp(theta(1)-theta(2)*gridx(j)+beta*(V0(j)-V1(j)))./(1+(exp(theta(1)-theta(2)*gridx(j)+beta*(V0(j)-V1(j)))));
        if P1(j)==0;
            P1(j)=0.0000000000000000001;
        end
    end
    
    pdiff=zeros(NN,1);
    pdiffnfx=zeros(NN,1);
    for xc = 1:NN
        
        lnp1=log(P1);
        sump=lnp1.*F0(XT(xc),:)';
        pdiff(xc)=beta*(log(P1(1))-sum(sump));
        
    end
    
    pdiff=pdiff';
    
    loglik=@(theta) -sum((1-y).*log((exp(theta(1)-theta(2)*XT+pdiff))./(1+(exp(theta(1)-theta(2)*XT+pdiff))))+y.*log(1-(exp(theta(1)-theta(2)*XT+pdiff))./(1+(exp(theta(1)-theta(2)*XT+pdiff)))));
    
    theta=fminunc(loglik,theta0);
end

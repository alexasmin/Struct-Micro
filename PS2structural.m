% P2 structural 
%Instead of creating function files for the value function iteration and
%the loglikelihood I wrote everything in one mfile. 
%% First step estimation of transition probabilities 
clear;clc;
Data = load('-mat', 'cleandata');
vars = fieldnames(Data);
for i = 1:length(vars)
    assignin('base', vars{i}, Data.(vars{i}));
end

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


%% value function iteration 
tol=0.001;
maxiter=1000;
maxiterout=100;
beta=0.99;
% first guesses 

v0=zeros(n,1);%-ones(n,1)*10;
v1=zeros(n,1);%-ones(n,1)*10;
tauMgues=1;
tauRgues=1;

% run the whole PS in one mfile
%[ V0, V1, s ] = valuefunction(tauRgues, tauMgues, beta, v0, v1, F0, gridx, tol, maxiter);

%% Outer loop
ome=0.3; %updating weight
nn=size(gridx,2);
R=tauRgues;
M=tauMgues;
gridx=gridx./5000;
for q=1:maxiterout
    %updating the thetas
    if q>1
        R=ome*R+(1-ome)*theta(1);
        M=ome*M+(1-ome)*theta(2);
    end

%% value function iteration 
%expacted values of both value fucntions 
V0=zeros(nn,1);
V1=zeros(nn,1);
gam=double(eulergamma);
v0=zeros(n,1);%-ones(n,1)*10;
v1=zeros(n,1);%-ones(n,1)*10;
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
    end
    %update value functions
    v0=V0;
    v1=V1;
    
end

%% log likelhood 

%create one large vactor over all observations for the difference between
%the value functions, the miles and the decisions
w=1;
for i = 1:T
    for j = 1:N
       if no_obs(i,j) == 0
           vdiff(w)=beta*(V0(disEExt(i,j))-V1(disEExt(i,j)));
           XT(w)=disEExt(i,j);
           y(w)=drep(i,j);
           w=w+1;
       end
        
        NN=w-1;
    end
end

%% Estimate the loge likelihood
loglik=@(theta) -sum((1-y).*log((exp(theta(1)-theta(2)*XT+vdiff))./(1+(exp(theta(1)-theta(2)*XT+vdiff))))+y.*log(1-(exp(theta(1)-theta(2)*XT+vdiff))./(1+(exp(theta(1)-theta(2)*XT+vdiff)))));

%theta(1) is R 
theta0=[R,M];

theta=fminunc(loglik,theta0);

%stopping creterion
if abs(R-theta(1))<tol && abs(M-theta(2))<tol
    display('dope')
    break
end


end

p1nfx=zeros(90,1);
for j=1:90
    p1nfx(j)=1-exp(theta(1)-theta(2)*gridx(j)+beta*(V0(j)-V1(j)))./(1+(exp(theta(1)-theta(2)*gridx(j)+beta*(V0(j)-V1(j)))));
end

save p1nfx.txt p1nfx -ascii -double 
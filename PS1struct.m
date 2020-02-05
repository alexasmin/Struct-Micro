% PS 1 Structural 
%Exercise 1 in do file 

%% EX 2

clc;

relwag=zeros(140,1);
relL=zeros(140,1);
s=1;

for i=1:140
    relwag(i)=(incwage(s+1)/uhrsworkt(s+1))/(incwage(s)/uhrsworkt(s));
    relL(i)=(uhrsworkt(s+1)/uhrsworkt(s));
    s=s+2;
end

lnrelw=log(relwag);
lnL=log(relL);

X=[ones(140,1) lnL];

theta=OLS(X,lnrelw);


lam=1/(1+exp(theta(1)));
phi=theta(2)+1;

%%
Lijt=zeros(140,1);
s=1;
for i=1:140
    Lijt(i)=(lam*(uhrsworkt(s)*52)^phi+(1-lam)*(uhrsworkt(s+1)*52)^phi)^(1/phi);
    s=s+2;
end

%%
exp1=[];
exp2=[];
exp3=[];
exp4=[];
exp5=[];
exp6=[];
exp7=[];
exp11=[];
exp21=[];
exp31=[];
exp41=[];
exp51=[];
exp61=[];
exp71=[];
dumedu=[];
s=1;
for i=1:2:280
    if experience(i)==1 && edu(i)==0
        exp1(s)=1;
    else
        exp1(s)=0;
    end
    if experience(i)==2 && edu(i)==0
        exp2(s)=1;
    else
        exp2(s)=0;
    end
    if experience(i)==3 && edu(i)==0
        exp3(s)=1;
    else
        exp3(s)=0;
    end
    if experience(i)==4 && edu(i)==0
        exp4(s)=1;
    else
        exp4(s)=0;
    end
    if experience(i)==5 && edu(i)==0
        exp5(s)=1;
    else
        exp5(s)=0;
    end
    
    if experience(i)==6 && edu(i)==0
        exp6(s)=1;
    else
        exp6(s)=0;
    end
    if experience(i)==7 && edu(i)==0
        exp7(s)=1;
    else
        exp7(s)=0;
    end
    if experience(i)==1 && edu(i)==1
        exp11(s)=1;
    else
        exp11(s)=0;
    end
    if experience(i)==2 && edu(i)==1
        exp21(s)=1;
    else
        exp21(s)=0;
    end
    if experience(i)==3 && edu(i)==1
        exp31(s)=1;
    else
        exp31(s)=0;
    end
    if experience(i)==4 && edu(i)==1
        exp41(s)=1;
    else
        exp41(s)=0;
    end
    if experience(i)==5 && edu(i)==1
        exp51(s)=1;
    else
        exp51(s)=0;
    end
    
    if experience(i)==6 && edu(i)==1
        exp61(s)=1;
    else
        exp61(s)=0;
    end
    if experience(i)==7 && edu(i)==1
        exp71(s)=1;
    else
        exp71(s)=0;
    end 
         
     dumedu(s)=edu(i);
     s=s+1;
        
        
end

edut9=zeros(140,1);
edut10=zeros(140,1);
edut11=zeros(140,1);
edut12=zeros(140,1);
edut13=zeros(140,1);
edut14=zeros(140,1);
edut15=zeros(140,1);
edut16=zeros(140,1);
edut17=zeros(140,1);
edut18=zeros(140,1);
s=1;
for i=1:2:28
        edut9(s)=edu(i);
        s=s+1;
end
for i=29:2:56
    edut10(s)=edu(i);
    s=s+1;
end
 for i=57:2:84
        edut11(s)=edu(i);
        s=s+1;
end
for i=85:2:112
    edut12(s)=edu(i);
    s=s+1;
end   
for i=113:2:140
        edut13(s)=edu(i);
        s=s+1;
end
for i=141:2:168
    edut14(s)=edu(i);
    s=s+1;
end
for i=169:2:196
        edut15(s)=edu(i);
        s=s+1;
end
for i=197:2:224
    edut16(s)=edu(i);
    s=s+1;
end
for i=225:2:252
        edut17(s)=edu(i);
        s=s+1;
end
for i=253:2:280
    edut18(s)=edu(i);
    s=s+1;
end

%%
t9=[ones(14,1); zeros(126,1)];
t10=[zeros(14,1); ones(14,1); zeros(112,1)];
t11=[zeros(28,1); ones(14,1); zeros(98,1)];
t12=[zeros(42,1); ones(14,1); zeros(84,1)];
t13=[zeros(56,1); ones(14,1); zeros(70,1)];
t14=[zeros(70,1); ones(14,1); zeros(56,1)];
t15=[zeros(84,1); ones(14,1); zeros(42,1)];
t16=[zeros(98,1); ones(14,1); zeros(28,1)];
t17=[zeros(112,1); ones(14,1); zeros(14,1)];
t18=[zeros(126,1); ones(14,1)];






%%
s=1;

for i=1:140
    wijt(i)=(incwage(s)+incwage(s+1))/Lijt(i);
    s=s+2;
end

lnwijt=log(wijt);

%%
X=[ t10 t11 t12 t13 t14 t15 t16 t17 t18 edut10 edut11 edut12 edut13 edut14 edut15 edut16 edut17 edut18 exp1' exp2' exp3' exp4' exp5' exp6' exp7' exp11' exp21' exp31' exp41' exp51' exp61' exp71' Lijt];




theta2=OLS(X,lnwijt');

gam=exp(theta2(19:end-1));

eta=theta2(33)+1;

%% normilazation 

gamnorm1=zeros(7,1);
%x=(sum(gam(1:6))-1)/5;
gamnorm2=zeros(7,1);
for i=1:7
    gamnorm1(i)=gam(i)/sum(gam(1:7));
    gamnorm2(i)=gam(7+i)/sum(gam(8:end));
end

gamnorm=[gamnorm1; gamnorm2];
%% Lit
s=1;
A=zeros(140,1);
Lit=zeros(20,1);
for i=1:10
    B=Lijt(s:s+13).^eta;
    A(s:s+13)=B.*gamnorm;
    s=s+14;
end
s=1;
for i=1:20
    Lit(i)=sum(A(s:s+6))^eta;
    s=s+7;
end


%% Third step 


edutime0=[0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1];
edutime1=[1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0];
t9=zeros(20,1);
t9(1:2)=[1 1];
t10=zeros(20,1);
t10(3:4)=[1 1];
t11=zeros(20,1);
t11(5:6)=[1 1];
t12=zeros(20,1);
t12(7:8)=[1 1];
t13=zeros(20,1);
t13(8:9)=[1 1];
t14=zeros(20,1);
t14(10:11)=[1 1];
t15=zeros(20,1);
t15(12:13)=[1 1];
t16=zeros(20,1);
t16(14:15)=[1 1];
t17=zeros(20,1);
t17(16:17)=[1 1];
t18=zeros(20,1);
t18(18:19)=[1 1];

X=[ t10 t11 t12 t13 t14 t15 t16 t17 t18 edutime0' edutime1' Lit];
%%
s=1;
for i=1:20
    wit(i)=sum(wijt(s:s+6))/Lit(i);
    s=s+7;
end
lnwit=log(wit);   
theta3=OLS(X,lnwit');  

tetta=exp(theta3(10:11))./sum(exp(theta3(10:11)));
rho=theta3(12)+1;









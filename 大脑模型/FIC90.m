load scma
%load ROICorrelation_sub_003
load ROICorrelation_FisherZ_sub_003
sc=avama;
sc=sc-diag(diag(sc));
%FC_emp=ROICorrelation(1:90,1:90);
FC_emp=ROICorrelation_FisherZ(1:90,1:90);

ds   = 100;    % BOLD downsampling rate
f_c  = 0.25;   % Frequency cut-off to remove spurious BOLD
we=1;
wi=0.7;
iz=0.382;
wplus=1.4;
Jnmda=0.15;
ge=310;
gi=615;
de=0.16;
di=0.087;
IthrE=0.403;
IthrI=0.288;
taonmda=100;
taogaba=10;
gammat=0.641/1000;
sigmat=0.01;

%脑区个数
N=size(sc,1);
JNFIC=ones(N,1);
JJJ=diag(JNFIC);
Isubdiag = find(tril(ones(N),-1)); % Values below the diagonal.
fc=FC_emp(Isubdiag);
 % Vector of all FC values below the diagonal 
%模拟时间长度
dt=0.1;
tmax=200000; % WARNING: tmax should be larger than 50000  !
tspan=0:dt:tmax;
tzhouqi=100000;
G=0.1:0.1:6;
numWg=length(G);
Coef     =zeros(1,numWg);
Maxrate  =zeros(1,numWg);
nTrials=1; % number of trials for each G
HE=@(x) ge*(x-IthrE)./(1-exp(-de*ge*(x-IthrE)));
HI=@(x) ge*(x-IthrI)./(1-exp(-di*gi*(x-IthrI)));


neuro_act=zeros(tmax,N);
ix=1;
Coef=zeros(1,numWg);
Maxrate=zeros(1,numWg);
%模拟过程
for g = G
    
     display(sprintf('gl=%g',g))
 
 
     MaxR=0;
     cb=zeros(length(Isubdiag),1);
     for tr=1:nTrials
        
        nn=1;
        j=0;
        
        %初值
        SE=0.001*ones(N,1);
        SI=0.001*ones(N,1);
        IE=zeros(N,1);
        II=zeros(N,1);
        rE=zeros(N,1);
        rI=zeros(N,1);
        
        for zhouqi=1:20
        if zhouqi>1
        set1=ave<-0.021;
        set2=ave>-0.031;
        set0=set1+set2;
        J(set0==2)=J(set0==2)-0.01;
        J(set0<2)=J(set0<2)+0.01;           
        JJ=diag(J);
        else
        J=JNFIC;
        JJ=diag(J);
        end
        ave=0;
        for i=1:tzhouqi
            
            IE=we*iz*ones(N,1)+wplus*Jnmda*SE+g*Jnmda*sc*SE-JJ*SI;
            II=wi*iz*ones(N,1)+Jnmda*SE-SI;
            rE=feval(HE,IE);
            rI=feval(HI,II);
            SE=SE+(-SE/taonmda+(1-SE)*gammat.*rE+sigmat*randn(N,1))*dt;
            SI=SI+(-SI/taogaba+rI+sigmat*randn(N,1))*dt;
            j=j+1;
            if j==(1/dt)
                
                neuro_act(nn,:)=rE';  % Down-sampling
                nn=nn+1;
                j=0;
            end
            ave=ave+(IE-125/310);
            
        end
        ave=ave/tzhouqi;
        
        end
        nn=nn-1;
        T=nn*0.001;
        
        B = BOLD(T,neuro_act(1:nn,1)'); % B=BOLD activity, bf=Foutrier transform, f=frequency range)
        BOLD_act = zeros(length(B),N);
        BOLD_act(:,1) = B;  
        for nnew=2:N
             B = BOLD(T,neuro_act(1:nn,nnew));
             BOLD_act(:,nnew) = B;
             
        end

        % Downsampling and reordering removing the first 500ms
        bds = BOLD_act(500:ds:end,:);
        bds(isnan(bds))=0;
        Cb  = corrcoef(bds);
        cb  = cb + atanh(Cb(Isubdiag))/nTrials;
        %cb  = cb + Cb(Isubdiag)/nTrials; % Vector containing all the FC values below the diagonal 
        
        
        clear BOLD BOLD_act bds
        MaxR = MaxR + max(mean(neuro_act(end-10000:end,:)))/nTrials;


    end
cb(isnan(cb))=0;
r_c = corrcoef(cb,fc);
Coef(ix)=r_c(2);
Maxrate(ix)=MaxR;
ix=ix+1;

end
plot(G,Coef,'ko-','markerfacecolor','k');
xlabel('Coupling strength (G)','FontSize',9)
ylabel('similarity','FontSize',9)
%saveas(gcf,'julei360FIC382','jpg')
%save('moxing2julei360FIC.mat');
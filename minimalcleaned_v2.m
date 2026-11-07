clear;

%Constants and units
MeV=1/140;
GeV=1000*MeV;
%Parameters
SVZ=0;
s0=(2*GeV)^2; % scale at which we are in high energy
power=2; % L_2 norm of vector
normaval=2; % L_2 norm of vector
%constraints
M=49;
nmax=49;
multnmax=nmax-1;
lmax=11;%lmax='uniform<>15';
%Number of points/resolution
Mline=2*2;
Mcircle=Mline;
numpoints=Mcircle-1 ;
nWatsonSteps=Mline/2-1; % sequential Watsonian reinforcement steps

%Load data
P = load('pointsnear.mat').Expression1; %New point %P=[0.0818378, -0.00292256];
disp(P);
vnu = importdata(join(['vnuvalsM=' num2str(M) 'nmax=' num2str(nmax) 'lmax=15.mat'],""));
H = importdata(join(['hmatyifeiM=' num2str(M) 'nmax=' num2str(nmax) 'lmax=15.mat'],""));
%rescaled amplitudes
reH = importdata(join(['rescaleReHM=' num2str(M) 'nmax=' num2str(nmax) 'lmax=15.mat'],""));
imH = importdata(join(['rescaleImHM=' num2str(M) 'nmax=' num2str(nmax) 'lmax=15.mat'],""));
imH2 = importdata(join(['rescaleImH2M=' num2str(M) 'nmax=' num2str(nmax) 'lmax=15.mat'],""));
%import veczz and vecsin
veczz=importdata(join(['veczzM=' num2str(M) 'nmaxtilde=' num2str(multnmax) '.mat'],""));
vecsin=importdata(join(['vecsinM=' num2str(M) 'nmaxtilde=' num2str(multnmax) '.mat'],""));
vecsin2=importdata("speccoeffparamlowM=49nmax=48.mat");

vecsin2=permute(vecsin2, [1, 3, 2]);
dimsbspline=size(vecsin);
dimsbspline=dimsbspline(2);
%First index is basis index, second is vnu evaluation.

%functional
F=importdata( join(['FvecyifeiM=' num2str(M) 'nmax=' num2str(nmax) 'lmax=15.mat'],""));
%functional unphysical
Fun=importdata(join(['FunvecyifeiM=' num2str(M) 'nmax=' num2str(nmax) 'lmax=15.mat'],""));

dims=size(Fun);
MunphySB=dims(1);
vsunphys=(1:MunphySB)/2;

%define theoretical point
funphytheo = importdata("funphytheo.mat");

const10P1chiraloneloop = importdata("const10P1chiraloneloop.mat");

%define factor for cursive F_1^1
FF3factor=importdata(join(["FF3factorM=" num2str(M) '.mat'],""));
QCDs=importdata(join(["SVZSRat2M=" num2str(M) '.mat'],""));
%int coefficients at 2 GeV
%S
intrho01=importdata( join(["intS01at2M=" num2str(M) '.mat'],""));
intrho00=importdata( join(["intS00at2M=" num2str(M) '.mat'],""));
intrho02=importdata( join(["intS02at2M=" num2str(M) '.mat'],""));
%P
intrho1m1=importdata( join(["intP1m1at2M=" num2str(M) '.mat'],""));
intrho10=importdata(  join(["intP10at2M=" num2str(M) '.mat'],""));
intrho11=importdata(  join(["intP11at2M=" num2str(M) '.mat'],""));
%D
intrho2m2=importdata( join(["intD0m2at2M=" num2str(M) '.mat'],""));
intrho2m1=importdata(  join(["intD0m1at2M=" num2str(M) '.mat'],""));
intrho20=importdata(  join(["intD00at2M=" num2str(M) '.mat'],""));

%Pinched integral kernels
pinchedintrho0=importdata( join(["pinchedintS00at2M=" num2str(M) '.mat'],""));
pinchedintrho1=importdata( join(["pinchedintP10at2M=" num2str(M) '.mat'],""));
%pinched integral values
pinchedQCD=importdata(  join(["pinchedSVZSR" '.mat'],""));


Amplcoef=importdata("Amplcoefnmax=49.mat");

%Define parameters
dims=size(vnu);
Mvnu=dims(1);
Mreg= 10;


Areg=1*10^(1);
Breg=1*10^(1);
Creg=1*10^(1);

anReg=2*10^(2);
bnReg=2*10^(2);
cnReg=2*10^(4);

annReg=2*10^(2);
bnnReg=2*10^(2);
cnnReg=2*10^(4);

MtolCSB=2*10^(-3);

MSVZ0=1*10^(-7);
MFF0=3*10^(-8);

MSVZ=4*10^(-6);
MFF=2*10^(-6);

MSVZ2=6*10^(-6);
MFF2=4*10^(-2);


disp(["M" "nmax" "lmax" "MtolCSB" "MSVZ1" "MFF1" "MSVZ0" "MFF0" "MSVZ2" "MFF2" "nmax FF" "Mreg" ; M nmax lmax MtolCSB MSVZ MFF MSVZ0 MFF0 MSVZ2 MFF2  multnmax Mreg*10^power]);
indexes=find(vnu>s0);
is0=indexes(1)-1

fstore=MtolCSB;

mpi=140/1000;
sqrt(vnu(is0))*mpi
vnu(is0)


dims=size(H);
i = dims(1); l = dims(2);m=dims(3); v=dims(4);


F00=1:v;
F00(:)=F(1,1,:);
F11=1:v;
F11(:)=F(2,1,:);
Funcomplement=Fun;

dims=size(Fun);

MunphySB=dims(1);



%define ratios CSB
rS0P1chi=3*(2*vsunphys-1)./(vsunphys-4);
rS2P1chi= 3*(2-vsunphys)./(vsunphys-4);



fS0unphystheo=funphytheo(1,:);
fS2unphystheo=funphytheo(2,:);
fP1unphystheo=funphytheo(3,:);


f2theo = importdata("f2theo.mat");
f2S0theo=0.0818378;
f2S2theo=transpose(f2theo(2,:));
f2P1theo=-0.00292256;


%define functional unphysical region coefficients
Fun00=zeros(MunphySB,v);
Fun00(:,:)=Fun(:,1,:);


Fun11=zeros(MunphySB,v);
Fun11(:,:)=Fun(:,3,:);

Fun02=zeros(MunphySB,v);
Fun02(:,:)=Fun(:,2,:);

Fun20=zeros(MunphySB,v);
Fun20(:,:)=Funcomplement(:,1,:);

Fun22=zeros(MunphySB,v);
Fun22(:,:)=Funcomplement(:,2,:);

Fun13=zeros(MunphySB,v);
Fun13(:,:)=Funcomplement(:,3,:);


%sqrt(4*pi/3)*(1/(8*pi^3))*


cons1=FF3factor(2,:).';
cons0=FF3factor(1,:).';
cons2=FF3factor(3,:).';







%define QCD sum rules values
%equation 2.17g
Nf=2;


QCD00=QCDs(1);
QCD01=QCDs(2);
QCD02=QCDs(3);
QCD1m1=QCDs(4);
QCD10=QCDs(5);
QCD11=QCDs(6);
QCD2m2=QCDs(7);
QCD2m1=QCDs(8);
QCD20=QCDs(9);



% define vectors for storing maximization results
sAA = zeros(Mline,5);
sAAsym= zeros(Mline,1);
% define vectors for storing f0 rho sigma
sfrs=zeros(v,Mline);
%define vector to compare with minimum distance to P



sFormfactor0=zeros(m,1);
sFormfactor1=zeros(m,1);
sFormfactor2=zeros(m,1);



sspectral0=zeros(m,1);
sspectral1=zeros(m,1);
sspectral2=zeros(m,1);

sFormfactorall0=zeros(Mline,m,1);
sFormfactorall1=zeros(Mline,m,1);
sFormfactorall2=zeros(Mline,m,1);



sspectralall0=zeros(Mline,m,1);
sspectralall1=zeros(Mline,m,1);
sspectralall2=zeros(Mline,m,1);

smatrix=zeros(6,m);

phase1=zeros(m,1);
phase2=zeros(m,1);
phase3=zeros(m,1);
phase4=zeros(m,1);
phase5=zeros(m,1);
phase6=zeros(m,1);


elasticity1=zeros(m,1);
elasticity2=zeros(m,1);
elasticity3=zeros(m,1);
elasticity4=zeros(m,1);
elasticity5=zeros(m,1);
elasticity6=zeros(m,1);

% Curves for the two extremal solutions (max first, min second).
% Channel order: S0, P1, D0, S2, F1, D2.
channelNames = {'S0','P1','D0','S2','F1','D2'};
partialWaveAll = complex(zeros(2,6,m));
phaseAll = zeros(2,6,m);
elasticityAll = zeros(2,6,m);

% The form-factor positivity matrices exist for S0, P1 and S2.
ffChannelNames = {'S0','P1','S2'};
formFactorAbs2All = zeros(2,3,m);
spectralAll = zeros(2,3,m);
ffSpectralRatioAll = nan(2,3,m);






%Start bootstrapping
for minmax=0:1



    cvx_clear
    cvx_begin  
    cvx_solver mosek

    variable f0
    variable rho1(nmax,nmax)
    variable rho2(nmax,nmax) symmetric
    variable sigma(2*nmax)
    variable ImFl1(multnmax)
    variable rhol1(dimsbspline)
    variable ImFl0(multnmax)
    variable rhol0(dimsbspline)
    variable ImFl2(multnmax)
    variable rhol2(dimsbspline)




    x1=reshape(rho1,nmax*nmax,1);
    x2=reshape(rho2,nmax*nmax,1);
    x=transpose(cat(2,f0 ,transpose(x1),transpose(x2),transpose(sigma)));

    % functional f2
    f2S0=F00*x;
    f2P1=F11*x;
    fS0unphys=Fun00*x;
    fP1unphys=Fun11*x;
    fS2unphys=Fun02*x;

    fD0unphys=Fun20*x;
    fF1unphys=Fun13*x;
    fD2unphys=Fun22*x;


    %functinoal at symmetric point:
    Asym=Amplcoef'*x;

    %SVZ coefficients l=1 I=1
    intrho10s=intrho10.'*(squeeze(vecsin2(2,1:is0,:))*rhol1);
    intrho1m1s=intrho1m1.'*(squeeze(vecsin2(2,1:is0,:))*rhol1);
    intrho11s=intrho11.'*(squeeze(vecsin2(2,1:is0,:))*rhol1);


    %SVZ coefficients l=0 I=0
    intrho00s=intrho00.'*(squeeze(vecsin2(1,1:is0,:))*rhol0);
    intrho01s=intrho01.'*(squeeze(vecsin2(1,1:is0,:))*rhol0);
    intrho02s=intrho02.'*(squeeze(vecsin2(1,1:is0,:))*rhol0);

    %SVZ coefficients l=0 I=2
    intrho2m2s=intrho2m2.'*(squeeze(vecsin2(3,1:is0,:))*rhol2);
    intrho2m1s=intrho2m1.'*(squeeze(vecsin2(3,1:is0,:))*rhol2);
    intrho20s=intrho20.'*(squeeze(vecsin2(3,1:is0,:))*rhol2);


    if minmax==1
        minimize(Asym);
    else
        maximize(Asym);
    end
    subject to


    % regularization
    norm(x1,4) <= Mreg*10^(power);
    norm(x2,4) <= Mreg*10^(power);

    rhol0(1)<=Areg;
    rhol0(1)>=1;

    rhol1(1)<=Breg;
    rhol1(1)>=1;

    rhol2(1)<=Creg;
    rhol2(1)>=1;

    abs(rhol0(2:floor(dimsbspline/2)))<=anReg;
    abs(rhol1(3:floor(dimsbspline/2)))<=bnReg;
    abs(rhol2(3:floor(dimsbspline/2)))<=cnReg;

    abs(rhol0(floor(dimsbspline/2)+1:dimsbspline))<=annReg;
    abs(rhol1(floor(dimsbspline/2)+1:dimsbspline))<=bnnReg;
    abs(rhol2(floor(dimsbspline/2)+1:dimsbspline))<=cnnReg;

    % quadratic cone
    tensor_mult(imH(: ,1:lmax,:,:),x,[4],[1]).^2+tensor_mult(reH(: ,1:lmax,:,:),x,[4],[1]).^2 <= 2*tensor_mult(imH2(:,1:lmax,:,:),x,[4],[1])

    %chiSB
    norm([fS0unphys-rS0P1chi'.*fP1unphys; fS2unphys-rS2P1chi'.*fP1unphys])<=MtolCSB;

    %Positive definite matrix
    nh2=zeros(m,v);
    nh3=zeros(m,v);
    nh4=zeros(m,v);
    nh2(:,:)=1i*H(2,1,:,:);
    nh3(:,:)=1i*H(1,1,:,:);
    nh4(:,:)=1i*H(1,2,:,:);
    for iterpo=1:is0
        FF1=1+veczz(iterpo,:)*ImFl1;
        spec1=squeeze(vecsin2(2,iterpo,:))'*rhol1;

        FF0=1+veczz(iterpo,:)*ImFl0;
        spec0=squeeze(vecsin2(1,iterpo,:))'*rhol0;

        FF2=1+veczz(iterpo,:)*ImFl2;
        spec2=squeeze(vecsin2(3,iterpo,:))'*rhol2;


        %positivesemidefinite l=1, I=1
        [1,1+nh2(iterpo,:)*x,cons1(iterpo)*FF1;1+conj(nh2(iterpo,:))*x,1,cons1(iterpo)*conj(FF1);cons1(iterpo)*conj(FF1),cons1(iterpo)*FF1,spec1]  == hermitian_semidefinite( 3 );
        %positivesemidefinite l=0, I=0
        [1,1+nh3(iterpo,:)*x,cons0(iterpo)*FF0;1+conj(nh3(iterpo,:))*x,1,cons0(iterpo)*conj(FF0);cons0(iterpo)*conj(FF0),cons0(iterpo)*FF0,spec0]  == hermitian_semidefinite( 3 );
        %positivesemidefinite l=0, I=2
        [1,1+nh4(iterpo,:)*x,cons2(iterpo)*FF2;1+conj(nh4(iterpo,:))*x,1,cons2(iterpo)*conj(FF2);cons2(iterpo)*conj(FF2),cons2(iterpo)*FF2,spec2]   == hermitian_semidefinite( 3 );


    end

    %FESR
    norm([intrho1m1s-QCD1m1; intrho10s-QCD10; intrho11s-QCD11].' )<= MSVZ;
    norm([intrho01s-QCD01; intrho00s-QCD00; intrho02s-QCD02].')<= MSVZ0;
    norm([intrho2m2s-QCD2m2; intrho2m1s-QCD2m1; intrho20s-QCD20].')<= MSVZ2;


    %High energy limit
    real((cons0(is0:m).')'.*(1+veczz(is0:m,:)*ImFl0)).^2+imag((cons0(is0:m).')'.*(1+veczz(is0:m,:)*ImFl0)).^2<= MFF0;
    real((cons1(is0:m).')'.*(1+veczz(is0:m,:)*ImFl1)).^2+imag((cons1(is0:m).')'.*(1+veczz(is0:m,:)*ImFl1)).^2<= MFF;
    real((cons2(is0:m).')'.*(1+veczz(is0:m,:)*ImFl2)).^2+imag((cons2(is0:m).')'.*(1+veczz(is0:m,:)*ImFl2)).^2<= MFF2;


    cvx_end

    %% Save observables for plotting
    solutionIndex = minmax + 1;

    % Partial waves and S-matrix observables.  The phase convention is the
    % one used in GTBwithlowspecparmax0Asym.m: S=1+i h=eta exp(2 i delta).
    channelMap = [1 1; 2 1; 1 2; 3 1; 2 2; 3 2];
    for channelIndex = 1:6
        nhPlot = reshape(H(channelMap(channelIndex,1), ...
                           channelMap(channelIndex,2),:,:),[m,v]);
        partialWave = nhPlot*x;
        sMatrix = 1 + 1i*partialWave;
        partialWaveAll(solutionIndex,channelIndex,:) = partialWave;
        phaseAll(solutionIndex,channelIndex,:) = ...
            unwrap(angle(sMatrix))*180/(2*pi);
        elasticityAll(solutionIndex,channelIndex,:) = abs(sMatrix);
    end

    % Spectral densities and squared (normalised) form factors.  Their
    % ratio directly tests the 2x2 principal-minor inequality |F|^2<=rho.
    ffCoefficients = {ImFl0,ImFl1,ImFl2};
    rhoCoefficients = {rhol0,rhol1,rhol2};
    ffNormalisations = {cons0,cons1,cons2};
    spectralIndices = [1 2 3];
    for channelIndex = 1:3
        formFactor = ffNormalisations{channelIndex}.* ...
            (1 + veczz*ffCoefficients{channelIndex});
        spectralDensity = squeeze(vecsin2(spectralIndices(channelIndex),:,:))* ...
            rhoCoefficients{channelIndex};
        formFactorAbs2 = abs(formFactor).^2;

        formFactorAbs2All(solutionIndex,channelIndex,:) = formFactorAbs2;
        spectralAll(solutionIndex,channelIndex,:) = spectralDensity;
        validDensity = spectralDensity > 100*eps;
        ratio = nan(m,1);
        ratio(validDensity) = formFactorAbs2(validDensity)./ ...
            spectralDensity(validDensity);
        ffSpectralRatioAll(solutionIndex,channelIndex,:) = ratio;
    end


    imoldH=tensor_mult(imH(:,:,1:is0,:),x,[4],[1]);

    reoldH=tensor_mult(reH(:,:,1:is0,:),x,[4],[1]);

    oldH=reoldH+1.0j*imoldH;
    array1=abs(reoldH(1,1,:)+1.0j*imoldH(1,1,:)).* reshape((1+veczz(1:is0,:)*ImFl0)./abs(1+veczz(1:is0,:)*ImFl0),[1 1 is0]);
    array2=abs(reoldH(1,2,:)+1.0j*imoldH(1,2,:)).* reshape((1+veczz(1:is0,:)*ImFl2)./abs(1+veczz(1:is0,:)*ImFl2),[1 1 is0]);
    array3=abs(reoldH(2,1,:)+1.0j*imoldH(2,1,:)).* reshape((1+veczz(1:is0,:)*ImFl1)./abs(1+veczz(1:is0,:)*ImFl1),[1 1 is0]);

    array1 = reshape(array1, [1 1 is0]);
    array2 = reshape(array2, [1 1 is0]);
    array3 = reshape(array3, [1 1 is0]);


    oldH(1,1,:) = array1;
    oldH(1,2,:) = array2;
    oldH(2,1,:) = array3;
    reoldH=real(oldH);
    imoldH=imag(oldH);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Watsonian unitarization
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sequentially maximize the tangent functional to the partial-wave
    % unitarity inequality, following
    % GTBWatsonianunitarizationreinforcesolutions.m.  For S0, P1 and S2
    % the tangent direction is first aligned with the corresponding form
    % factor phase (Watson's theorem).
    for watsonStep=1:nWatsonSteps
        cvx_clear
        cvx_begin 
        cvx_solver mosek

        variable f0
        variable rho1(nmax,nmax)
        variable rho2(nmax,nmax) symmetric
        variable sigma(2*nmax)
        variable ImFl1(multnmax)
        variable rhol1(dimsbspline)
        variable ImFl0(multnmax)
        variable rhol0(dimsbspline)
        variable ImFl2(multnmax)
        variable rhol2(dimsbspline)

        x1=reshape(rho1,nmax*nmax,1);
        x2=reshape(rho2,nmax*nmax,1);
        x=transpose(cat(2,f0,transpose(x1),transpose(x2),transpose(sigma)));

        f2S0=F00*x;
        f2P1=F11*x;
        fS0unphys=Fun00*x;
        fP1unphys=Fun11*x;
        fS2unphys=Fun02*x;
        fD0unphys=Fun20*x;
        fF1unphys=Fun13*x;
        fD2unphys=Fun22*x;
        Asym=Amplcoef'*x;

        intrho10s=intrho10.'*(squeeze(vecsin2(2,1:is0,:))*rhol1);
        intrho1m1s=intrho1m1.'*(squeeze(vecsin2(2,1:is0,:))*rhol1);
        intrho11s=intrho11.'*(squeeze(vecsin2(2,1:is0,:))*rhol1);
        intrho00s=intrho00.'*(squeeze(vecsin2(1,1:is0,:))*rhol0);
        intrho01s=intrho01.'*(squeeze(vecsin2(1,1:is0,:))*rhol0);
        intrho02s=intrho02.'*(squeeze(vecsin2(1,1:is0,:))*rhol0);
        intrho2m2s=intrho2m2.'*(squeeze(vecsin2(3,1:is0,:))*rhol2);
        intrho2m1s=intrho2m1.'*(squeeze(vecsin2(3,1:is0,:))*rhol2);
        intrho20s=intrho20.'*(squeeze(vecsin2(3,1:is0,:))*rhol2);

        trialImH=tensor_mult(imH(:,1:3,1:is0,:),x,[4],[1]);
        trialReH=tensor_mult(reH(:,1:3,1:is0,:),x,[4],[1]);
        trialImH2=tensor_mult(imH2(:,1:3,1:is0,:),x,[4],[1]);
        maximize(sum(sum(sum(imoldH(:,1:3,:).*trialImH + ...
            reoldH(:,1:3,:).*trialReH - trialImH2))));

        subject to
        norm(x1,4) <= Mreg*10^(power);
        norm(x2,4) <= Mreg*10^(power);

        rhol0(1)<=Areg;
        rhol0(1)>=1;
        rhol1(1)<=Breg;
        rhol1(1)>=1;
        rhol2(1)<=Creg;
        rhol2(1)>=1;
        abs(rhol0(2:floor(dimsbspline/2)))<=anReg;
        abs(rhol1(3:floor(dimsbspline/2)))<=bnReg;
        abs(rhol2(3:floor(dimsbspline/2)))<=cnReg;
        abs(rhol0(floor(dimsbspline/2)+1:dimsbspline))<=annReg;
        abs(rhol1(floor(dimsbspline/2)+1:dimsbspline))<=bnnReg;
        abs(rhol2(floor(dimsbspline/2)+1:dimsbspline))<=cnnReg;

        tensor_mult(imH(:,1:lmax,:,:),x,[4],[1]).^2 + ...
            tensor_mult(reH(:,1:lmax,:,:),x,[4],[1]).^2 <= ...
            2*tensor_mult(imH2(:,1:lmax,:,:),x,[4],[1]);

        norm([fS0unphys-rS0P1chi'.*fP1unphys; ...
              fS2unphys-rS2P1chi'.*fP1unphys])<=MtolCSB;

        nh2=zeros(m,v);
        nh3=zeros(m,v);
        nh4=zeros(m,v);
        nh2(:,:)=1i*H(2,1,:,:);
        nh3(:,:)=1i*H(1,1,:,:);
        nh4(:,:)=1i*H(1,2,:,:);
        for iterpo=1:is0
            FF1=1+veczz(iterpo,:)*ImFl1;
            spec1=squeeze(vecsin2(2,iterpo,:))'*rhol1;
            FF0=1+veczz(iterpo,:)*ImFl0;
            spec0=squeeze(vecsin2(1,iterpo,:))'*rhol0;
            FF2=1+veczz(iterpo,:)*ImFl2;
            spec2=squeeze(vecsin2(3,iterpo,:))'*rhol2;

            [1,1+nh2(iterpo,:)*x,cons1(iterpo)*FF1; ...
             1+conj(nh2(iterpo,:))*x,1,cons1(iterpo)*conj(FF1); ...
             cons1(iterpo)*conj(FF1),cons1(iterpo)*FF1,spec1] ...
                == hermitian_semidefinite(3);
            [1,1+nh3(iterpo,:)*x,cons0(iterpo)*FF0; ...
             1+conj(nh3(iterpo,:))*x,1,cons0(iterpo)*conj(FF0); ...
             cons0(iterpo)*conj(FF0),cons0(iterpo)*FF0,spec0] ...
                == hermitian_semidefinite(3);
            [1,1+nh4(iterpo,:)*x,cons2(iterpo)*FF2; ...
             1+conj(nh4(iterpo,:))*x,1,cons2(iterpo)*conj(FF2); ...
             cons2(iterpo)*conj(FF2),cons2(iterpo)*FF2,spec2] ...
                == hermitian_semidefinite(3);
        end

        norm([intrho1m1s-QCD1m1;intrho10s-QCD10;intrho11s-QCD11].')<=MSVZ;
        norm([intrho01s-QCD01;intrho00s-QCD00;intrho02s-QCD02].')<=MSVZ0;
        norm([intrho2m2s-QCD2m2;intrho2m1s-QCD2m1;intrho20s-QCD20].')<=MSVZ2;

        real(cons0(is0:m).*(1+veczz(is0:m,:)*ImFl0)).^2 + ...
            imag(cons0(is0:m).*(1+veczz(is0:m,:)*ImFl0)).^2<=MFF0;
        real(cons1(is0:m).*(1+veczz(is0:m,:)*ImFl1)).^2 + ...
            imag(cons1(is0:m).*(1+veczz(is0:m,:)*ImFl1)).^2<=MFF;
        real(cons2(is0:m).*(1+veczz(is0:m,:)*ImFl2)).^2 + ...
            imag(cons2(is0:m).*(1+veczz(is0:m,:)*ImFl2)).^2<=MFF2;
        cvx_end

        % Build the tangent point for the next reinforcement iteration.
        imoldH=tensor_mult(imH(:,:,1:is0,:),x,[4],[1]);
        reoldH=tensor_mult(reH(:,:,1:is0,:),x,[4],[1]);
        oldH=reoldH+1.0j*imoldH;
        watsonFF0=1+veczz(1:is0,:)*ImFl0;
        watsonFF1=1+veczz(1:is0,:)*ImFl1;
        watsonFF2=1+veczz(1:is0,:)*ImFl2;
        oldH(1,1,:)=abs(oldH(1,1,:)).*reshape(watsonFF0./abs(watsonFF0),[1 1 is0]);
        oldH(1,2,:)=abs(oldH(1,2,:)).*reshape(watsonFF2./abs(watsonFF2),[1 1 is0]);
        oldH(2,1,:)=abs(oldH(2,1,:)).*reshape(watsonFF1./abs(watsonFF1),[1 1 is0]);
        reoldH=real(oldH);
        imoldH=imag(oldH);
    end

    % Replace the pre-Watson curves by the reinforced solution.
    for channelIndex=1:6
        nhPlot=reshape(H(channelMap(channelIndex,1), ...
            channelMap(channelIndex,2),:,:),[m,v]);
        partialWave=nhPlot*x;
        sMatrix=1+1i*partialWave;
        partialWaveAll(solutionIndex,channelIndex,:)=partialWave;
        phaseAll(solutionIndex,channelIndex,:)= ...
            unwrap(angle(sMatrix))*180/(2*pi);
        elasticityAll(solutionIndex,channelIndex,:)=abs(sMatrix);
    end

    ffCoefficients={ImFl0,ImFl1,ImFl2};
    rhoCoefficients={rhol0,rhol1,rhol2};
    ffNormalisations={cons0,cons1,cons2};
    for channelIndex=1:3
        formFactor=ffNormalisations{channelIndex}.* ...
            (1+veczz*ffCoefficients{channelIndex});
        spectralDensity=squeeze(vecsin2(channelIndex,:,:))* ...
            rhoCoefficients{channelIndex};
        formFactorAbs2=abs(formFactor).^2;
        formFactorAbs2All(solutionIndex,channelIndex,:)=formFactorAbs2;
        spectralAll(solutionIndex,channelIndex,:)=spectralDensity;
        validDensity=spectralDensity>100*eps;
        ratio=nan(m,1);
        ratio(validDensity)=formFactorAbs2(validDensity)./spectralDensity(validDensity);
        ffSpectralRatioAll(solutionIndex,channelIndex,:)=ratio;
    end


end


%%%%%% Plot extremal solutions %%%%%%
energyGeV = sqrt(vnu(:))*mpi;
solutionNames = {'maximum','minimum'};
solutionStyles = {'-','--'};

% Phenomenological phase shifts, with the same channel/file assignment as
% GTBwithlowspecparmax0Asym.m.
phenoFiles = {'pheno.mat','phenop1.mat','phenod0.mat', ...
              'phenos2.mat','phenof1.mat','phenod2.mat'};
phenoPhase = cell(1,6);
for channelIndex = 1:6
    if isfile(phenoFiles{channelIndex})
        phenoPhase{channelIndex} = importdata(phenoFiles{channelIndex});
    else
        warning('Phenomenology file %s was not found.',phenoFiles{channelIndex});
        phenoPhase{channelIndex} = [];
    end
end

% Partial waves: real and imaginary parts for all six channels.
figure('Name','Partial waves','Color','w');
partialWaveLayout = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
for channelIndex = 1:6
    ax = nexttile(partialWaveLayout);
    hold(ax,'on');
    for solutionIndex = 1:2
        partialWave = squeeze(partialWaveAll(solutionIndex,channelIndex,:));
        plot(ax,energyGeV,real(partialWave),solutionStyles{solutionIndex}, ...
            'LineWidth',1.4);
        plot(ax,energyGeV,imag(partialWave),solutionStyles{solutionIndex}, ...
            'LineWidth',1.4);
    end
    yline(ax,0,':','HandleVisibility','off');
    grid(ax,'on'); box(ax,'off');
    title(ax,channelNames{channelIndex});
    xlabel(ax,'sqrt(s) [GeV]'); ylabel(ax,'partial wave');
    xlim(ax,[energyGeV(1) 3]);
    legend(ax,{'Re max','Im max','Re min','Im min'},'Location','best');
end
title(partialWaveLayout,'Partial waves');

% Phase shifts and elasticities.
figure('Name','Phase shifts and elasticities','Color','w');
phaseLayout = tiledlayout(2,6,'TileSpacing','compact','Padding','compact');
for channelIndex = 1:6
    ax = nexttile(phaseLayout,channelIndex);
    hold(ax,'on');
    if ~isempty(phenoPhase{channelIndex})
        plot(ax,phenoPhase{channelIndex}(:,1),phenoPhase{channelIndex}(:,2), ...
            'k-','LineWidth',1.5,'DisplayName','Pheno');
    end
    for solutionIndex = 1:2
        plot(ax,energyGeV,squeeze(phaseAll(solutionIndex,channelIndex,:)), ...
            solutionStyles{solutionIndex},'LineWidth',1.4, ...
            'DisplayName',solutionNames{solutionIndex});
    end
    grid(ax,'on'); box(ax,'off'); title(ax,channelNames{channelIndex});
    xlabel(ax,'sqrt(s) [GeV]'); ylabel(ax,'delta [deg]');
    xlim(ax,[energyGeV(1) 3]);
    legend(ax,'show','Location','best');

    ax = nexttile(phaseLayout,6+channelIndex);
    hold(ax,'on');
    for solutionIndex = 1:2
        plot(ax,energyGeV,squeeze(elasticityAll(solutionIndex,channelIndex,:)), ...
            solutionStyles{solutionIndex},'LineWidth',1.4);
    end
    yline(ax,1,':','elastic limit','HandleVisibility','off');
    grid(ax,'on'); box(ax,'off');
    xlabel(ax,'sqrt(s) [GeV]'); ylabel(ax,'eta');
    xlim(ax,[energyGeV(1) 3]);
end
title(phaseLayout,'Phase shifts (top) and elasticities (bottom)');

% Spectral density versus |form factor|^2.  Equality means saturation of
% the corresponding principal-minor positivity inequality.
figure('Name','Spectral density and form-factor saturation','Color','w');
saturationLayout = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
for channelIndex = 1:3
    ax = nexttile(saturationLayout,channelIndex);
    hold(ax,'on');
    for solutionIndex = 1:2
        plot(ax,energyGeV,squeeze(spectralAll(solutionIndex,channelIndex,:)), ...
            solutionStyles{solutionIndex},'LineWidth',1.5);
        plot(ax,energyGeV,squeeze(formFactorAbs2All(solutionIndex,channelIndex,:)), ...
            solutionStyles{solutionIndex},'LineWidth',1.1);
    end
    grid(ax,'on'); box(ax,'off'); title(ax,ffChannelNames{channelIndex});
    xlabel(ax,'sqrt(s) [GeV]'); ylabel(ax,'rho and |F|^2');
    xlim(ax,[energyGeV(1) 3]);
    legend(ax,{'rho max','|F|^2 max','rho min','|F|^2 min'},'Location','best');

    ax = nexttile(saturationLayout,3+channelIndex);
    hold(ax,'on');
    for solutionIndex = 1:2
        plot(ax,energyGeV(1:is0), ...
            squeeze(ffSpectralRatioAll(solutionIndex,channelIndex,1:is0)), ...
            solutionStyles{solutionIndex},'LineWidth',1.5);
    end
    yline(ax,1,':','saturation','HandleVisibility','off');
    ylim(ax,[0 1.05]); grid(ax,'on'); box(ax,'off');
    xlabel(ax,'sqrt(s) [GeV]'); ylabel(ax,'|F|^2/rho');
    xlim(ax,[energyGeV(1) 3]);
end
title(saturationLayout,'Form-factor/spectral-density inequality');


%%%%%Export results %%%%%



beep on; beep;

%%%%%% Definition of functions %%%%%%%%
function C = tensor_mult(A, B, sum_idx_A, sum_idx_B)

sum_idx_A = reshape(sum_idx_A, 1, numel(sum_idx_A));
sum_idx_B = reshape(sum_idx_B, 1, numel(sum_idx_B));

size_A = size(A);
size_B = size(B);

num_size_A = numel(size_A);
num_size_B = numel(size_B);

perm_A = 1:num_size_A;
perm_B = 1:num_size_B;

num_idx = size(sum_idx_A, 2);

assert(isequal(size(sum_idx_A) , size(sum_idx_B) ));
assert(isequal(size_A(sum_idx_A), size_B(sum_idx_B)));

sum_dim = prod(size_A(sum_idx_A));

for i = 1:num_idx
    perm_A = perm_A(perm_A~=sum_idx_A(i));
    perm_B = perm_B(perm_B~=sum_idx_B(i));
end

perm_A = [perm_A, sum_idx_A];
perm_B = [sum_idx_B, perm_B];

size_A(sum_idx_A) = [];
size_B(sum_idx_B) = [];

if isempty(size_A)
    size_A = 1;
end
if isempty(size_B)
    size_B = 1;
end

C = squeeze(reshape(...
    reshape(permute(A, perm_A), [prod(size_A), sum_dim]) * ...
    reshape(permute(B, perm_B), [sum_dim, prod(size_B)]), ...
    [size_A,size_B]));

end






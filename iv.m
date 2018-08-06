clear all;
  % Contains (all MKS, except energy in eV)
  hbar=1.06e-34;
  q=1.6e-19;
  m=0.25*9.1e-31;
  IE=(q*q)/(2.0*pi*hbar);
  Ef=0.1;
  kT=0.025;

% inputs
  a=3.0e-10;
  t0=hbar*hbar/(2.0*m*a*a*q);

% Hamiltonian Matrix
  NS=15;
  NC=16;
  ND=15;
  Np=NS+NC+ND;
% RT Barrier
  UB=[zeros(NS,1);0.4*ones(4,1);zeros(NC-8,1);0.4*ones(4,1);zeros(ND,1)]; %RT Barrier
%  UB=[zeros(NS,1);0.4*ones(NC,1);zeros(ND,1);];%tunneling barrier
  T=(2.0*t0*diag(ones(1,Np)))-(t0*diag(ones(1,Np-1),1))-(t0*diag(ones(1,Np-1),-1));
  T=T+diag(UB);
  %Bias
  NV=26;
  VV=linspace(0,0.1,NV);
  for iV=1:NV
    V=VV(iV);
    mu1=Ef+(V/2);
    mu2=Ef-(V/2);
    U1=V*[0.5*ones(1,NS) linspace(0.5,-0.5,NC) -0.5*ones(1,ND)];
    U1=U1';%Applied potential profile

 %Energy grid for Green’s function method
    NE=101;E=linspace(-.2,.8,NE);zplus=i*1e-12;dE=E(2)-E(1);
    f1=1./(1+exp((E-mu1)./kT));
    f2=1./(1+exp((E-mu2)./kT));
    %For infinite 2-D cross-section
    %f1=(2*m*kT*q/(2*pi*hbarˆ2)).*log(1+exp((mu1-E)./kT));
    %f2=(2*m*kT*q/(2*pi*hbarˆ2)).*log(1+exp((mu2-E)./kT));

    %Transmission

    I=0;%Current
    for k=1:NE
      sig1=zeros(Np);sig2=zeros(Np);sig3=zeros(Np);
      ck=1-((E(k)+zplus-U1(1)-UB(1))/(2*t0));ka=acos(ck);
      sig1(1,1)=-t0*exp(i*ka);
      gam1=i*(sig1-sig1');
      ck=1-((E(k)+zplus-U1(Np)-UB(Np))/(2*t0));
      ka=acos(ck);
      sig2(Np,Np)=-t0*exp(i*ka);
      gam2=i*(sig2-sig2');
      sig3(Np/2,Np/2)=-i*0.00025;
      gam3=i*(sig3-sig3'); %B¨uttiker probe
      G=inv(((E(k)+zplus)*eye(Np))-T-diag(U1)-sig1-sig2-sig3);
      T12=real(trace(gam1*G*gam2*G'));
      T13=real(trace(gam1*G*gam3*G'));
      T23=real(trace(gam2*G*gam3*G'));
      TM(k)=T12+(T13*T23/(T12+T23));
      I=I+(dE*IE*TM(k)*(f1(k)-f2(k)));
    end
    II(iV)=I;
  end
  XX=a*1.0e+09*[1:1:Np];
  XS=XX([1:NS-4]);
  XD=XX([NS+NC+5:Np]);
  hold on
  h=plot(XX,U1+UB,'b');
  %h=plot(VV,II,'b');
%  h=plot(XS,mu1*ones(1,NS-4),'b--');
%  h=plot(XD,mu2*ones(1,NS-4),'b--');
 % axis([0 .5 0 3.5e-7])
  axis([0 15 -.3 .7])
  set(h,'linewidth',[2.0]);
  set(gca,'Fontsize',[14]);
%  xlabel('Voltage (V)');
  xlabel('z (nm)');
  ylabel('Energy (eV)');
%  ylabel('Current (A)');
  grid on

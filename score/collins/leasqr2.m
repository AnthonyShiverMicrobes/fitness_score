function [f,p,kvg,iter,corp,covp,covr,stdresid,Z,r2]= leasqr2(x,y,pin,F,extraP,stol,niter,wt,dp,dFdp,options)
% LEASQR2 is a modification of LEASQR that can take an extra variable
% (extraP) of additional parameters (that are not fit).

%if (sscanf(version,'%f') >= 4),
vernum= sscanf(version,'%f');
if vernum(1) >= 4,
  global verbose
  plotcmd='plot(x(:,1),y,''+'',x(:,1),f); figure(gcf)';
else
  plotcmd='plot(x(:,1),y,''+'',x(:,1),f); shg';
end;
if (exist('OCTAVE_VERSION'))
  global verbose
  plotcmd='plot(x(:,1),y,"+;data;",x(:,1),f,";fit;");';
end;

if(exist('verbose')~=1), %If verbose undefined, print nothing
	verbose=0;       %This will not tell them the results
end;

if (nargin <= 9), dFdp='dfdp2'; end;
if (nargin <= 8), dp=.001*(pin*0+1); end; %DT
if (nargin <= 7), wt=ones(length(y),1); end;	% SMB modification
if (nargin <= 6), niter=20; end;
if (nargin == 5), stol=.0001; end;
%

y=y(:); wt=wt(:); pin=pin(:); dp=dp(:); x=x(:); %change all vectors to columns
% check data vectors- same length?
m=length(y); n=length(pin); p=pin;[m1,m2]=size(x);
if m1~=m ,error('input(x)/output(y) data must have same number of rows ') ,end;

if (nargin <= 10), 
  options=[zeros(n,1), Inf*ones(n,1)];
  nor = n; noc = 2;
else
  [nor, noc]=size(options);
  if (nor ~= n),
    error('options and parameter matrices must have same number of rows'),
  end;
  if (noc ~= 2),
    options=[options(noc,1), Inf*ones(noc,1)];
  end;
end;
pprec=options(:,1);
maxstep=options(:,2);
%

% set up for iterations
%
f=feval(F,x,p,extraP); fbest=f; pbest=p;
if (size(f) ~= size(y))
    f = f';
end
r=wt.*(y-f);
sbest=r'*r;
ss = sbest;             %I added this line
nrm=zeros(n,1);
chgprev=Inf*ones(n,1);
kvg=0;
epsLlast=1;
epstab=[.1, 1, 1e2, 1e4, 1e6];

% do iterations
%
for iter=1:niter,
  pprev=pbest;
  prt=feval(dFdp,x,fbest,pprev,dp,F,extraP);
  r=wt.*(y-fbest);
  sprev=sbest;
  sgoal=(1-stol)*sprev;
  for j=1:n,
    if dp(j)==0,
      nrm(j)=0;
    else
      prt(:,j)=wt.*prt(:,j);
      nrm(j)=prt(:,j)'*prt(:,j);
      if nrm(j)>0,
        nrm(j)=1/sqrt(nrm(j));
      end;
    end
    prt(:,j)=nrm(j)*prt(:,j);
  end;
% above loop could ? be replaced by:
% prt=prt.*wt(:,ones(1,n)); 
% nrm=dp./sqrt(diag(prt'*prt)); 
% prt=prt.*nrm(:,ones(1,m))';
  [prt,s,v]=svd(prt,0);
  s=diag(s);
  g=prt'*r;
  for jjj=1:length(epstab),
    epsL = max(epsLlast*epstab(jjj),1e-7);
    se=sqrt((s.*s)+epsL);
    gse=g./se;
    chg=((v*gse).*nrm);
%   check the change constraints and apply as necessary
    ochg=chg;
    for iii=1:n,
      if (maxstep(iii)==Inf), break; end;
      chg(iii)=max(chg(iii),-abs(maxstep(iii)*pprev(iii)));
      chg(iii)=min(chg(iii),abs(maxstep(iii)*pprev(iii)));
    end;
    if (verbose & any(ochg ~= chg)),
      disp(['Change in parameter(s): ', ...
         sprintf('%d ',find(ochg ~= chg)), 'were constrained']);
    end;
    aprec=abs(pprec.*pbest);       %---
% ss=scalar sum of squares=sum((wt.*(y-f))^2).
    if (any(abs(chg) > 0.1*aprec)),%---  % only worth evaluating function if
      p=chg+pprev;                       % there is some non-miniscule change
      f=feval(F,x,p,extraP);
      r=wt.*(y-f);
      ss=r'*r;
      if ss<sbest,
        pbest=p;
        fbest=f;
        sbest=ss;
      end;
      if ss<=sgoal,
        break;
      end;
    end;                          %---
  end;
  epsLlast = epsL;
  if (verbose),
    eval(plotcmd);
  end;
  if ss<eps,
    break;
  end
  aprec=abs(pprec.*pbest);
%  [aprec, chg, chgprev]
  if (all(abs(chg) < aprec) & all(abs(chgprev) < aprec)),
    kvg=1;
    if (verbose),
      fprintf('Parameter changes converged to specified precision\n');
    end;
    break;
  else
    chgprev=chg;
  end;
  if ss>sgoal,
    break;
  end;
end;

% set return values
%
p=pbest;
f=fbest;
ss=sbest;
kvg=((sbest>sgoal)|(sbest<=eps)|kvg);
if kvg ~= 1 , disp(' CONVERGENCE NOT ACHIEVED! '), end;

% CALC VARIANCE COV MATRIX AND CORRELATION MATRIX OF PARAMETERS
% re-evaluate the Jacobian at optimal values
jac=feval(dFdp,x,f,p,dp,F,extraP);
msk = dp ~= 0;
n = sum(msk);           % reduce n to equal number of estimated parameters
jac = jac(:, msk);	% use only fitted parameters

%% following section is Ray Muzic's estimate for covariance and correlation
%% assuming covariance of data is a diagonal matrix proportional to
%% diag(1/wt.^2).  
%% cov matrix of data est. from Bard Eq. 7-5-13, and Row 1 Table 5.1 

if vernum(1) >= 4,
  Q=sparse(1:m,1:m,(0*wt+1)./(wt.^2));  % save memory
  Qinv=inv(Q);
else
  Q=diag((0*wt+1)./(wt.^2));
  Qinv=diag(wt.*wt);
end;
resid=y-f;                                    %un-weighted residuals
covr=resid'*Qinv*resid*Q/(m-n);                 %covariance of residuals
Vy=1/(1-n/m)*covr;  % Eq. 7-13-22, Bard         %covariance of the data 

jtgjinv=inv(jac'*Qinv*jac);			%argument of inv may be singular
covp=jtgjinv*jac'*Qinv*Vy*Qinv*jac*jtgjinv; % Eq. 7-5-13, Bard %cov of parm est
d=sqrt(abs(diag(covp)));
corp=covp./(d*d');

covr=diag(covr);                 % convert returned values to compact storage
stdresid=resid./sqrt(diag(Vy));  % compute then convert for compact storage
Z=((m-n)*jac'*Qinv*jac)/(n*resid'*Qinv*resid);

%%% alt. est. of cov. mat. of parm.:(Delforge, Circulation, 82:1494-1504, 1990
%%disp('Alternate estimate of cov. of param. est.')
%%acovp=resid'*Qinv*resid/(m-n)*jtgjinv

%Calculate R^2 (Ref Draper & Smith p.46)
%
r=corrcoef(y,f);
if (exist('OCTAVE_VERSION'))
  r2=r^2;
else
  r2=r(1,2).^2;
end

% if someone has asked for it, let them have it
%
if (verbose), 
  eval(plotcmd);
  disp(' Least Squares Estimates of Parameters')
  disp(p')
  disp(' Correlation matrix of parameters estimated')
  disp(corp)
  disp(' Covariance matrix of Residuals' )
  disp(covr)
  disp(' Correlation Coefficient R^2')
  disp(r2)
  sprintf(' 95%% conf region: F(0.05)(%.0f,%.0f)>= delta_pvec''*Z*delta_pvec',n,m-n)
  Z
%   runs test according to Bard. p 201.
  n1 = sum((f-y) < 0);
  n2 = sum((f-y) > 0);
  nrun=sum(abs(diff((f-y)<0)))+1;
  if ((n1>10)&(n2>10)), % sufficent data for test?
    zed=(nrun-(2*n1*n2/(n1+n2)+1)+0.5)/(2*n1*n2*(2*n1*n2-n1-n2)...
      /((n1+n2)^2*(n1+n2-1)));
    if (zed < 0),
      prob = erfc(-zed/sqrt(2))/2*100;
      disp([num2str(prob),'% chance of fewer than ',num2str(nrun),' runs.']);
    else,
      prob = erfc(zed/sqrt(2))/2*100;
      disp([num2str(prob),'% chance of greater than ',num2str(nrun),' runs.']);
    end;
  end;
end;

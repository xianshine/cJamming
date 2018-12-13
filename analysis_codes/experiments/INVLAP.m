% INVLAP – Numerical Inversion of Laplace Transforms 
% Copyright (c) 2011, Juraj Valsa All rights reserved.

% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:

%    * Redistributions of source code must retain the above copyright 
%     notice, this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright 
%     notice, this list of conditions and the following disclaimer in 
%      the documentation and/or other materials provided with the distribution
%    * Neither the name of the BUT Brno, Czech Republic nor the names 
%      of its contributors may be used to endorse or promote products derived 
%      from this software without specific prior written permission.
      
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.

function [radt,ft]=INVLAP(Fs,tini,tend,nnt,a,ns,nd);
% Fs is formula for F(s) as a string
% tini, tend are limits of the solution interval
% nnt is total number of time instants
% a, ns, nd are parameters of the method 
% if not given, the method uses implicit values a=6, ns=20, nd=19
% it is recommended to preserve a=6
% increasing ns and nd leads to lower error
% an example of function calling  
% [t,ft]=INVLAP('s/(s^2+4*pi^2)',0,10,1001);
% to plot the graph of results write plot(t,ft), grid on, zoom on
FF=strrep(strrep(strrep(Fs,'*','.*'),'/','./'),'^','.^');
if nargin==4
  a=6; ns=20; nd=19;  end;    % implicit parameters
radt=linspace(tini,tend,nnt); % time vector
if tini==0  radt=radt(2:1:nnt);  end;  % t=0 is not allowed
tic					% measure the CPU time
for n=1:ns+1+nd               % prepare necessary coefficients
   alfa(n)=a+(n-1)*pi*j;
   beta(n)=-exp(a)*(-1)^n;
end;
n=1:nd;
bdif=fliplr(cumsum(gamma(nd+1)./gamma(nd+2-n)./gamma(n)))./2^nd;
beta(ns+2:ns+1+nd)=beta(ns+2:ns+1+nd).*bdif;
beta(1)=beta(1)/2;
for kt=1:nnt                  % cycle for time t
   tt=radt(kt);
   s=alfa/tt;                 % complex frequency s
   bt=beta/tt;
   btF=bt.*eval(FF);          % functional value F(s)
   ft(kt)=sum(real(btF));     % original f(tt)
end;
toc

function [p_n_calc SSE] =  contact(omega_max, z,r, pc, display_result,gene)
p_n_calc = zeros(length(z),length(r));
for ii = 1:length(z)
    for jj = 1:length(r)
        p_n_calc(ii,jj) =p_neighbor_contact(z(ii,1),r, omega_max, r(jj,1),pc);
    end
end
p_n_calc = sum(p_n_calc.*(ones(length(z),1)*(r(:,2))'),2)*mean(diff(r(:,1)));
SSE = sum((p_n_calc-z(:,2)).^2);

if  display_result
    plot(p_n_calc); hold on
    plot(z(:,2),'.r');hold off
    title([gene,' p=',num2str(pc)]);
    xlabel('# of contact');
    ylabel('pdf');
    legend('fit','data');
    drawnow
end
function p = p_neighbor_contact(z,r,omega_max,rc,pc)
omega = 0.1:0.1:pi*2;
f = p_angle_contact(r, rc, omega);
p = inverse_laplace_contact(f,z,omega_max,omega,pc);



% INVLAP – Numerical Inversion of Laplace Transforms 
function ft=inverse_laplace_contact(f,z,omega_max,omega, pc, a,ns,nd)
% f: angle distribution given rc
% neighbor: the number of neighbors given rc
% omega_max: maximum angle around a particle rc
% omega: range of angle
% a, ns, nd are parameters of the method 
% if not given, the method uses implicit values a=6, ns=20, nd=19
% it is recommended to preserve a=6
% increasing ns and nd leads to lower error

if nargin==5
  a=6; ns=20; nd=19;  end;    % implicit parameters

alfa = zeros(1,ns+1+nd);
beta = alfa;
for n=1:ns+1+nd               % prepare necessary coefficients
   alfa(n)=a+(n-1)*pi*1i;
   beta(n)=-exp(a)*(-1)^n;
end
n=1:nd;
bdif=fliplr(cumsum(gamma(nd+1)./gamma(nd+2-n)./gamma(n)))./2^nd;
beta(ns+2:ns+1+nd)=beta(ns+2:ns+1+nd).*bdif;
beta(1)=beta(1)/2;

s=alfa/omega_max;                 % complex frequency s
bt=beta/omega_max;
btF=bt.*laplace_angle_contact(f,s,z,omega,pc);          % functional value F(s)
ft=sum(real(btF));     % original f(tt)


function F = laplace_angle_contact(f,s,z,omega,pc)
% perform step 2 and 3 of the fitting program, calculating the laplace
% transform of the angle distribution and F for a given n, s
F = zeros(size(s));
for ii = 1:length(s)
    a = real(s(ii));
    b = imag(s(ii));
    F(ii) = mean(diff(omega))*(sum(cos(b*omega).*exp(-a*omega).*f)-1i*sum(sin(b*omega).*exp(-a*omega).*f));
end
F = (1-F)./s.*(pc*F).^z./(1-(1-pc)*F).^(z+1);


function f = p_angle_contact(r, rc, omega)
% perform step 1 of the fitting program, calculating the angle distribution
% given radius distribution and central particle radius
r = [r; 1e10 0];
f = zeros(size(omega));
for ii = 1:length(omega)
    omega_scaled = omega(ii)/2/pi;
    current_r = rc*(omega_scaled*2-omega_scaled^2+sqrt(omega_scaled*2-omega_scaled^2))/(1-omega_scaled)^2;
    index=find(r(:,1)>current_r,1);
    p = ((current_r-r(index-1,1))*r(index,2)+(r(index,1)-current_r)*r(index-1,2))/(r(index,1)-r(index-1,1));
    f(ii) = rc * (1+omega_scaled*(2-omega_scaled)+2*sqrt(omega_scaled*(2-omega_scaled)))/(1-omega_scaled)^3/sqrt(omega_scaled*(2-omega_scaled))*p/2/pi;
end





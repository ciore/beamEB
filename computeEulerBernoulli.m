% This file is part of beamEB, a code to compute the Euler-Bernoulli
% solution for various common beam problems.
% 
% Copyright (C) 2020 Ciarán O'Reilly <ciaran@kth.se>
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

function beam=computeEulerBernoulli(model,plotcurves)

if nargin<2
  plotcurves=0;
end
  
names=fieldnames(model);
for i=1:numel(names)
  eval([names{i},'=model.',names{i},';'])
end

x=0:L/1e3:L;

switch loadcase
  
  case 'simple_pt'
    
    Q = -P*(L - xP)/L*ones(size(x));
    Q(x>xP) = P*xP/L;
    M = -(P*x*(L - xP))/L;
    M(x>xP) = -(P*xP*(L - x(x>xP)))/L;
    t = -(P*(L - xP)*(3*x.^2 + xP^2 - 2*L*xP))/(6*E*I*L);
    t(x>xP) = (P*xP*(2*L^2 - 6*L*x(x>xP) + 3*x(x>xP).^2 + xP^2))/(6*E*I*L);
    w = -(P*x.*(L - xP).*(x.^2 + xP^2 - 2*L*xP))/(6*E*I*L);
    w(x>xP) = -(P*xP*(L - x(x>xP)).*(x(x>xP).^2 - 2*L*x(x>xP) + xP^2))/(6*E*I*L);
    beta = (1:5)*pi/L;
    f=beta.^2*sqrt(model.E*model.I/(model.rho*model.A))/2/pi;
     
  case 'cantilever_ptend'
    
    Q = -P*ones(size(x));
    M = P*(L - x);
    t = P*x.*(2*L-x)/(2*E*I);
    w = P*x.^2.*(3*L - x)/(6*E*I);
    betaL = [1.875104068711961e+00     4.694091132974175e+00     7.854757438237613e+00     1.099554073487547e+01];
    f=(betaL/L).^2*sqrt(model.E*model.I/(model.rho*model.A))/2/pi;

  case 'cantilever_dist'
    
    Q = -P*(L - x);
    M = P*(L^2 - 2*L*x + x.^2)/2;
    t = P*x.*(3*L^2 - 3*L*x + x.^2)/(6*E*I);
    w = P*x.^2.*(6*L^2 - 4*L*x + x.^2)/(24*E*I);
    betaL = [1.875104068711961e+00     4.694091132974175e+00     7.854757438237613e+00     1.099554073487547e+01];
    f=(betaL/L).^2*sqrt(model.E*model.I/(model.rho*model.A))/2/pi;

  otherwise
    
    warning('Case not defined')

end

if plotcurves
  subplot(2,2,1)
  plot(x,w)
  xlabel('x'),ylabel('w')
  subplot(2,2,2)
  plot(x,t)
  xlabel('x'),ylabel('t')
  subplot(2,2,3)
  plot(x,M)
  xlabel('x'),ylabel('M')
  subplot(2,2,4)
  plot(x,Q)
  xlabel('x'),ylabel('Q')
end

if exist('f')
  beam=struct('x',x,'w',w,'t',t,'M',M,'Q',Q,'f',f);
else
  beam=struct('x',x,'w',w,'t',t,'M',M,'Q',Q);
end
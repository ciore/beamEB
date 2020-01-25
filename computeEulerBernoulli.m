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

x=0:1e-3:L;

switch loadcase
  
  case 'simple_pt'
    
    Q = -F*(L - xF)/L*ones(size(x));
    Q(x>xF) = F*xF/L;
    M = -(F*x*(L - xF))/L;
    M(x>xF) = -(F*xF*(L - x(x>xF)))/L;
    t = -(F*(L - xF)*(- 2*L^2*xF + 4*L*xF^2 - 3*xF^2 + 3*x.^2))/(6*E*I*L);
    t(x>xF) = (F*xF*(2*L^3 - 6*L^2*xF + 4*L*xF^2 + 6*L*xF - 6*L*x(x>xF) - 3*xF^2 + 3*x(x>xF).^2))/(6*E*I*L);
    w = -(F*x*(L - xF).*(- 2*L^2*xF + 4*L*xF^2 - 3*xF^2 + x.^2))/(6*E*I*L);
    w(x>xF) = -(F*xF*(L - x(x>xF)).*(2*L^3 - 6*L^2*xF - 2*L^2 + 4*L*xF^2 + 6*L*xF - 2*L*x(x>xF) - 3*xF^2 + x(x>xF).^2))/(6*E*I*L);
    
  case 'cantilever_ptend'
    
    Q = -F*ones(size(x));
    M = F*(L - x);
    t = F*x.*(2*L-x)/(2*E*I);
    w = F*x.^2.*(3*L - x)/(6*E*I);

  case 'cantilever_dist'
    
    Q = -F*(L - x);
    M = F*(L^2 - 2*L*x + x.^2)/2;
    t = F*x.*(3*L^2 - 3*L*x + x.^2)/(6*E*I);
    w = F*x.^2.*(6*L^2 - 4*L*x + x.^2)/(24*E*I);

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

beam=struct('x',x,'w',w,'t',t,'M',M,'Q',Q);
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

clear

%% define the model
loadcase='simple_pt';
L=2;


switch loadcase

  case 'simple_ptc'
    
    
    
  case 'simple_pt'
    
    f = @(x,L) sin(x*L);
    df = @(x,L) L*cos(L*x);
    
    beta=[];
    neigs=5;
    x0=1;
    while numel(beta)<neigs
      x=x0;
      d=1;
      while  d>1e-10
        d=x;
        x = x - f(x,L)/df(x,L);
        d = abs(d-x)/abs(d);
      end
      if not(ismember(x,beta))
        beta = [beta x];
      end
      x0=x0+1;
    end
    
  case 'cantilever_ptend'
    
    f = @(x,L) cos(x*L) + 1/cosh(x*L);
    df = @(x,L) - L*sin(L*x) - (L*sinh(L*x))/cosh(L*x)^2;
    
    beta=[];
    neigs=5;
    x0=1;
    while numel(beta)<neigs
      x=x0;
      d=1;
      while  d>1e-10
        d=x;
        x = x - f(x,L)/df(x,L);
        d = abs(d-x)/abs(d);
      end
      if not(ismember(x,beta))
        beta = [beta x];
      end
      x0=x0+1;
    end
    
  case 'cantilever_dist'
    
    
    
  otherwise
    
    warning('Case not defined')
    
end
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

syms x P xP
syms E I L positive
casename='simple_dist';

switch casename

  case 'simple_ptc'
    
    w=[P*x.*(3*L^2-4*x.^2)/(48*E*I)
      P*(L-x).*(-L^2+8*L*x-4*x.^2)/(48*E*I)];
    t=[P*(L^2-4*x.^2)/(16*E*I)
      P*(3*L^2-8*L*x+4*x.^2)/(16*E*I)];
    M=[-P*x/2
      -P*(L-x)/2];
        
    Q=[-P/2
      P/2]
    
    M=[simplify(int(Q(1),x)-subs(int(Q(1),x),x,0))
      simplify(int(Q(2),x)-subs(int(Q(2),x),x,L))]
    
    t=[simplify(int(M(1)/E/I,x)+subs(int(M(2)/E/I,x),x,L/2))
      simplify(int(M(2)/E/I,x)+subs(int(M(1)/E/I,x),x,L/2))]
    
    w=[simplify(int(t(1),x)-subs(int(t(1),x),x,0))
      simplify(int(t(2),x)-subs(int(t(2),x),x,L))]
    
  case 'simple_pt'
    
    w=[P*(L-xP)*x.*(L^2-(L-xP)^2-x.^2)/(6*L*E*I)
     P*xP*(L-x).*(2*L*x-xP^2-x.^2)/(6*E*I*L)];
    t=[P*(L-xP)*(L^2-(L-xP)^2-3*x.^2)/(6*L*E*I)
      P*xP*(2*L^2-6*L*x+xP^2+3*x.^2)/(6*L*E*I)];
    M=[-P*(L-xP)*x/L
      -P*xP*(L-x)/L];
    
    Q=[-P*(L-xP)/L;
      P*xP/L];
    
    M=[simplify(int(Q(1),x)-subs(int(Q(1),x),x,0))
      simplify(int(Q(2),x)-subs(int(Q(2),x),x,L))]
    
    t=[simplify(int(M(1)/E/I,x)-subs(int(M(1)/E/I,x),x,xP))
      simplify(int(M(2)/E/I,x)-subs(int(M(2)/E/I,x),x,xP))];
    dt=(int(t(1),0,xP)+int(t(2),xP,L))/L;
    t=simplify(t-dt)
    
    w=[simplify(int(t(1),x)-subs(int(t(1),x),x,0))
      simplify(int(t(2),x)-subs(int(t(2),x),x,L))]

  case 'simple_dist'
    
    Q=-P/2*(L-2*x)

    M=simplify(int(Q(1),x)-subs(int(Q(1),x),x,0))

    t=simplify(int(M(1)/E/I,x)-subs(int(M(1)/E/I,x),x,L/2))

    w=simplify(int(t(1),x)-subs(int(t(1),x),x,0))
    
  case 'cantilever_ptend'
    
    w=P*x.^2.*(3*L-x)/(6*E*I);
    t=P*x.*(2*L-x)/(2*E*I);
    M=P*(L-x);
    Q=-P*ones(size(x));
    
    M=simplify(int(Q,x)-subs(int(Q,x),x,L))
    
    t=simplify(int(M/E/I,x)-subs(int(M/E/I,x),x,0))
    
    w=int(t,x)-subs(int(t,x),x,0)
    
  case 'cantilever_dist'
    
    w=P*x.^2.*(6*L^2-4*L*x+x.^2)/(24*E*I);
    t=P*x.*(3*L^2-3*L*x+x.^2)/(6*E*I);
    M=P*(L^2-2*L*x+x.^2);
    Q=-P*(L-x);
    
    M=simplify(int(Q,x)-subs(int(Q,x),x,L))
    
    t=simplify(int(M/E/I,x)-subs(int(M/E/I,x),x,0))
    
    w=int(t,x)-subs(int(t,x),x,0)
    
  otherwise
    
    warning('Case not defined')
    
end
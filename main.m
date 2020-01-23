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

%% main
clear

%% define the model
model.F=-1e4;
model.xF=0.25;
model.L=1;
model.B=0.1;
model.H=0.1;
model.tb=0.1;
model.th=0.01;
model.E=2.1e11;
model.nu=0.285;
model.rho=7.8e3;
model.A=model.B*model.H;
model.I=model.tb*(model.H-2*model.th)^3/12+model.B*(model.H^3-(model.H-2*model.th)^3)/12;
model.casename='simple_pt';

%% compute analytical solution from Euler-Bernoulli theory
beam=computeEulerBernoulli(model);

%% compute with comsol
comsol=runCOMSOLBeam(model,0);
wc=mpheval(comsol,'v','edim',1);
tc=mpheval(comsol,'dtang(v, x)','edim',1);
Mc=mpheval(comsol,'beam.Mzl','edim',1);
Qc=mpheval(comsol,'beam.Tyl','edim',1);

%% compare results
subplot(2,2,1)
plot(beam.x,beam.w,wc.p(1,:),wc.d1,'*')
xlabel('x'),ylabel('w')
subplot(2,2,2)
plot(beam.x,beam.t,tc.p(1,:),tc.d1,'*')
xlabel('x'),ylabel('t')
subplot(2,2,3)
plot(beam.x,beam.M,Mc.p(1,:),Mc.d1,'*')
xlabel('x'),ylabel('M')
subplot(2,2,4)
plot(beam.x,beam.Q,Qc.p(1,:),Qc.d1,'*')
xlabel('x'),ylabel('Q')

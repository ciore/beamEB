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
model.xF=1.5;
model.L=2;
model.B=[0.1 0.1 0.1];
model.H=[0.01 0.08 0.01];
model.El=[2.1e11 2.1e11 2.1e11];
model.nul=[0.285 0.285 0.285];
model.rhol=[7.8e3 7.8e3 7.8e3];
model.xsection='ibeam';
model.loadcase='simple_pt';
model.A=computeXArea(model);
model.I=computeSecondMomentArea(model);
model.E=computeEeq(model);
model.nu=computeNueq(model);
model.rho=computeRhoeq(model);

%% compute analytical solution from Euler-Bernoulli theory
beam=computeEulerBernoulli(model,1);

%% compute with comsol
comsol=runCOMSOLBeam(model);
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

%% 
function A=computeXArea(model)
  switch model.xsection
    case 'ibeam'
      A=sum(model.B.*model.H);
  end
end

%% 
function I=computeSecondMomentArea(model)
  switch model.xsection
    case 'ibeam'
      I0=model.B.*model.H.^3/12;
      A=model.B.*model.H;
      d=([0 cumsum(model.H(1:end-1))]+model.H/2-sum(model.H)/2);
      I=sum(I0+A.*d.^2);
  end
end

%% 
function Eeq=computeEeq(model)
  switch model.xsection
    case 'ibeam'
      I0=model.B.*model.H.^3/12;
      A=model.B.*model.H;
      d=([0 cumsum(model.H(1:end-1))]+model.H/2-sum(model.H)/2);
      I=(I0+A.*d.^2);
      Eeq=sum(model.El.*I)/sum(I);
  end
end

%% 
function nueq=computeNueq(model)
  switch model.xsection
    case 'ibeam'
      A=model.B.*model.H;
      nueq=sum(model.nul.*A)/sum(A);
  end
end

%% 
function rhoeq=computeRhoeq(model)
  switch model.xsection
    case 'ibeam'
      A=model.B.*model.H;
      rhoeq=sum(model.rhol.*A)/sum(A);
  end
end
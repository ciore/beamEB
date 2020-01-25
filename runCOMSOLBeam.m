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

function out = model(inputs,update)

import com.comsol.model.*
import com.comsol.model.util.*

if nargin<2
  update=0;
end

if update
  
  model = mphload('beam.mph');
  
  model.param.set('L', num2str(inputs.L,12), 'Length');
  model.param.set('A', num2str(inputs.A,12), 'X-Section Area');
  model.param.set('I', num2str(inputs.I,12), '2nd Moment Area');
  model.param.set('F', num2str(inputs.F,12), 'Force');
  model.param.set('xF', num2str(inputs.xF,12), 'x of Force');
  model.param.set('Young', num2str(inputs.E,12), 'Young''s mod');
  model.param.set('Poisson', num2str(inputs.nu,12), 'Poisson''s ratio');
  model.param.set('Rho', num2str(inputs.rho,12), 'density');
   
  model.sol('sol1').runAll;
  
%   model.sol('sol2').runAll;

  mphsave(model,'tmp')
  
else

  model = ModelUtil.create('Model');
  
  model.param.set('L', num2str(inputs.L,12), 'Length');
  model.param.set('A', num2str(inputs.A,12), 'X-Section Area');
  model.param.set('I', num2str(inputs.I,12), '2nd Moment Area');
  model.param.set('F', num2str(inputs.F,12), 'Force');
  model.param.set('xF', num2str(inputs.xF,12), 'x of Force');
  model.param.set('Young', num2str(inputs.E,12), 'Young''s mod');
  model.param.set('Poisson', num2str(inputs.nu,12), 'Poisson''s ratio');
  model.param.set('Rho', num2str(inputs.rho,12), 'density');
  
  model.component.create('comp1', true);

  model.component('comp1').geom.create('geom1', 2);
  model.component('comp1').geom('geom1').create('ls1', 'LineSegment');
  model.component('comp1').geom('geom1').feature('ls1').set('specify1', 'coord');
  model.component('comp1').geom('geom1').feature('ls1').set('specify2', 'coord');
  model.component('comp1').geom('geom1').feature('ls1').set('coord2', {'L' '0'});
  model.component('comp1').geom('geom1').run;

  model.component('comp1').physics.create('beam', 'HermitianBeam', 'geom1');
  model.component('comp1').physics('beam').field('displacement').field('u');
  model.component('comp1').physics('beam').field('displacement').component({'u' 'v' 'w'});
  switch inputs.loadcase
    case 'simple_pt'
      model.component('comp1').geom('geom1').create('pt1', 'Point');
      model.component('comp1').geom('geom1').feature('pt1').set('p', {'xF' '0'});
      model.component('comp1').geom('geom1').run;
      model.component('comp1').physics('beam').create('pdr1', 'DispRot0', 0);
      model.component('comp1').physics('beam').feature('pdr1').selection.set([1]);
      model.component('comp1').physics('beam').feature('pdr1').set('Direction', [1; 1; 0]);
      model.component('comp1').physics('beam').create('pdr2', 'DispRot0', 0);
      model.component('comp1').physics('beam').feature('pdr2').selection.set([3]);
      model.component('comp1').physics('beam').feature('pdr2').set('Direction', [0; 1; 0]);
      model.component('comp1').physics('beam').create('pl1', 'PointLoad', 0);
      model.component('comp1').physics('beam').feature('pl1').selection.set([2]);
      model.component('comp1').physics('beam').feature('pl1').set('Fp', {'0'; 'F'; '0'});
    case 'cantilever_ptend'
      model.component('comp1').physics('beam').create('fix1', 'Fixed', 0);
      model.component('comp1').physics('beam').feature('fix1').selection.set([1]);
      model.component('comp1').physics('beam').create('pl1', 'PointLoad', 0);
      model.component('comp1').physics('beam').feature('pl1').selection.set([2]);
      model.component('comp1').physics('beam').feature('pl1').set('Fp', {'0'; 'F'; '0'});
    case 'cantilever_dist'
      model.component('comp1').physics('beam').create('fix1', 'Fixed', 0);
      model.component('comp1').physics('beam').feature('fix1').selection.set([1]);
      model.component('comp1').physics('beam').create('el1', 'EdgeLoad', 1);
      model.component('comp1').physics('beam').feature('el1').selection.set([1]);
      model.component('comp1').physics('beam').feature('el1').set('FeperLength', {'0'; 'F'; '0'});
    otherwise 
      warning('Case not defined')
  end
  model.component('comp1').physics('beam').feature('emm1').set('E_mat', 'userdef');
  model.component('comp1').physics('beam').feature('emm1').set('E', 'Young');
  model.component('comp1').physics('beam').feature('emm1').set('nu_mat', 'userdef');
  model.component('comp1').physics('beam').feature('emm1').set('nu', 'Poisson');
  model.component('comp1').physics('beam').feature('emm1').set('rho_mat', 'userdef');
  model.component('comp1').physics('beam').feature('emm1').set('rho', 'Rho');
  model.component('comp1').physics('beam').feature('csd1').set('area', 'A');
  model.component('comp1').physics('beam').feature('csd1').set('Izz', 'I');

  model.component('comp1').mesh.create('mesh1');
  model.component('comp1').mesh('mesh1').create('auto_f1', 'Edge');
  model.component('comp1').mesh('mesh1').feature('auto_f1').selection.remaining;
  model.component('comp1').mesh('mesh1').feature('size').set('custom', 'on');
  model.component('comp1').mesh('mesh1').feature('size').set('hmax', 0.125);
  model.component('comp1').mesh('mesh1').run;

  model.study.create('std1');
  model.study('std1').create('stat', 'Stationary');
  model.sol.create('sol1');
  model.sol('sol1').study('std1');
  model.sol('sol1').attach('std1');
  model.sol('sol1').create('st1', 'StudyStep');
  model.sol('sol1').create('v1', 'Variables');
  model.sol('sol1').create('s1', 'Stationary');
  model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
  model.sol('sol1').runAll;
 
%   model.study.create('std2');
%   model.study('std2').create('eig', 'Eigenfrequency');
%   model.study('std2').feature('eig').set('neigs', 1);
%   model.sol.create('sol2');
%   model.sol('sol2').study('std2');
%   model.sol('sol2').attach('std2');
%   model.sol('sol2').create('st1', 'StudyStep');
%   model.sol('sol2').create('v1', 'Variables');
%   model.sol('sol2').create('e1', 'Eigenvalue');
%   model.sol('sol2').runAll;
  
  mphsave(model,'beam')
  
end

out = model;


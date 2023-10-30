% This script is written and read by pdetool and should NOT be edited.
% There are two recommended alternatives:
% 1) Export the required variables from pdetool and create a MATLAB script
%    to perform operations on these.
% 2) Define the problem completely using a MATLAB script. See
%    https://www.mathworks.com/help/pde/examples.html for examples
%    of this approach.
function pdemodel
[pde_fig,ax]=pdeinit;
pdetool('appl_cb',1);
set(ax,'DataAspectRatio',[1/100 0.97058823529411775/100 1/100]);
set(ax,'PlotBoxAspectRatio',[1.4999999999999998/100 1/100 1.7647058823529411/100]);
set(ax,'XLim',[-0.29999999999999999/100 1.3999999999999999/100]);
set(ax,'YLim',[-0.10000000000000001/100 1/100]);
set(ax,'XTickMode','auto');
set(ax,'YTickMode','auto');
pdetool('gridon','on');

% Geometry description:
pderect([0.0010000000000000001 0.011000000000000001 0.00050000000000000003 0.0055000000000000004],'R1');
pderect([0 0.0010000000000000001 0.0045000000000000007 0.0060000000000000009],'R2');
pderect([0.0040000000000000002 0.01 0.0044999999999999996 0.0014999999999999999],'R6');
pderect([0.01 0.0090000000000000002 0.0034999999999999998 0.0025],'R5');
pderect([0.0059999999999999998 0.0069999999999999996 0.0044999999999999996 0.0014999999999999999],'R7');
pdeellip(0.0010000000000000001,0.0029999999999999999,0.0010000000000000001,0.0025,...
0,'E1');
pdeellip(0.0040000000000000002,0.0029999999999999999,0.0010000000000000001,0.0025,...
0,'E2');
pderect([0 0.0010000000000000001 0 0.0015000000000000002],'R3');
pderect([0.0010000000000000001 0.0040000000000000002 0.0044999999999999996 0.0014999999999999999],'R4');
set(findobj(get(pde_fig,'Children'),'Tag','PDEEval'),'String','((((R1+R3+R2+E1)-R4)+E2)-R6)+R7+R5')
% Boundary conditions:
pdetool('changemode',0)
pdetool('removeb',[45 44 50 49 48 47 43 42 22 20 10 ]);
pdesetbd(40,...
'dir',...
1,...
'1',...
'0')
pdesetbd(39,...
'dir',...
1,...
'1',...
'0')
pdesetbd(38,...
'dir',...
1,...
'1',...
'0')
pdesetbd(37,...
'dir',...
1,...
'1',...
'0')
pdesetbd(35,...
'dir',...
1,...
'1',...
'0')
pdesetbd(34,...
'dir',...
1,...
'1',...
'0')
pdesetbd(32,...
'dir',...
1,...
'1',...
'0')
pdesetbd(31,...
'dir',...
1,...
'1',...
'0')
pdesetbd(30,...
'dir',...
1,...
'1',...
'0')
pdesetbd(28,...
'dir',...
1,...
'1',...
'0')
pdesetbd(27,...
'dir',...
1,...
'1',...
'0')
pdesetbd(25,...
'dir',...
1,...
'1',...
'0')
pdesetbd(24,...
'dir',...
1,...
'1',...
'0')
pdesetbd(23,...
'dir',...
1,...
'1',...
'0')
pdesetbd(22,...
'dir',...
1,...
'1',...
'0')
pdesetbd(21,...
'dir',...
1,...
'1',...
'0')
pdesetbd(20,...
'dir',...
1,...
'1',...
'0')
pdesetbd(19,...
'dir',...
1,...
'1',...
'0')
pdesetbd(18,...
'dir',...
1,...
'1',...
'0')
pdesetbd(17,...
'dir',...
1,...
'1',...
'0')
pdesetbd(15,...
'dir',...
1,...
'1',...
'0')
pdesetbd(14,...
'dir',...
1,...
'1',...
'0')
pdesetbd(12,...
'dir',...
1,...
'1',...
'0')
pdesetbd(11,...
'dir',...
1,...
'1',...
'0')
pdesetbd(10,...
'dir',...
1,...
'1',...
'0')
pdesetbd(9,...
'dir',...
1,...
'1',...
'0')
pdesetbd(8,...
'dir',...
1,...
'1',...
'0')
pdesetbd(7,...
'dir',...
1,...
'1',...
'0')
pdesetbd(6,...
'dir',...
1,...
'1',...
'0')
pdesetbd(5,...
'dir',...
1,...
'1',...
'0')
pdesetbd(4,...
'dir',...
1,...
'1',...
'0')
pdesetbd(3,...
'dir',...
1,...
'1',...
'0')
pdesetbd(2,...
'dir',...
1,...
'1',...
'0')
pdesetbd(1,...
'dir',...
1,...
'1',...
'0')

% Mesh generation:
setappdata(pde_fig,'Hgrad',1.3);
setappdata(pde_fig,'refinemethod','regular');
setappdata(pde_fig,'jiggle',char('on','mean',''));
setappdata(pde_fig,'MesherVersion','preR2013a');
pdetool('initmesh')
pdetool('refine')

% PDE coefficients:
pdeseteq(1,...
'1.0',...
'0.0',...
'10.0',...
'1.0',...
'0:10',...
'0.0',...
'0.0',...
'[0 100]')
setappdata(pde_fig,'currparam',...
['1.0 ';...
'0.0 ';...
'10.0';...
'1.0 '])

% Solve parameters:
setappdata(pde_fig,'solveparam',...
char('0','1584','10','pdeadworst',...
'0.5','longest','0','1E-4','','fixed','Inf'))

% Plotflags and user data strings:
setappdata(pde_fig,'plotflags',[1 1 1 1 1 1 1 1 0 0 0 1 1 0 0 0 0 1]);
setappdata(pde_fig,'colstring','');
setappdata(pde_fig,'arrowstring','');
setappdata(pde_fig,'deformstring','');
setappdata(pde_fig,'heightstring','');
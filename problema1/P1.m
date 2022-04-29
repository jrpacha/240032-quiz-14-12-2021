clearvars
close all

% ==================================== P1 ==================================
fileSols = 'solutionP1.txt';
fOut = fileSols;

%Real constants: Materials and sections area
Area=3250.0;     %section Area (in mm^2);
Y=2.0e5;         %Young modulus, in N/mm^2 (1GP=1.0 kN/mm^2)

%F=[29705;-20515;0.0]; %force applied to node 6
F=[300000;-200000;0.0];

% Data
fOut=fopen(fileSols,'w');
%fprintf(fOut,'Length.............................. L = %e (mm)\n',      D);
%fprintf(fOut,'Force............................... F = %e (N)\n',       F);
fprintf(fOut,'Area............................. Area = %e (mm^2)\n', Area);
fprintf(fOut,'Young modulus....................... Y = %e (N/mm^2)\n',  Y);
%fprintf(fOut,'Num. nod. displ............... nodDisp = %d\n',     nodDisp);
fprintf(fOut,'\n');

%Goemetry
H=1800.0*sqrt(3.0);      %mm
z=3536.0;                %mm

% nodes=[0.0,0.0,0.0;
%        3600.0,0.0,0.0;
%        7200.0,0.0,0.0;
%        10800.0,0.0,0.0;
%        1800.0,H,0.0;
%        5400.0,H,z;
%        9000.0,H,0.0];

nodes = [
    0 0 0;
    3600 0 0;
    7200 0 0;
    10800 0 0;
    1800 3118 0;
    5400 3118 z;
    9000 3118 0;
    ];

elem = [ 
    1 2;
    2 3;
    3 4;
    1 5;
    5 2;
    5 6;
    2 6;
    6 3;
    6 7;
    3 7;
    7 4;
    ];
      
numNod=size(nodes,1);
numElem=size(elem,1);
ndim=size(nodes,2);

%Real constants
A=Area*ones(1,numElem);
E=Y*ones(1,numElem);

%Assembly
u=zeros(ndim*numNod,1);
Q=zeros(ndim*numNod,1);
K=zeros(ndim*numNod);

%%% Part (A)
clc
fprintf('\t\t ** Problema 1 **\n')
Ke=spatialLinkStiffMatrix(nodes,elem,7,E,A);
fprintf('(a) The value of the entry (3,4) of the local stiff matrix\n')
fprintf('    for the seventh element is K7(3,4) = %.5e\n',Ke(3,4))
fprintf('    Hint. The value of K7(4,5) = %.5e\n',Ke(4,5))     

for e=1:numElem
    Ke=spatialLinkStiffMatrix(nodes,elem,e,E,A);
    rows=[ndim*elem(e,1)-2,ndim*elem(e,1)-1,ndim*elem(e,1),...
          ndim*elem(e,2)-2,ndim*elem(e,2)-1,ndim*elem(e,2)];
    cols=rows;
    K(rows,cols)=K(rows,cols)+Ke; %Assembly
end

%Loads
%node 6:
nod=6;
Q(ndim*nod-2)=F(1);     %N
Q(ndim*nod-1)=F(2);     %N
Q(ndim*nod)=F(3);       %N 

%Boundary Conditions
fixedNods=[];

%node 1:
nod=1;
fixedNods=[fixedNods,ndim*nod-2];     %(u1_x=0);
fixedNods=[fixedNods,ndim*nod-1];     %(u1_y=0);
fixedNods=[fixedNods,ndim*nod];       %(u1_z=0);
u(ndim*nod-2)=0.0;
u(ndim*nod-1)=0.0;
u(ndim*nod)=0.0;

%node 2:
nod=2;
fixedNods=[fixedNods,ndim*nod-2];     %(u2_x=0);
fixedNods=[fixedNods,ndim*nod-1];     %(u2_y=0);
fixedNods=[fixedNods,ndim*nod];       %(u2_z=0);
u(ndim*nod-2)=0.0;
u(ndim*nod-1)=0.0;
u(ndim*nod)=0.0;

%node 3:
nod=3;
fixedNods=[fixedNods,ndim*nod-2];     %(u3_x=0);
fixedNods=[fixedNods,ndim*nod-1];     %(u3_y=0);
fixedNods=[fixedNods,ndim*nod];       %(u3_z=0);
u(ndim*nod-2)=0.0;
u(ndim*nod-1)=0.0;
u(ndim*nod)=0.0;

%node 4:
nod=4;
fixedNods=[fixedNods,ndim*nod-2];     %(u3_x=0);
fixedNods=[fixedNods,ndim*nod-1];     %(u3_y=0);
fixedNods=[fixedNods,ndim*nod];       %(u3_z=0);
u(ndim*nod-2)=0.0;
u(ndim*nod-1)=0.0;
u(ndim*nod)=0.0;

%node 5:
nod=5;
fixedNods=[fixedNods,ndim*nod];       %(u3_z=0);
u(ndim*nod)=0.0;

%node 7:
nod=7;
fixedNods=[fixedNods,ndim*nod];       %(u3_z=0);
u(ndim*nod)=0.0;

%Reduced system
freeNods=setdiff(1:ndim*numNod,fixedNods);
Qm=Q(freeNods,1)-K(freeNods,fixedNods)*u(fixedNods);
Km=K(freeNods,freeNods);

%Solve the reduced system
um=Km\Qm;
u(freeNods)=um;

%%% Part (b)
fprintf('(b) The maximum of the displacements of the z-component of all\n')
fprintf('    nodes is: %.5e\n',max(u(3:3:end)))
fprintf('    Hint. The maximum of the x-displacements is %.5e\n',...
    max(u(1:3:end)))

%%% Part (c)
longElem = zeros(numElem,1);
U = [u(1:3:end),u(2:3:end),u(3:3:end)];
for e = 1:numElem
    r1 = nodes(elem(e,1),:) + U(elem(e,1),:);
    r2 = nodes(elem(e,2),:) + U(elem(e,2),:);
    longElem(e) = norm(r1-r2);
end
[maxLongElem,e] = max(longElem);
rows = [ndim*elem(e,1)-2;ndim*elem(e,1)-1;ndim*elem(e,1)];
r2 = nodes(elem(e,1),:)' + u(rows);

fprintf('(c) The final length of the most deformed bar is: %.5e\n',...
    maxLongElem)
fprintf('    Hint. The y-component of the 1st. vertex of the bar\n')
fprintf('          of maximum deformation is %.5e\n',r2(2))


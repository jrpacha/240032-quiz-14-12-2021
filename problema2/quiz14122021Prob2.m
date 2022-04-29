clearvars
close all

fileName='solution.xlsx';

kc= 1;

temp = 10.0;
N = 1240;
tempN = 150.0;

elemE = 1111;
nodeHintA = 555;

thresholdTempB = 48.0;
thresholdTempHintB = 70.0;

%eval('mesh2x2Quad');
eval('wirecoathang');

numNod=size(nodes,1);
numElem=size(elem,1);
%numbering=1;
numbering=0;
%figure()
plotElementsOld(nodes,elem,numbering);
[indNodBd, indElemBd, indLocalEdgBd, edges] = boundaryNodes(nodes, elem);
hold on
plot(nodes(indNodBd,1),nodes(indNodBd,2),'o','markerFaceColor','red',...
    'markerSize',6)
plot(nodes(N,1),nodes(N,2),'o','markerFaceColor','green',...
    'markerSize',6)
hold off

% a11=1;
% a12=1;
% a21=a12;
% a22=a11;
% a00=1;
% f=1;

a11=kc;    %Coefficients for Exercise 2
a12=0.0;
a21=a12;
a22=a11;
a00=0.0;
f=0.0;

coeff=[a11,a12,a21,a22,a00,f];
K=zeros(numNod);
F=zeros(numNod,1);
for e=1:numElem
    [Ke,Fe]=bilinearQuadElement(coeff,nodes,elem,e);
    %
    % Assemble the elements
    %
    rows=[elem(e,1); elem(e,2); elem(e,3); elem(e,4)];
    colums= rows;
    K(rows,colums)=K(rows,colums)+Ke; %assembly
    if (coeff(6) ~= 0)
        F(rows)=F(rows)+Fe;
    end
end 

%Boundary Conditions (BC)
%Impose Boundary Conditions for this example. In this case only essential 
%and natural BC are considered. We will do a thermal example for convection 
%BC (mixed).
%indRigth=[3;6;9];
%indLeft=[1;4;7];

fixedNodes=[indNodBd',N]; %Fixed Nodes (global num.)
freeNodes = setdiff(1:numNod,fixedNodes); %Complementary of fixed nodes

% ------------ Essential BC
u=zeros(numNod,1); %initialize u vector
u(indNodBd)=temp; %all of them are zero
u(N)=tempN;
Fm=F(freeNodes)-K(freeNodes,fixedNodes)*u(fixedNodes);%here u can be 
                                                      %different from zero 
                                                      %only for fixed nodes
% ------------ Natural BC
Q=zeros(numNod,1); %initialize Q vector
Q(freeNodes)=0; %all of them are zero
% modify the linear system, this is only valid if BC = 0.
Km=K(freeNodes,freeNodes);
Fm=Fm+Q(freeNodes);

%Compute the solution
%solve the System

format short e; %just to a better view of the numbers
um=Km\Fm;
u(freeNodes)=um;

%PostProcess: Compute secondary variables and plot results
Q=K*u-F;

%%% Part A
clc
fprintf('\t\t ** Problema 2 **\n')
nods = elem(elemE,:);
meanValue = sum(u(nods))/size(elem,2);
fprintf('(a) Mean value of u for the element e = %d, <u> = %.5e\n',...
    elemE,meanValue)
fprintf('    Hint. The solution for node N = %d, u(%d) = %.5e\n',...
    nodeHintA,nodeHintA,u(nodeHintA))
    
%%% Part B
nods = find(u > thresholdTempB);
fprintf('(b) Num. of nodes with u strictly greater than %.2f: %d\n',...
    thresholdTempB, length(nods))
nods = find(u < thresholdTempHintB);
fprintf('    Hint. There are %d nodes below %.2f\n',...
    length(nods),thresholdTempHintB)

%%% Part C
nodsTop = find(nodes(:,2) > 9.99);
figure()
plotElementsOld(nodes,elem,numbering);
hold on 
plot(nodes(nodsTop,1),nodes(nodsTop,2),'o','markerFaceColor','red',...
    'markerSize',6)
hold off
fprintf('(c) Sum of the flows Q(i) over all of the nodes in the top\n')
fprintf('    side of the boundary, sum(Q(i),i=1...nodsTop) = %.5e\n',...
    sum(Q(nodsTop)))
fprintf('    Hint: the minimal value of Q among the nodes in the top\n')
fprintf('          side of the boundary is %.5e\n',min(Q(nodsTop)))


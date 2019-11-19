% [globalSystem,boundStruct,meshStruct]=INPUTDATA(meshStruct)
% This file defines the element properties, prescribed displacements 
% (essential boundary conditions) and applied loads for the TRUSS2D code. 

% The portions of the code that might change for each new
% problem are clearly indicated. 
function [globalSystem,boundStruct,meshStruct]=InputData(meshStruct) 

% unpack necessary input
numEls=meshStruct.numEls;
numDOF=meshStruct.numDOF;
numNodes=meshStruct.numNodes;

%% Element properties (DEFINE THIS FOR EACH PROBLEM)
% These properties may be the same or different for each element in the
% mesh (be sure to use consistent units). 
elArea 	= 1e-4*ones(numEls,1);   	% Elements area  
elYM    = 2e11*ones(numEls,1);   	% Young's Modulus 
 
%% Applied forces (DEFINE THIS FOR EACH PROBLEM)
% Each row is the global node number, the DOF, and the value for any applied loads
appF=1000;
% for example, appForces=[3 2 20e3]; means that global node number 3 has an
% applied load in the y direction with magnitude 20e3
% appForces=[2 1 appF]; 
% appForces=[2 2 -appF];  
appForces = [3 3 -appF];
%% Prescribed displacement boundary conditions (DEFINE THIS FOR EACH PROBLEM)
% Each row is the global node number, the DOF, and the value for any
% essential BCs. for example, essBCs=[3 2 0;] means that global node number 3 has a 
% required displacement of 0 in the y direction
% essBCs=[1 1 0; 1 2 0; 2 2 0];
% essBCs=[1 1 0; 1 2 0; 3 1 0; 3 2 0];
essBCs=[1 1 0; 1 2 0; 1 3 0; 3 1 0; 3 2 0; 3 3 0; 4 1 0; 4 2 0; 4 3 0];

% initialize global system of equations
numEq=numNodes*numDOF;
F=zeros(numEq,1);
k_hat=zeros(numEq,1);
K=zeros(numEq);

% initialize initial displacement vector
% NOTE: in order for the Newton-Raphson method to iterate to the (one)
% physically-sensible solution, the user needs to determine if this should
% be a negative or positive initialization
d= -1*ones(numEq,1);
% needed 
% initialize Newton-Raphson solver parameters
iter_max = 50;
tol = 1E-8;

% Map the applied loads to the proper location in the global force vector
for frc=1:size(appForces,1)
    gnn=appForces(frc,1);        % global node number for this applied load
    gdof=(gnn-1)*numDOF+appForces(frc,2); % global DOF for the applied load
    val=appForces(frc,3);        % value of the applied load
    F(gdof)=val; % populate the global force vector with applied loads
end

% Package variables into the output structs
globalSystem.k_hat=k_hat;
globalSystem.K=K;
globalSystem.F=F;
globalSystem.d=d;
globalSystem.iter_max=iter_max;
globalSystem.tol=tol;

boundStruct.essBCs   =essBCs;
boundStruct.appForces=appForces;

meshStruct.numEq =numEq;
meshStruct.elArea=elArea;
meshStruct.elYM  =elYM;


 
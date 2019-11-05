function [globalSystem,boundStruct,meshStruct]=InputData(meshStruct) 
% [globalSystem,boundStruct,meshStruct]=INPUTDATA(meshStruct)
% This file defines the element properties, prescribed
% displacements (essential boundary conditions) and applied loads for the
% TRUSS2D code. 
% The portions of the code that might change for each new
% problem are clearly indicated. 

% unpack necessary input
numEls=meshStruct.numEls;
numDOF=meshStruct.numDOF;
numNodes=meshStruct.numNodes;

% Element properties (be sure to use consistent units). These properties
% may be the same or different for each element in the mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE THIS FOR EACH PROBLEM
elArea 	= 1e-4*ones(numEls,1);   	% Elements area  
elYM    = 2e11*ones(numEls,1);   	% Young's Modulus 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Applied forces. Each row is the global node number, the DOF, and
% the value for any applied loads
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE THIS FOR EACH PROBLEM
myf=1000;
% for example, appForces=[3 2 20e3]; means that global node number 3 has an
% applied load in the y direction with magnitude 20e3
appForces=[2 2 -myf];  
% appForces=[2 1 myf];     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prescribed displacement boundary conditions. Each row is the global node
% number, the DOF, and the value for any essential BCs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE THIS FOR EACH PROBLEM
% for example, essBCs=[3 2 0;] means that global node number 3 has a 
% required displacement of 0 in the y direction
essBCs=[1 1 0; 1 2 0; 3 1 0; 3 2 0];
% essBCs=[1 1 0; 1 2 0; 1 2 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize global system of equations
numEq=numNodes*numDOF;
F=zeros(numEq,1);
k_hat=zeros(numEq,1);
K=zeros(numEq);

% initialize initial displacement vector 
d=rand(numEq,1);
% initialize Newton-Raphson solver parameters
iter_max = 50;
tol = 1E-12;

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


 
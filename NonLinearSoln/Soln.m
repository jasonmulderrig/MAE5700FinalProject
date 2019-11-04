function globalSystem = Soln(globalSystem,meshStruct,boundStruct)
% globalSystem = SOLN(globalSystem,meshStruct,boundStruct)
% Apply the essential boundary conditions and solve the global system for
% the nodal displacements for TRUSS2D. 
% last edit: 10 July 2015 H. Ritz

% unpack necessary input
K=globalSystem.K;
k_hat=globalSystem.k_hat;
F=globalSystem.F;
d=globalSystem.d;
essBCs=boundStruct.essBCs;
numDOF=meshStruct.numDOF;
numEq =meshStruct.numEq;

% Initial guess for {d} will be the same as initialization 
% First calculation of residual
res = K*d*k_hat'*d - F;
% Iteration for Newton Raphson Method


% % partition the matrix K, vectors f and d first we want to find the global
% % DOFs with essential boundary conditions
% numEBC=size(essBCs,1); % number of essential boundary conditions
% % initialize array of global indices to essential boundary conditions
% essDOF=zeros(numEBC,1);
% for ebc=1:numEBC
%     % essDOF stores the index to the degrees of freedom with essential
%     % boundary conditions
%     essDOF(ebc)=(essBCs(ebc,1)-1)*numDOF+essBCs(ebc,2); 
% end
% indF=setdiff(1:numEq,essDOF); % this returns the indices to the DOF that 
%                               % do NOT have essential boundary conditions
%                               % (free DOF)
% K_E	= K(essDOF,essDOF);       % Extract K_E matrix 
% K_F	= K(indF,indF);           % Extract K_F matrix 
% K_EF    = K(essDOF,indF);     % Extract K_EF matrix
% f_F  	= F(indF);            % Extract f_F vector
% d_E  	= essBCs(:,3);        % Extract d_E vector
% f_E     = F(essDOF);          % Extract f_E vector 
%     
% % solve for d_F
% d_F	=K_F\( f_F - K_EF'* d_E);
%  
% % reconstruct the global solution d
% d(essDOF)=d_E;                
% d(indF)=d_F;
% 
% % compute the reaction forces on the DOF with essential BCs
% reactionVec = K_E*d_E+K_EF*d_F-f_E;
% 
% % Package variables into the output structs
% globalSystem.d=d;
% globalSystem.reactionVec=reactionVec;


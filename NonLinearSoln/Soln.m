% globalSystem = SOLN(globalSystem,meshStruct,boundStruct)
% Apply the essential boundary conditions and solve the global system for
% the nodal displacements for TRUSS2D. 
function globalSystem = Soln(globalSystem,meshStruct,boundStruct)

% unpack necessary input
K=globalSystem.K;
k_hat=globalSystem.k_hat;
F=globalSystem.F;
d=globalSystem.d;
iter_max=globalSystem.iter_max;
tol=globalSystem.tol;
essBCs=boundStruct.essBCs;
numDOF=meshStruct.numDOF;
numEq =meshStruct.numEq;

% find the global DOFs with essential boundary conditions
numEBC=size(essBCs,1); % number of essential boundary conditions
% initialize array of global indices to essential boundary conditions
indE=zeros(numEBC,1);
for ebc=1:numEBC
    % indE stores the index to the degrees of freedom with essential
    % boundary conditions
    indE(ebc)=(essBCs(ebc,1)-1)*numDOF+essBCs(ebc,2); 
end
% this returns the indices to the DOF that do NOT have essential boundary 
% conditions (free DOF)
indF=setdiff(1:numEq,indE); 
% impose the essential boundary conditions on the displacement solution 
% vector for each Newton-Raphson iteration
d(indE) = essBCs(:,3);

% Use the Newton-Raphson method to determine the displacement solution
% vector to the nonlinear system. Note that the initial displacement
% solution vector will simply be the initialized zero vector.
iter = 1;
% Array to store residual values for output in PresentResults.m
storeRes = []; 
while iter < iter_max
    % calculate the residual
    R = K*d*k_hat'*d - F;
    % partition the residual
    R_F = R(indF);
    % Store residual values 
    storeRes(iter) = max(abs(R_F));
    if max(abs(R_F)) <= tol % Newton-Raphson method converged on an individual DOF basis
        break;
    else % Newton-Raphson method did not converge on an individual DOF basis
        % calculate the Jacobian of the residual
        J = K*(k_hat'*d) + K*(d*k_hat');
        % partition the J matrix
        J_F	 = J(indF,indF); % Extract J_F matrix
        if (J_F == zeros(size(J_F,1),size(J_F,2)))
            error('Define a new initial guess for the displacements so that the Jacobian is not zero');
        end
        % invert the partitioned J matrix
        J_inv_F	 = inv(J_F);
        % update the displacement solution vector
        d(indF) = d(indF) - J_inv_F*R_F;
    end
    iter = iter+1;
end   

F_int = K*d*k_hat'*d;
f_int = F_int(indE);
f_ext = F(indE);

reactionVec = f_int - f_ext;

% Package variables into the output structs
globalSystem.d=d;
globalSystem.reactionVec=reactionVec;
globalSystem.iter=iter;
globalSystem.storeRes = storeRes; 
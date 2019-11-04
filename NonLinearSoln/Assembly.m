function globalSystem=Assembly(globalSystem,meshStruct)
% globalSystem=ASSEMBLY(globalSystem,meshStruct)
% Assemble global stiffness matrix K for the TRUSS2D3D code. 
% last edit: 30 July 2015 H. Ritz

% unpack necessary inputs
nnpe=meshStruct.nnpe;
numDOF=meshStruct.numDOF;
numEls=meshStruct.numEls;
gatherMat=meshStruct.gatherMat;
K=globalSystem.K;
k_hat=globalSystem.k_hat;

% for each element, make the local stiffness matrix and then assemble
% "ke" into our global stiffness matrix K
for e=1:numEls % for each element
    [K_elem, k_hat_elem] = TrussElem(e,meshStruct); % make the local stiffness matrix
    for Lrw = 1 : (nnpe*numDOF)
        Grw = gatherMat(e,Lrw); % global row index
        % Assemble local k_hat vector into the global k_hat vector here
        k_hat(Grw) = k_hat(Grw) + k_hat_elem(Lrw);
        for Lcl = 1 : (nnpe*numDOF)
            Gcl = gatherMat(e,Lcl); % global column index
            % Assemble local stiffness matrix into the global
            % stiffness matrix here
            K(Grw,Gcl) = K(Grw,Gcl) + K_elem(Lrw,Lcl);
        end
    end
end

% Package variables into the output structs
globalSystem.K=K;
globalSystem.k_hat = k_hat; 

function [K_elem, k_hat_elem] = TrussElem(elmID,meshStruct)
% localstiffnessmatrix = TrussElem(elementnumber,meshStruct)
% generate the local stiffness matrix for use with TRUSS2D3D code.

% unpack necessary input
elCon  =meshStruct.elCon;
nCoords=meshStruct.nCoords;
elArea =meshStruct.elArea;
elYM   =meshStruct.elYM;
numDim =meshStruct.numDim;
nnpe   =meshStruct.nnpe;

gn1 = elCon(elmID,1); % Extract the global node numbers 
gn2 = elCon(elmID,2); % for the current element

x1=nCoords(gn1,1); y1=nCoords(gn1,2); % extract the nodal coordinates for
x2=nCoords(gn2,1); y2=nCoords(gn2,2); % each node in the current element
          
switch numDim % the local stiffness matrix has a different form
              % depending on the number of spatial dimensions
    case 2 % 2D problems        
        L=sqrt((x2-x1)^2+(y2-y1)^2); % L = length of the element
        a=(x2-x1)/L; % cosine of the angle of the element with the X axis
        b=(y2-y1)/L; % cosine of the angle of the element with the Y axis
        
        const = elArea(elmID)*elYM(elmID)/(L^2); % constant coefficient for the element
        % Rotation matrix for 2D elements
        Re = [a b 0 0;
              0 0 a b];
        % stiffness matrix for 2D truss element in the global coordinate system
        K_elem = const*Re'*[1 -1; -1 1]*Re;
        k_hat_elem = Re'*[-1; 1];
    case 3 % 3D problems
        z1=nCoords(gn1,3);      
        z2=nCoords(gn2,3); 
        L=sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2); % L = length of the element
        a=(x2-x1)/L; % cosine of the angle of the element with the X axis
        b=(y2-y1)/L; % cosine of the angle of the element with the Y axis
        c=(z2-z1)/L; % cosine of the angle of the element with the Z axis
        % constant coefficient for the element
        const = elArea(elmID)*elYM(elmID)/(L^2); 
        % Stiffness matrix for 3D truss element in the global coordinate
        % system
        Re = [a b c 0 0 0; 
              0 0 0 a b c]; 
        K_elem = const*Re'*[1 -1; -1 1]*Re;
        k_hat_elem = Re'*[-1; 1];
end

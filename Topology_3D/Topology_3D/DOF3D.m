function edofMat = DOF3D(nx,ny,nz)
ne = nx*ny*nz; 
nodegrd = reshape ( 1:(nx+1)*(ny+1),ny+1,nx+1);
nodeids = reshape ( nodegrd ( 1:end-1 ,1:end-1),nx*ny,1);
nodeidz = 0 : (ny +1)*( nx +1) : ( nz-1)*( ny+1)*( nx +1);
nodeids = repmat ( nodeids , size ( nodeidz ) )+repmat ( nodeidz , size ( nodeids ) );
edofVec = 3* nodeids (:)+1;
edofMat = repmat ( edofVec,1,24) + repmat ( [ 0 1 2 3*ny + [ 3 4 5 0 1 2 ] -3 -2 -1 3*( ny +1)*( nx +1)+[0 1 2 3*ny+[ 3 4 5 0 1 2 ] -3 -2 -1]],ne,1 );
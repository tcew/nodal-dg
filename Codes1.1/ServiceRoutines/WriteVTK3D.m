function WriteVTK3D(filename, Nout, vnames, varargin)
  nfields = nargin - 3;

  Globals3D;

  Npout = (Nout+1)*(Nout+2)*(Nout+3)/6;
  [rout,sout,tout] = EquiNodes3D(Nout);


  % build matrix to interpolate field data to equally spaced nodes
  Vout = Vandermonde3D(N, rout, sout, tout);
  interp = Vout*invV;


  % form symmetric tetrahedralization of local nodes
  startrow = zeros(Nout+1, Nout+1);

  sk = 1;
  for i=0:Nout
    for j=0:Nout-i
      startrow(j+1, i+1) = sk;
      sk = sk + Nout+1-i-j;
    end
  end

  % contruct tetrahedralization
  tet = zeros(Nout*Nout*Nout,4);

  sk = 1;
  for i=0:Nout
    for j=0:Nout-i
      for k=0:Nout-i-j-1
        % Add Tet 1
        tet(sk,1) = startrow(j+1, i+1)+k;
        tet(sk,2) = startrow(j+1, i+1)+k+1;
        tet(sk,3) = startrow(j+2, i+1)+k;
        tet(sk,4) = startrow(j+1, i+2)+k;
        sk = sk+1;

        if(k < Nout-i-j-1)
          % Add Tet 2
          tet(sk,1) = startrow(j+1, i+1)+k+1;
          tet(sk,2) = startrow(j+2, i+1)+k;
          tet(sk,3) = startrow(j+1, i+2)+k;
          tet(sk,4) = startrow(j+1, i+2)+k+1;
          sk = sk+1;

          % Add Tet 3
          tet(sk,1) = startrow(j+1, i+1)+k+1;
          tet(sk,2) = startrow(j+2, i+1)+k+1;
          tet(sk,3) = startrow(j+2, i+1)+k;
          tet(sk,4) = startrow(j+1, i+2)+k+1;
          sk = sk+1;

          % Add Tet 4
          tet(sk,1) = startrow(j+1, i+2)+k;
          tet(sk,2) = startrow(j+2, i+2)+k;
          tet(sk,3) = startrow(j+1, i+2)+k+1;
          tet(sk,4) = startrow(j+2, i+1)+k;
          sk = sk+1;

          % Add Tet 5
          tet(sk,1) = startrow(j+2, i+1)+k;
          tet(sk,2) = startrow(j+2, i+1)+k+1;
          tet(sk,3) = startrow(j+2, i+2)+k;
          tet(sk,4) = startrow(j+1, i+2)+k+1;
          sk = sk+1;
        end

        if(k < Nout-i-j-2)
          % Add Tet 6
          tet(sk,1) = startrow(j+1, i+2)+k+1;
          tet(sk,2) = startrow(j+2, i+1)+k+1;
          tet(sk,3) = startrow(j+2, i+2)+k;
          tet(sk,4) = startrow(j+2, i+2)+k+1;
          sk = sk+1;
        end
      end
    end
  end

  Ntet = sk-1;

  % create global tet mesh
  TET = [];
  for k=1:K
    TET = [TET; tet+(k-1)*Npout];
  end

  % interpolate node coordinates and field to equally spaced nodes
  xout = interp*x;
  yout = interp*y;
  zout = interp*z;

  uout = cell(1, nfields);
  for n=1:nfields
    uout{n} = interp*varargin{n};
  end

  Ntotal = length(xout(:));
  [nTET, nVERT] = size(TET);

  fid = fopen(filename, 'w');
  fprintf(fid, '# vtk DataFile Version 2');
  fprintf(fid, '\nNUDG simulation');
  fprintf(fid, '\nASCII');
  fprintf(fid, '\nDATASET UNSTRUCTURED_GRID\n');
  fprintf(fid, '\nPOINTS %d double', Ntotal);
  fprintf(fid, '\n%25.16e %25.16e %25.16e', [xout(:) yout(:) zout(:)]');
  fprintf(fid, '\nCELLS %d %d', nTET, nTET*5);
  fprintf(fid, '\n4 %10d %10d %10d %10d', (TET-1)');
  fprintf(fid, '\nCELL_TYPES %d', nTET);
  fprintf(fid, '\n%d', repmat([10], nTET, 1));
  fprintf(fid, '\nPOINT_DATA %d', Ntotal);
  for n=1:nfields
    fprintf(fid, '\nSCALARS %s double 1', vnames{n});
    fprintf(fid, '\nLOOKUP_TABLE default');
    fprintf(fid, '\n%25.16e', uout{n});
  end
  fclose(fid);
end

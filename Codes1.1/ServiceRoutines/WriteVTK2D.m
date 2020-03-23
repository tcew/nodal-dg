function WriteVTK2D(filename, Nout, vnames, varargin)
  nfields = nargin - 3;

  Globals2D;

  % build equally spaced grid on reference triangle
  Npout = (Nout+1)*(Nout+2)/2;
  rout = zeros(Npout,1); sout = zeros(Npout,1);
  sk = 1;
  for n=1:Nout+1
    for m=1:Nout+2-n
      rout(sk) = -1 + 2*(m-1)/Nout;
      sout(sk) = -1 + 2*(n-1)/Nout;
      counter(n,m) = sk; sk = sk+1;
    end
  end

  % build matrix to interpolate field data to equally spaced nodes
  interp = InterpMatrix2D(rout, sout);

  % build triangulation of equally spaced nodes on reference triangle
  tri = [];
  for n=1:Nout+1
    for m=1:Nout+1-n,
      v1 = counter(n,m);   v2 = counter(n,m+1);
      v3 = counter(n+1,m); v4 = counter(n+1,m+1);
      if(v4)
        tri = [tri;[[v1 v2 v3];[v2 v4 v3]]];
      else
        tri = [tri;[[v1 v2 v3]]];
      end
    end
  end

  % build triangulation for all equally spaced nodes on all elements
  TRI = [];
  for k=1:K
    TRI = [TRI; tri+(k-1)*Npout];
  end

  % interpolate node coordinates and field to equally spaced nodes
  xout = interp*x;
  yout = interp*y;
  zout = zeros(size(xout));

  uout = cell(1, nfields);
  for n=1:nfields
    uout{n} = interp*varargin{n};
  end

  Ntotal = length(xout(:));
  [nTRI, nVERT] = size(TRI);

  fid = fopen(filename, 'w');
  fprintf(fid, '# vtk DataFile Version 2');
  fprintf(fid, '\nNUDG simulation');
  fprintf(fid, '\nASCII');
  fprintf(fid, '\nDATASET UNSTRUCTURED_GRID\n');
  fprintf(fid, '\nPOINTS %d double', Ntotal);
  fprintf(fid, '\n%25.16e %25.16e %25.16e', [xout(:) yout(:) zout(:)]');
  fprintf(fid, '\nCELLS %d %d', nTRI, nTRI*4);
  fprintf(fid, '\n3 %10d %10d %10d', (TRI-1)');
  fprintf(fid, '\nCELL_TYPES %d', nTRI);
  fprintf(fid, '\n%d', repmat([5], nTRI, 1));
  fprintf(fid, '\nPOINT_DATA %d', Ntotal);
  for n=1:nfields
    fprintf(fid, '\nSCALARS %s double 1', vnames{n});
    fprintf(fid, '\nLOOKUP_TABLE default');
    fprintf(fid, '\n%25.16e', uout{n});
  end
  fclose(fid);
end

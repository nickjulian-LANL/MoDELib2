clc
clear
close all

format long

alignToSlipSystem0=1;
lattice='fcc';

if(alignToSlipSystem0)
switch lattice
    case 'bcc'
        s=[1 1 -1];
        n=[1 0 1];
    case 'fcc'
        s=[0 1 1];
        n=[-1 1 -1];
        n=n/norm(n);
    otherwise
end
        s=s/norm(s);
        n=n/norm(n);
        C2G1=[s;cross(n,s);n]
        %C2G1=[s;n;cross(s,n)]

    else
    C2G1=eye(3)
end

switch lattice
    case 'bcc'
        A=1/sqrt(3)*[-1 1 1;1 -1 1;1 1 -1]; % bcc
        material='W';
    case 'fcc'
        A=1/sqrt(2)*[0 1 1;1 0 1;1 1 0]; % fcc
%        material='AlMg15';
%        material='Al';
        material='Ni';
    otherwise
end
A=C2G1*A;
invA=inv(A);

% compute the three vectors that we determine to be along the crystal directions x1 x2 x3
% Find unit cell vectors along global axes
for i=1:3
Bi=invA(:,i);
Bi=Bi/max(abs(Bi));
[N,D] = rat(Bi,0.1);
Nr=N./D*prod(D);
g=gcd(gcd(Nr(1),Nr(2)),Nr(3));
Nr=Nr/g;
L(:,i)=A*Nr;
Ln(i)=norm(L(:,i));
end

%scaling=1000;
%T=round(scaling*diag(1./Ln))

% if you want to apply a shear to the box, then you have to populate the following matrix to transform the box
% The first column are the integer coefficients of vectors x1 x2 x3 of the box (the rows of L)
%T=[1000 0 0
%    0 500 0
%    0 0 1000]; % integer combinations of L columns

% T can do scaling or tilting, obtained as linear combinations of L
%
T = [1 0 0;
   5 1 0;
   0 0 1];
% T_ij is the integer coefficient that multiplies vector L_i in determining B_j

B=L*T % box , this scales the unit mesh,. T is a matrix of integers. this multiplication produces a linear combination of (column) vectors in L.
%B contains in columns the three box vectors

x0=[0 0 0];


fID=fopen('polycrystal.txt','w');
fprintf(fID,['materialFile=../../../MaterialsLibrary/' material '.txt; \n']);
fprintf(fID,'enablePartials=0; \n');
fprintf(fID,'absoluteTemperature = 300; # [K] simulation temperature \n');
fprintf(fID,'dislocationMobilityType=default; # default or FCC,BCC,HEXbasal,HEXprismatic,HEXpyramidal \n');
fprintf(fID,'meshFile=../../../MeshLibrary/unitCube.msh; \n');
fprintf(fID,'C2G1=');
fprintf(fID,'%1.15f ',C2G1(1,:));
fprintf(fID,'\n ');
fprintf(fID,'%1.15f ',C2G1(2,:));
fprintf(fID,'\n ');
fprintf(fID,'%1.15f ',C2G1(3,:));
fprintf(fID,';\n ');
fprintf(fID,'\n ');
fprintf(fID,'\n ');

fprintf(fID,'A=');
fprintf(fID,'%1.15f ',B(1,:));
fprintf(fID,'\n ');
fprintf(fID,'%1.15f ',B(2,:));
fprintf(fID,'\n ');
fprintf(fID,'%1.15f ',B(3,:));
fprintf(fID,';\n ');
fprintf(fID,'\n ');
fprintf(fID,'\n ');

fprintf(fID,'x0=');
fprintf(fID,'%1.15f ',x0(1,:));
fprintf(fID,';\n ');

fprintf(fID,'periodicFaceIDs=');
fprintf(fID,'%d ',[0 1 2 3 4 5]);
fprintf(fID,';\n ');

fprintf(fID,'\n#################\n');
fprintf(fID,'# GlidePlaneNoise\n');
fprintf(fID,'gridSize=');
fprintf(fID,'%d ',[256 256]);
fprintf(fID,'; # size of grid on the glide plane\n ');
fprintf(fID,'gridSpacing_SI=');
fprintf(fID,'%d ',[1 1]*1e-10);
fprintf(fID,'; # [m] spacing of grid on the glide plane\n ');

fprintf(fID,'solidSolutionNoiseMode=0; # 0=no noise, 1= read noise, 2=compute noise\n');
fprintf(fID,'solidSolutionNoiseFile_xz=../../../NoiseLibrary/noise_xz.vtk;\n');
fprintf(fID,'solidSolutionNoiseFile_yz=../../../NoiseLibrary/noise_yz.vtk;\n');
fprintf(fID,'stackingFaultNoiseMode=0; # 0=no noise\n');
fprintf(fID,'spreadLstress_A=1; # add comment\n');
fprintf(fID,'a_cai_A=1; # add comment\n');
fprintf(fID,'seed=0; # add comment\n');





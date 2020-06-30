function  RMass = Inertia(dof,Ine,Ainf);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function to insert the inertia array (mass and moment of inertia) % 
%    computed with AQUADYN code for a 7 DoF Wave Energy Converter    %
%                                                                    %
%  MARCO ALVES                                                       %
%  WAVEC                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ii=length(dof);
for l=1:1:ii
    kl=dof(l);
    for j=1:1:ii
        kj=dof(j); 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % NOTE: For more than 6DoF is adevisable %
        % to check the dimensionless method      %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        x=rem(l,6); if x==0; x=6; end
        y=rem(j,6); if y==0; y=6; end
        switch (x & y)
              case (x<=3 & y<=3), n=3;
              case (x<=3 & y>3)|(x>3 & y<=3), n=4;
              case (x>3 & y>3), n=5;
        end    
        DimAM=1000*(Ine(1,1))^n;
        mass(l,j)=(Ine(kl+1,2*kj));
        iadm(l,j)=(Ainf(kl,2*kj)*DimAM);
    end
end
RMass=mass;%+iadm;
function  Hst = HydroStatic(dof,HySt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Function to insert the Hydrostatic Coefficients computed   %            
%       with AQUADYN code for a 7 DoF Wave Energy Converter       %
%                                                                 %
%  MARCO ALVES                                                    %
%  WAVEC                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lr=HySt(1,1);
g=9.81;
wd=1000;
ii=length(dof);
%% The numerical code AQUADYN gives the hydrostatic coefficients array  %%
%% with dimensions or dimensionless. In the first case the value of     %%
%% constant "units" should be 1 and in the second case it should be 2.  %%
units=1;
if units==1, n1=1; n2=1; n3=1; 
   elseif units==2 
          n1=wd*g*Lr^2;
          n2=wd*g*Lr^3;
          n3=wd*g*Lr^4;
end
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
             case (x<3|x>5)&(y<3|y>5), dim=1;
             case (x==3 & y==3), dim=n1;
             case (x>2 & y>2) & (x==3 | y==3) & abs(x-y)<=2, dim=n2;
             case (x>3 & y>3) & (x<6 & y<6) & abs(x-y)<=1, dim=n3;
        end    
        HstC(l,j)=(HySt(kl+1,2*kj)*dim);
    end
end
Hst=[HstC];
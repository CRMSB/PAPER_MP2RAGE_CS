%
% [struct_out] = MP2RAGE_LookUpTable(IMs4D,struct_MP2RAGE)
%
% Author:   Aurelien TROTIER  (a.trotier@gmail.com)
% Date:     2020-03-18
% Partner: none
% Institute: CRMSB (Bordeaux)
%
% Function description:
%
%
% Input:
%       IMs4D : TI1 and TI2 images of the MP2RAGE with the 2 TIs along dimension 4
%       struct_MP2RAGE : MP2RAGE structure containing informations about sequences
%               .ETL = echo train length
%               .TI1 = inversion time 1 (ms)
%               .TI2 = inversion time 2 (ms)
%               .alpha1 = flip angle of the first readout (degree)
%               .alpha2 = flip angle of the second readout (degree)
%               .MP2RAGE_TR = delay of the sequence between 2 inversions pulses (ms)
%               .TR = delay between 2 excitation pulses during the readout (ms)
%       
% Output:
%
%
%
%
% Algorithm & bibliography:
%
% See also :
%
% To do :
%

function [struct_out] = MP2RAGE_LookUpTable(IMs4D,struct_MP2RAGE)

%% reco MP2RAGE 1 : combined channeld before images
    IM2MP=real(conj(IMs4D(:,:,:,1)).*IMs4D(:,:,:,2)./(abs(IMs4D(:,:,:,1)).^2+abs(IMs4D(:,:,:,2)).^2));




%% T1 map lookup-table
T1 = 0:10:10000;
lookuptable = param2LUT(T1,struct_MP2RAGE);

[maxValue] = max(lookuptable);
[~, I1] = find(lookuptable == maxValue,1,'first');
%indice du premier max
lookuptable(1:I1-1)=lookuptable(I1);    %we replace all the value before the max by max
[minValue] = min(lookuptable);            
[~,I2]=find(lookuptable==minValue,1,'first');          %indice of the first min after the max

T1 = T1(I1:I2);                         %on prend les T1 entre ces deux pics
lookuptable = lookuptable(I1:I2);

%% use lookup table on each pixel
T1map=0*IM2MP;

temp=IM2MP(:);
for i=1:length(temp)
    [~,I]=find(lookuptable < temp(i),1,'first');
    
    if isnan(I)==0
        T1map(i)=T1(I);
    end
end

struct_out.T1map = T1map;
struct_out.IM2MP = IM2MP;
struct_out.simu.lookUpTable = lookuptable;
struct_out.simu.T1 = T1;
end


%% sub function
function [lookuptable] = param2LUT(T1,struct_MP2RAGE)

    n = struct_MP2RAGE.ETL;
    TI1 = struct_MP2RAGE.TI1;
    TI2 = struct_MP2RAGE.TI2;
    alpha1 = struct_MP2RAGE.alpha1;
    alpha2 = struct_MP2RAGE.alpha2;
    MP2RAGETR = struct_MP2RAGE.MP2RAGE_TR;
    TR = struct_MP2RAGE.TR;

% check parameter
disp('-------------------------------------------------');
disp(['alpha1 :', num2str(alpha1), ' alpha2 :', num2str(alpha2)]);
disp(['TI1 :', num2str(TI1), ' TI2 :', num2str(TI2)]);
disp(['TR :', num2str(TR), ' MP2RAGETR :', num2str(MP2RAGETR), ' size of echo trains :', num2str(n)]);
disp('-------------------------------------------------');

E1=exp(-TR./T1);
EA=exp(-(TI1-(n./2-1).*TR)./T1);
EB=exp(-(TI2-TI1-n.*TR)./T1);
EC=exp(-(MP2RAGETR-(TI2+(n./2).*TR))./T1);

eff=1;

% compute mzss =[(B).*(cosd(alpha2).*E1).^n+A].*EC+(1-EC);

B = ((1-EA).*(cosd(alpha1).*E1).^n+(1-E1).*(1-(cosd(alpha1).*E1).^n)./(1-cosd(alpha1).*E1)).*EB+(1-EB);
A = (1-E1).*((1-(cosd(alpha2).*E1).^n)./(1-cosd(alpha2).*E1));

mzss_num=((B).*(cosd(alpha2).*E1).^n+A).*EC+(1-EC);
mzss_denom=(1+eff.*(cosd(alpha1).*cosd(alpha2)).^n .* exp(-MP2RAGETR./T1));

mzss=mzss_num./mzss_denom;

% compute GRE1= sind(alpha1).*(A.*(cosd(alpha1).*E1).^(n./2-1)+B)

A=-eff.*mzss.*EA+(1-EA);
B=(1-E1).*(1-(cosd(alpha1).*E1).^(n./2-1))./(1-cosd(alpha1).*E1);

GRE1=sind(alpha1).*(A.*(cosd(alpha1).*E1).^(n./2-1)+B);

% compute GRE2= sind(alpha2).*(A-B)

A=(mzss-(1-EC))./(EC.*(cosd(alpha2).*E1).^(n./2));
B=(1-E1).*((cosd(alpha2).*E1).^(-n./2)-1)./(1-cosd(alpha2).*E1);

GRE2=sind(alpha2).*(A-B);

lookuptable=GRE2.*GRE1./(GRE1.*GRE1+GRE2.*GRE2);

figure;plot(T1,GRE2.*GRE1./(GRE1.*GRE1+GRE2.*GRE2));title('MP2RAGE T1 bijection');
pause(1)

end

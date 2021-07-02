
function handle_connection_MP2RAGE_CS_allCase(connection)
tic
disp("handle_connection was called.")

next_acquisition = @connection.next;

sRO = connection.header.encoding.reconSpace.matrixSize.x;
sPE = connection.header.encoding.reconSpace.matrixSize.y;
sSE = connection.header.encoding.reconSpace.matrixSize.z;


nContrast = connection.header.encoding.encodingLimits.contrast.maximum;
nSet = connection.header.encoding.encodingLimits.set.maximum;

maxDim = max(nContrast, nSet) + 1; % Sequence of pierre use nSet and my sequence use nContrast
acquisition = next_acquisition(); % Call input function to produce the next acquisition.

nCh = size(acquisition.data.data,2);
kdata=zeros(sRO,nCh,sPE*sSE,maxDim);

row = acquisition.data.header.kspace_encode_step_1 + 1;
col = acquisition.data.header.kspace_encode_step_2 + 1;
%rep = acquisition.data.header.kspace_encode_step_2 + 1;
if nContrast > nSet
    TI_idx = acquisition.data.header.contrast + 1;
else
    TI_idx = acquisition.data.header.set + 1;
end

kdata(:,:,sub2ind([sPE,sSE,maxDim],row,col,TI_idx))=acquisition.data.data;

if (acquisition.ref.count > 0)
    kref=zeros(sRO,nCh,sPE*sSE,maxDim);
    row = acquisition.ref.header.kspace_encode_step_1 + 1;
    col = acquisition.ref.header.kspace_encode_step_2 + 1;
    
    if nContrast > nSet
        TI_idx = acquisition.ref.header.contrast + 1;
    else
        TI_idx = acquisition.ref.header.set + 1;
    end
    kref(:,:,sub2ind([sPE,sSE,maxDim],row,col,TI_idx))=acquisition.ref.data;
end

%% prepare parameter for reconstruction
if strcmp(connection.header.sequenceParameters.sequence_type,'Unknown') % sequence a_MP2RAGEPhase (developped by CRMSB)    
    kdata = permute(kdata,[1 3 2 4]);  
    kdata=reshape(kdata,sRO,sPE,sSE,nCh,1,[]); %% buffering the echo according to bart convetion [RO,E1,E2,CHA,MAP,CON]  
    
    struct_MP2RAGE.calibSize = connection.header.userParameters.userParameterDouble(3).value;
    struct_MP2RAGE.ETL = connection.header.userParameters.userParameterLong(6).value;
    struct_MP2RAGE.TI1 = connection.header.userParameters.userParameterLong(2).value;
    struct_MP2RAGE.TI2 = connection.header.userParameters.userParameterLong(3).value;
    struct_MP2RAGE.alpha1 = connection.header.sequenceParameters.flipAngle_deg(1);
    struct_MP2RAGE.alpha2 = connection.header.sequenceParameters.flipAngle_deg(1);
    struct_MP2RAGE.MP2RAGE_TR = connection.header.userParameters.userParameterLong(4).value;
    struct_MP2RAGE.TR = connection.header.sequenceParameters.TR;
    
else % sequence MP2RAGE developped by LIRYC (Pierre Bour)
    kdata = permute(kdata,[1 3 2 4]);
    kdata=reshape(kdata,sRO,sPE,sSE,nCh,1,[]); %% buffering the echo according to bart convetion [RO,E1,E2,CHA,MAP,CON]
    
    struct_MP2RAGE.calibSize = connection.header.userParameters.userParameterDouble(3).value;
    struct_MP2RAGE.ETL = connection.header.encoding.echoTrainLength;
    struct_MP2RAGE.TI1 = connection.header.sequenceParameters.TI(1);
    struct_MP2RAGE.TI2 = connection.header.sequenceParameters.TI(2);
    struct_MP2RAGE.alpha1 = connection.header.sequenceParameters.flipAngle_deg(1);
    struct_MP2RAGE.alpha2 = connection.header.sequenceParameters.flipAngle_deg(1);
    struct_MP2RAGE.MP2RAGE_TR = connection.header.sequenceParameters.TR;
    struct_MP2RAGE.TR = connection.header.sequenceParameters.echo_spacing;    
end

if exist('kref')
    kref = permute(kref,[1 3 2 4]);
    kref=reshape(kref,sRO,sPE,sSE,nCh,1,[]); %% buffering the echo according to bart convetion [RO,E1,E2,CHA,MAP,CON]
else
    sy1 = size(kdata,2)/2-struct_MP2RAGE.calibSize/2+1;
    sy2 = size(kdata,2)/2+struct_MP2RAGE.calibSize/2;
    sz1 = size(kdata,3)/2-struct_MP2RAGE.calibSize/2+1;
    sz2 = size(kdata,3)/2+struct_MP2RAGE.calibSize/2;
    
    kref(:,sy1:sy2,sz1:sz2,:,:,:) = kdata(:,sy1:sy2,sz1:sz2,:,:,:);
end
toc
if 0 % plot mask
    mask = zeros(size(kref,1),size(kref,2),size(kref,3));
    mask(abs(kdata(:,:,:,1)) > 0) = 1;
    imshow(squeeze(mask(1,:,:)),[]);
end


%% Estimate coil sensitivity
sensitivity_coil_map=bart(['caldir ' num2str(struct_MP2RAGE.calibSize)], kdata(:,:,:,:,:,2));
%
%sensitivity_coil_map=bart(['ecalib -m1 -r ' num2str(struct_MP2RAGE.calibSize)], kref(:,:,:,:,:,2));

%% Parallel + CS
struct.bg_mult=0;

disp('-------------------------------------------------');
disp('********** bart(pics...)  **********');

img_combined_CS = zeros(size(kdata,1),size(kdata,2),size(kdata,3),maxDim);
for i = 1:maxDim
    img_combined_CS(:,:,:,i) = bart('pics -S -R W:7:0:0.01', kdata(:,:,:,:,:,i), sensitivity_coil_map);
end

[struct_T1map] = gadgetron.custom.MP2RAGE_CS.MP2RAGE_LookUpTable(img_combined_CS,struct_MP2RAGE); %% comb then MP2RAGE
disp('-----------RECO WITH correction background---------------');

multiFactor=struct.bg_mult*mean(mean(mean(abs(img_combined_CS(1:end,end-10:end,end-10:end,2)))));
MP2RAGE_CS=real((conj(img_combined_CS(:,:,:,1)).*img_combined_CS(:,:,:,2)-multiFactor)./(abs(img_combined_CS(:,:,:,1)).^2+abs(img_combined_CS(:,:,:,2)).^2+2*multiFactor));


%% Prepare image to send
img_to_send{1}=scale_image_for_gt_and_scanner(permute(abs(img_combined_CS(:,:,:,1)),[4, 1, 2, 3]));
img_to_send{2}=scale_image_for_gt_and_scanner(permute(abs(img_combined_CS(:,:,:,2)),[4, 1, 2, 3]));
img_to_send{3}=scale_MP2RAGE_for_gt_and_scanner(permute(MP2RAGE_CS,[4, 1, 2, 3]));
img_to_send{4}=permute(abs(struct_T1map.T1map),[4, 1, 2, 3]); % no scaling but T1 range limitated to 0 to 4095 (dicom format)

%% send image

for ii=1:length(img_to_send)
    image = gadgetron.types.Image.from_data(single(img_to_send{ii}), reference_header(acquisition));
    image.header.image_type = gadgetron.types.Image.MAGNITUDE;
    image = gadgetron.custom.utils.correct_image_siemens(image,acquisition,connection);
    image.header.image_series_index = (ii-1)*1000+5000;

    disp("Sending image to client.");
    connection.send(image);
end
end

%% suppport functions
function reference = reference_header(acquisition)
    % We pick the first header from the header arrays - we need it to initialize the image meta data.    
    reference = structfun(@(arr) arr(:, 1)', acquisition.data.header, 'UniformOutput', false);
end

function [ image_scale  ] = scale_MP2RAGE_for_gt_and_scanner( im_in )  
  disp([' before scaling max : ' , num2str(max(im_in(:))), ' min : ' , num2str(min(im_in(:)))]);
  
  image_scale=im_in*4095+4095/2;
  
  disp([' after scaling max : ' , num2str(max(image_scale(:))), ' min : ' , num2str(min(image_scale(:)))]);
  
end

function [ image_scale  ] = scale_image_for_gt_and_scanner( im_in )  
    disp([' before scaling max : ' , num2str(max(im_in(:))), ' min : ' , num2str(min(im_in(:)))]);
    
    image_scale=(im_in-min(im_in(:)))*(2^12-1)/(max(im_in(:))-min(im_in(:)));
    
    disp([' after scaling max : ' , num2str(max(image_scale(:))), ' min : ' , num2str(min(image_scale(:)))]);
    
  end
  


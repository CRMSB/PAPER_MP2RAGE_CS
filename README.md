# Compressed Sensing-MP2RAGE

This repository contains the reconstruction code of this [paper](https://journals.lww.com/investigativeradiology/Abstract/9000/The_Compressed_Sensing_MP2RAGE_as_a_Surrogate_to.98641.aspx) accepted in *Investigative Radiology* :

If you use this code please cite:

>Trotier AJ, Dilharreguy B, Anandra S, et al. The Compressed Sensing MP2RAGE as a Surrogate to the MPRAGE for Neuroimaging at 3 T. Investigative Radiology 2022 doi: 10.1097/RLI.0000000000000849.


## What is this repository for?
This repository contains the code necessary to reconstruct the images and also the corresponding T1 map from the MP2RAGE sequence accelerated with CS.
This repository contains the matlab code, the dedicated gadgetron gadgets and configuration files as well as a test script.

## How to use it ?
We wrap our tool on a docker image, so there is no need to install any library dependencies, the only requirement is to have docker, nvidia-docker and CUDA installed.
It is possible to modify the dockerfile/launch_docker script in order to remove the gpu requirements.

### Prerequisites:
- CUDA
- Docker (For running on CPU) (https://docs.docker.com/install/)
- NVIDIA-Docker (For running on GPU ) (https://github.com/nvidia/nvidia-docker/wiki)
- Matlab need to be installed on the Host

### Tested Configuration:

* Ubuntu 18.04
* MATLAB R2019b

# Installation

* Clone/download the git repository to your local machine
* Change the line in the **Dockerfile** : `ARG MATLAB_VERSION=R2019b` for your matlab version
* Create the docker image : 
```bash
cd PAPER_MP2RAGE_CS/docker
sudo docker build -t gadgetronpaper .
```
* Edit the **launch_docker.sh** according to your preferences :
  *   `-v /data/dumpSiemens:/data/dumpSiemens \` corresponds to the volume mounted on the host to the container
  *    ` -v /usr/local/MATLAB:/usr/local/MATLAB \` mounts the matlab folder from the host in the container
* launch the script :`sh launch_docker.sh` 

**in the created container**
* check if matlab is supported for gadgetron with :`gadgetron --info` 
* launch matlab and increase the 'java heap memory' in Preferences/General

## Data test:

* Knee / phantom datasets are available here : https://zenodo.org/record/5059979
* Initial datasets used in the paper can be obtained on request.


# Reconstruction procedure:

### Retrospective reconstruction from .h5
* Copy the dataset to your mounted folder 
* launch the container from a terminal : `docker exec -it <container> bash`
* start gadgetron : `gadgetron -p9002`
* From another terminal : `docker exec -it <container> bash`
* Run the reconstruction : `gadgetron_ismrmrd_client -c MP2RAGE_CS_Bucket.xml -f /data/dumpSiemens/<dataset.h5> -o /data/dumpSiemens/RECO.h5`
* You can visualize the reconstruction with ismrmrdviewer : https://github.com/ismrmrd/ismrmrdviewer or by loading the reconstruction under matlab (don't forget to include the matlab folder of the ismrmrd library):

```matlab
function [img] = ismrmrd_read_img(struct)
%This function read ismrmrd out.h5 image
%   input :
%           struct.path : string

if nargin < 1 || ~isfield(struct,'path')
    [file,path]=uigetfile('*.h5');
    struct.path = fullfile(path,file);
end

S=hdf5info(struct.path);

if exist(struct.path, 'file')
    dset = ismrmrd.Dataset(struct.path, 'dataset');
else
    error(['File ' struct.path ' does not exist.  Please generate it.'])
end

%% old way
disp(dset.fid.identifier)

S=hdf5info(struct.path);

    for i=1:length(S.GroupHierarchy(1).Groups(1).Groups)

        attributes=S.GroupHierarchy(1).Groups(1).Groups(i).Datasets(1).Name;
        dataset=S.GroupHierarchy(1).Groups(1).Groups(i).Datasets(2).Name;
        header=S.GroupHierarchy(1).Groups(1).Groups(i).Datasets(3).Name;
        tmp=hdf5read(struct.path,dataset);

        img(:,:,:,i)=squeeze(tmp);
    end
end
```
### From the scanner
If you have access to the sequence,
* launch the image : `docker exec -it <container> bash`
* start gadgetron on the dedicated port (generally 9002, see gadgetron.ini on scanner) : `gadgetron -p9002`
* If you want to dump the dataset, it will be stored on : `/data/dumpSiemens` according to the path in `/gadgetron-gadgets-MP2RAGE/MP2RAGE_CS_Bucket_DUMP.xml`


# Folder Structure:


`/docker/` - contains the dockerfile necessary to generate the docker image and a script to start a container with the appropriate parameters
`/gadgetron-gadgets-MP2RAGE/` - gadgetron gadgets / config which is installed along gadgetron in the image.
`/gadgetron-matlab-MP2RAGE/` - gadgetron-matlab interface from : https://github.com/gadgetron/gadgetron-matlab + custom matlab functions for MP2RAGE. Custom functions are stored under : `+gadgetron/+custom/`
# License

See LICENSE.txt



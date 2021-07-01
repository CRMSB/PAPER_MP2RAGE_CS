FROM nvidia/cuda:11.0-cudnn8-devel-rc

# Proxy issue ->
# ENV http_proxy="http://192.168.10.1:3128"

RUN apt-get update --quiet && \
    apt-get install --no-install-recommends --no-install-suggests --yes  \
    software-properties-common apt-utils wget build-essential cython3  python3-dev python3-pip libhdf5-serial-dev cmake git-core libboost-all-dev libfftw3-dev h5utils jq hdf\
5-tools liblapack-dev libatlas-base-dev libxml2-dev libfreetype6-dev pkg-config libxslt-dev libarmadillo-dev libace-dev gcc-multilib libgtest-dev liblapacke-dev libplplot-dev libdcmtk-dev sup\
ervisor supervisor net-tools cpio libpugixml-dev jove libopenblas-base libopenblas-dev libpugixml-dev vim cmake-curses-gui

RUN apt-get install \
    make \
    libhdf5-dev

RUN pip3 install -U pip setuptools
RUN pip3 install numpy scipy Cython tk-tools matplotlib scikit-image opencv_python pydicom scikit-learn sympy
RUN apt-get install --no-install-recommends --no-install-suggests --yes python3-psutil python3-pyxb python3-lxml python3-pil python3-h5py
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends --no-install-suggests --yes python3-tk

RUN pip3 install torch==1.4.0+cpu torchvision==0.5.0+cpu -f https://download.pytorch.org/whl/torch_stable.html

# for embedded python plot, we need agg backend
RUN mkdir -p /root/.config/matplotlib && touch /root/.config/matplotlib/matplotlibrc 
RUN echo "backend : agg" >> /root/.config/matplotlib/matplotlibrc

# Install MKL  
RUN cd /tmp && \
   wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB && \
   apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB && \
   sh -c 'echo deb https://apt.repos.intel.com/mkl all main > /etc/apt/sources.list.d/intel-mkl.list' && \
   apt-get update && \
   apt-get install --yes intel-mkl-64bit-2018.2-046 && \
   update-alternatives --install /usr/lib/x86_64-linux-gnu/libblas.so     libblas.so-x86_64-linux-gnu      /opt/intel/mkl/lib/intel64/libmkl_rt.so 150 && \
   update-alternatives --install /usr/lib/x86_64-linux-gnu/libblas.so.3   libblas.so.3-x86_64-linux-gnu    /opt/intel/mkl/lib/intel64/libmkl_rt.so 150 && \
   update-alternatives --install /usr/lib/x86_64-linux-gnu/liblapack.so   liblapack.so-x86_64-linux-gnu    /opt/intel/mkl/lib/intel64/libmkl_rt.so 150 && \
   update-alternatives --install /usr/lib/x86_64-linux-gnu/liblapack.so.3 liblapack.so.3-x86_64-linux-gnu  /opt/intel/mkl/lib/intel64/libmkl_rt.so 150 && \
   echo "/opt/intel/lib/intel64"     >  /etc/ld.so.conf.d/mkl.conf && \
   echo "/opt/intel/mkl/lib/intel64" >> /etc/ld.so.conf.d/mkl.conf && \
   ldconfig && \
   echo "MKL_THREADING_LAYER=GNU" >> /etc/environment

#ZFP
RUN mkdir /opt/code
RUN cd /opt/code && \
    git clone https://github.com/hansenms/ZFP.git && \
    cd ZFP && \
    mkdir lib && \
    make && \
    make shared && \
    make -j $(nproc) install

#Set more environment variables in preparation for Gadgetron installation
ENV GADGETRON_HOME=/usr/local \
    ISMRMRD_HOME=/usr/local

ENV PATH=$PATH:$GADGETRON_HOME/bin:$ISMRMRD_HOME/bin:/usr/local/cuda/bin \
    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ISMRMRD_HOME/lib:$GADGETRON_HOME/lib:/usr/local/lib64 \
    LIBRARY_PATH=$LIBRARY_PATH:/usr/local/cuda/lib64/stubs

#ISMRMRD
ARG ISMRMRD_URL=https://github.com/ismrmrd/ismrmrd.git
ARG ISMRMRD_COMMIT=27a7c58

RUN echo ${ISMRMRD_URL}
RUN cd /opt/code && \
    pwd && ls -l /opt/ && \
    git clone ${ISMRMRD_URL} && \
    echo ${ISMRMRD_COMMIT} && \
    pwd && ls -l /opt/code && \
    cd ismrmrd && \ 
    git checkout ${ISMRMRD_COMMIT} && \
    mkdir build && \
    cd build && \
    cmake ../ && \
    make -j $(nproc) && \
    make install

#GADGETRON
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda:/usr/local/cuda/lib64 
RUN ls -l /usr/local && ls -l /usr/local/cuda/
RUN nvcc -V
RUN echo ${LD_LIBRARY_PATH}
RUN ls -l /usr/local
RUN echo ${PATH}
ARG GADGETRON_URL=https://github.com/gadgetron/gadgetron.git
ARG GADGETRON_COMMIT=7bbf633
RUN cd /opt/code && \
    git clone ${GADGETRON_URL} && \
    cd gadgetron && \
    git checkout ${GADGETRON_COMMIT} && \
    mkdir build && \
    cd build && \
    cmake .. -DUSE_MKL=ON -DBUILD_SUPPRESS_WARNINGS=ON && \
    make -j $(nproc) && \
    make install && \
    /opt/code/gadgetron/docker/manifest --key .io.gadgetron.gadgetron.sha1 --value `git rev-parse HEAD` && \
    cp /opt/code/gadgetron/docker/start_supervisor /opt/ && \
    cp /opt/code/gadgetron/docker/supervisord.conf /opt/

#SIEMENS_TO_ISMRMRD
ARG SIEMENS_TO_ISMRMRD_URL=https://github.com/ismrmrd/siemens_to_ismrmrd.git
ARG SIEMENS_TO_ISMRMRD_COMMIT=11cf24b
RUN cd /opt/code && \
    git clone ${SIEMENS_TO_ISMRMRD_URL} && \
    cd siemens_to_ismrmrd && \
    git checkout ${SIEMENS_TO_ISMRMRD_COMMIT} && \
    mkdir build && \
    cd build && \
    cmake ../ && \
    make -j $(nproc) && \
    make install && \
    /opt/code/gadgetron/docker/manifest --key .io.gadgetron.siemens_to_ismrmrd.sha1 --value `git rev-parse HEAD` 

#BART
ARG BART_URL=https://github.com/mrirecon/bart
ARG BART_COMMIT=e1a4303
RUN cd /opt/code && \
    git clone ${BART_URL} && \
    cd bart && \
    git checkout ${BART_COMMIT} && \
    cd Makefiles && \
    > Makefile.local && \
    echo "CUDA=1" >> Makefile.local && \
    echo "CUDA_BASE = /usr/local/cuda/" >> Makefile.local && \
    echo "MKL=1" >> Makefile.local && \
    echo "MKL_BASE=/opt/intel/mkl" >> Makefile.local && \
    cd .. && \
    make && \
    make install

# Mount matlab folder into docker MATLAB 
ARG MATLAB_VERSION=R2019b
RUN ln -sf /usr/local/MATLAB/${MATLAB_VERSION}/bin/matlab /usr/local/bin/matlab
#ENV PATH=$PATH:/usr/local/MATLAB/${MATLAB_VERSION}/bin
RUN apt-get install --yes libxt6 default-jre



# clone repository
RUN cd /opt/code/  && \
    git clone https://github.com/CRMSB/PAPER_MP2RAGE_CS.git
    

RUN mkdir -p /data/dumpSiemens/ && \ 
    ln -s /opt/intel/mkl/lib/intel64/*.so /usr/local/lib/ && \
    cd /opt/code/PAPER_MP2RAGE_CS/gadgetron-gadgets-MP2RAGE/ && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make -j $(nproc) && \
    make install


# GADGETRON-MATLAB (branch docker)
RUN cd /opt/code && \
    mkdir -p /root/.matlab/${MATLAB_VERSION} && \
    cd /root/.matlab/${MATLAB_VERSION} && \
    touch javaclasspath.txt && \ 
    echo "/opt/code/PAPER_MP2RAGE_CS/gadgetron-matlab-MP2RAGE/java" >> javaclasspath.txt

# Gadgetron-gadget (required gadgetron installation)

RUN mkdir -p /data/dumpSiemens/ && \ 
    cd /opt/code/PAPER_MP2RAGE_CS/gadgetron-gadgets-MP2RAGE && \
    mkdir -p build && \
    cd build && \
    cmake .. && \
    make -j $(nproc) && \
    make install




# Clean up packages.
RUN  apt-get clean && \
   rm -rf /var/lib/apt/lists/*

#CMD ["/opt/start_supervisor"]
CMD ["/bin/bash"]

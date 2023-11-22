FROM ubuntu:22.04

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
#RUN apt-get install --no-install-recommends --no-install-suggests --yes python3-psutil python3-pyxb python3-lxml python3-pil python3-h5py
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends --no-install-suggests --yes python3-tk

#RUN pip3 install torch==1.4.0+cpu torchvision==0.5.0+cpu -f https://download.pytorch.org/whl/torch_stable.html

# for embedded python plot, we need agg backend
RUN mkdir -p /root/.config/matplotlib && touch /root/.config/matplotlib/matplotlibrc 
RUN echo "backend : agg" >> /root/.config/matplotlib/matplotlibrc

# Install MKL  

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

ENV PATH=$PATH:$GADGETRON_HOME/bin:$ISMRMRD_HOME/bin \
    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ISMRMRD_HOME/lib:$GADGETRON_HOME/lib:/usr/local/lib64

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
    cmake ..  -DBUILD_SUPPRESS_WARNINGS=ON && \
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

ENV MATLABPATH=$MATLABPATH:/opt/code/PAPER_MP2RAGE_CS/gadgetron-matlab-MP2RAGE:/opt/code/bart/matlab
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

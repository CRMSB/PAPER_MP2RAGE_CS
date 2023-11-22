#!/bin/bash
docker run --name=gt_paper_mp2rage_cs \
  --tty --interactive \
  --net=host \
  --privileged \
  -v /data/dumpSiemens:/data/dumpSiemens \
  -v /usr/local/MATLAB:/usr/local/MATLAB \
  -v /tmp/.X11-unix:/tmp/.X11-unix \
  --volume="$HOME/.Xauthority:/root/.Xauthority:rw" \
  -e DISPLAY \
  gadgetronpaper
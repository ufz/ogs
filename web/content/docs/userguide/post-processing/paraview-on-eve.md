+++
date = "2018-11-14T14:00:13+01`:00"
title = "ParaView on Eve frontends / envinf1"
author = "Lars Bilke"
weight = 1

[menu]
  [menu.userguide]
    parent = "post-processing"

aliases = ['paraview-on-envinf1']
+++

## Eve frontends

Connect to a frontend node and start `pvserver` via singularity:

```bash
ssh frontend1.eve.ufz.de -L 11111:frontend1.eve.ufz.de:11111
# OR: ssh frontend2.eve.ufz.de -L 11111:frontend2.eve.ufz.de:11111
# frontend2 has a NVidia K80 vs a P4000 in frontend1
singularity exec --nv -B /data:/data /data/ogs/images/pv-v5.8.0-egl-py2.sif /opt/paraview/bin/pvserver
```

On your local machine:

* [Install ParaView 5.8.0](https://www.paraview.org/download/)
* Connect to server
  * Add Server
    * Host: `localhost`
    * Port: `11111`
* Under `Settings / Render View` set `Remote Render Threshold` to a small value (e.g. 1) to ensure remote rendering

The port tunneling with ssh is required as these ports are blocked from the firewall. If port `11111` is already in use by another user just try a different port, e.g. `11112`: `pvserver -sp=11112`. Do not forget to tunnel this port with SSH too!

## `envinf1`

`envinf1` has a NVidia K20m. Same instructions as for eve but port tunneling is not necessary:

```bash
ssh 141.65.34.100
singularity exec --nv /data/shared/container/pv-v5.8.0-egl-py2.sif /opt/paraview/bin/pvserver
```

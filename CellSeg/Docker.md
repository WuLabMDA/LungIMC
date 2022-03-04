## 1. Docker image preparation
* Build from Dockerfile
```
$ docker build -t deepcell_imc:cpath .
```

## 2. Docker container running
```
$ docker run -it --rm --user $(id -u):$(id -g)  \
  -v /rsrch1/ip/pchen6/Codes/LungIMC:/App/LungIMC \
  -v /rsrch1/ip/pchen6/LungIMCData:/Data \
  --shm-size=128G --gpus '"device=2"' --cpuset-cpus=20-39 \
  --name imc_seg deepcell_imc:cpath
```
```
$ docker exec -it imc_seg /bin/bash
```

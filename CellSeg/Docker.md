## 1. Docker image preparation
* Build from Dockerfile
```
$ docker build -t deepcell_imc:lung .
```

## 2. Docker container running
```
$ docker run -it --rm --user $(id -u):$(id -g)  \
  -v /rsrch1/ip/pchen6/Codes/LungIMC:/App/LungIMC \
  -v /rsrch1/ip/pchen6/LungIMCData:/Data \
  --shm-size=128G --gpus '"device=0"' --cpuset-cpus=0-9 \
  --name deepcell_imc_lung deepcell_imc:lung
```
```
$ docker exec -it deepcell_imc_lung /bin/bash
```

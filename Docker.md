## 1. Docker image preparation
* Build from Dockerfile
```
$ docker build -t lungimc:ping .
```
* Tag built image for the Docker Hub
```
$ docker tag lungimc:ping pingjunchen/lungimc:ping
```


## 2. Docker container running
```
$ docker run -it --rm --user $(id -u):$(id -g) \
  -v /rsrch1/ip/pchen6/Codes/LungIMC:/App/LungIMC \
  -v /rsrch1/ip/pchen6/LungIMCData:/Data \
  --shm-size=240G --gpus '"device=1"' --cpuset-cpus=0-39 \
  --name lungimc_ping pingjunchen/lungimc:ping 
```

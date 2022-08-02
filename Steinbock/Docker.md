## 1. Docker image preparation
* Build from Dockerfile
```
$ docker build -t lungimc:r .
```

## 2. Docker container running
```
$ docker run -it --rm --user $(id -u):$(id -g) \
  -v /rsrch1/ip/pchen6/Codes/LungIMC:/App/LungIMC \
  -v /rsrch1/ip/pchen6/LungIMCData:/Data \
  --shm-size=896G --cpuset-cpus=100-255 \
  --name lungimc_r lungimc:r
```

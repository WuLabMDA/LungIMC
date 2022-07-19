## 1. Docker image preparation
* Build from Dockerfile
```
$ docker build -t lungimc:ping .
```

## 2. Docker container running
```
$ docker run -it --rm --user $(id -u):$(id -g) \
  -v /rsrch1/ip/pchen6/Codes/LungIMC:/App/LungIMC \
  -v /rsrch1/ip/pchen6/LungIMCData:/Data \
  --shm-size=228G --cpuset-cpus=0-39 \
  --name lungimc_ping lungimc:ping
```

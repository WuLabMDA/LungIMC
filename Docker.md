## 1. Docker image preparation
* Build from Dockerfile
```
$ docker build -t lungimc:ping .
```

## 2. Docker container running
```172.30.205.155
$ docker run -it --rm --user $(id -u):$(id -g) \
  -p 14321:14321 --hostname 1mcprddgx05 \
  -v /rsrch1/ip/pchen6/Codes/LungIMC:/App/LungIMC \
  -v /rsrch1/ip/pchen6/LungIMCData:/Data \
  --shm-size=640G --cpuset-cpus=150-255 \
  --name lungimc_ping lungimc:ping
```


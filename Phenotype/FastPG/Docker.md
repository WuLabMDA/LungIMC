### Run FastPG Docker Container
```
$ docker run -it --rm --user $(id -u):$(id -g) \
  -v /rsrch1/ip/pchen6/LungIMCData:/Data \
  -v /rsrch1/ip/pchen6/Codes/LungIMC:/Code/LungIMC \
  -w /Code/LungIMC/Phenotype/FastPG \
  --shm-size=240G --cpuset-cpus=0-39 \
  --name fastpg_ping jefferys/fastpg:latest
```

# Docker

Note: You might need to do this with `sudo`!

## Create docker container

```
$ ocker build -t <name[:tag]> - < Dockerfile
$ docker images
$ docker tag <IMAGE ID> <username or organization>/<repository>:<tag>
```

Example (uefbic/bic_gwas):
```
$ docker build -t bic_gwas - < bic_gwas.dockerfile
$ docker images

REPOSITORY           TAG                IMAGE ID            CREATED             SIZE
bic_gwas             latest             ce58de56bd33        2 seconds ago       460MB

$ docker tag ce58de56bd33 bicuef/bic_gwas:latest
$ docker images

REPOSITORY           TAG                IMAGE ID            CREATED             SIZE
bic_gwas             latest             ce58de56bd33        2 seconds ago       460MB
bicuef/bic_gwas      latest             ce58de56bd33       12 seconds ago       460MB
```

## Push image to Docker Hub (hub.docker.com)

```
$ docker images
$ docker login --username <username> docker.io
$ docker push <username or organization>/<repository>:<tag>
```

Example (uefbic/bic_gwas):
```
$ docker images

REPOSITORY           TAG                IMAGE ID            CREATED             SIZE
bic_gwas             latest             ce58de56bd33        2 seconds ago       460MB
bicuef/bic_gwas      latest             ce58de56bd33       12 seconds ago       460MB

$ docker login --username kuulasmaa docker.io

Password:
WARNING! Your password will be stored unencrypted in /home/teemu/.docker/config.json.
Configure a credential helper to remove this warning. See
https://docs.docker.com/engine/reference/commandline/login/#credentials-store

Login Succeeded

$ docker push bicuef/bic_gwas:latest
$ docker tag ce58de56bd33 bicuef/bic_gwas:latest
```

## Pull image from Docker Hub (hub.docker.com)

```
# required if pulling from private repository
$ docker login --username <username> docker.io
$ docker pull <username or organization>/<repository>:<tag>
```

Example (uefbic/bic_gwas):
```
# required if pulling from private repository
$ docker login --username kuulasmaa docker.io

Password:
WARNING! Your password will be stored unencrypted in /home/teemu/.docker/config.json.
Configure a credential helper to remove this warning. See
https://docs.docker.com/engine/reference/commandline/login/#credentials-store

Login Succeeded

$ docker pull bicuef/bic_gwas:latest

latest: Pulling from bicuef/bic_gwas
5c939e3a4d10: Already exists
c63719cdbe7a: Already exists
19a861ea6baf: Already exists
651c9d2d6c4f: Already exists
4bebf4e39cf6: Pull complete
869603b42c12: Pull complete
e5cb5341effa: Pull complete
78ef512cae80: Pull complete
254f1861480d: Pull complete
2ba0243c4fb7: Pull complete
209fd9efad73: Pull complete
84b0c200af1e: Pull complete
bd8a3bc91330: Pull complete
aef5ad0ca86e: Pull complete
459cba61a1fa: Pull complete
2601c96efb9d: Pull complete
93b7f12a17bb: Pull complete
95873088eb8b: Pull complete
859c2f1464bb: Pull complete
4b1b0ba452cd: Pull complete
418af688d307: Pull complete
9c6772bad5b8: Pull complete
83adc924bdea: Pull complete
76d748e3b622: Pull complete
f6a2ff49e318: Pull complete
00ea609722bb: Pull complete
760e8c2778ad: Pull complete
f333a9857443: Pull complete
d97b9eb418e1: Pull complete
39429f8a02e2: Pull complete
cee3f31ab4fb: Pull complete
d55a3e36d23e: Pull complete
1966a8969e0f: Pull complete
2d93e2b819ef: Pull complete
2044c72cf398: Pull complete
e991d0d08965: Pull complete
7dd7468b6fb4: Pull complete
Digest: sha256:01cc56e7c8dfb8bf22dc73fa1c4cb2cdd1235035a3bd5744bac87d1971f50108
Status: Downloaded newer image for bicuef/bic_gwas:latest
docker.io/bicuef/bic_gwas:latest
```

## Execute command from container

```
$ docker run -i --rm <image> <command>
```


## Convert docker to singularity

```
$ singularity build <Singularity Image File .sif> docker://<username or organization>/<repository>:<tag>
```

Example (uefbic/bic_gwas):

```
$ singularity build bic_gwas.sif docker://uefbic/bic_gwas:latest
```

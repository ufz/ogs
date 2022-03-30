+++
date = "2018-02-26T11:00:13+01:00"
title = "Docker"
author = "Lars Bilke"
weight = 1050

aliases = ["/docs/devguide/advanced/docker"]

[menu]
  [menu.devguide]
    parent = "environments"
+++

## Introduction

> Docker containers wrap up a piece of software in a complete filesystem that contains everything it needs to run: code, runtime, system tools, system libraries – anything you can install on a server. This guarantees that it will always run the same, regardless of the environment it is running in.
>
> <cite>[www.docker.com/whatisdocker](https://www.docker.com/whatisdocker)</cite>

See the [docs](https://docs.docker.com) for installation instructions (if you are on Windows we highly recommend the [Docker Desktop WSL 2 backend](https://docs.docker.com/docker-for-windows/wsl/)).

## Images

Docker images can be created by [Dockerfiles](https://docs.docker.com/reference/builder/) which can be stacked and thus form a directed graph. OGS-6 image definitions are created with [ufz/ogs-container-maker](https://github.com/ufz/ogs-container-maker). Built images can be found at the [GitLab Container Registry](https://gitlab.opengeosys.org/ogs/ogs/container_registry). You can also [create images from your local source code](https://github.com/ufz/ogs-container-maker#build-ogs-from-local-git-repo).

To build an image by yourself create a Dockerfile:

```bash
FROM ubuntu:20.04

RUN ...
```

Run the `build` command:

```bash
docker build -t repo/image_name path/to/directory
```

- `-t` specifies a name for the image, can be arbitrary chosen (but should match the corresponding image on Docker Hub if there is one)
- The path should specify the directory where the Dockerfile is located

Now you can see your build image with `$ docker images`.

## Run a container

To run commands inside a container:

```bash
docker run --rm image_name command_to_run
```

- `--rm` Cleans up after exiting the container

To run an interactive shell add the `-it`-switch:

```bash
docker run --rm -it image_name
```

It is useful to mount folders from the host operating system in the Docker container, e.g. to edit source code on your host with your favorite editor:

```bash
docker run --rm -it -v /host/directory:/container/directory image_name
```

## Prebuilt OGS-6 Docker images

There are docker images provided on the [GitLab Container Registry](https://gitlab.opengeosys.org/ogs/ogs/container_registry) which include everything necessary to build OGS-6. They are used by the CI but you can also use them for development.

Even better for development is the usage of [Singularity container]({{< relref "singularity.md" >}}) because they offer a transparent mapping of the user to container space.

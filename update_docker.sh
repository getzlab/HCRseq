#!/bin/bash

docker buildx build --platform=linux/amd64 -t gcr.io/broad-getzlab-fmhcrsparc/hcrseq:v0.9 .
docker push gcr.io/broad-getzlab-fmhcrsparc/hcrseq:v0.9


#!/bin/bash

docker buildx build --platform=linux/amd64 -t gcr.io/broad-getzlab-fmhcrsparc/hcrseq:v1.0 .
docker push gcr.io/broad-getzlab-fmhcrsparc/hcrseq:v1.0


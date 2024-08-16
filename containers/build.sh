#! /usr/bin/env bash
cd ${BASH_SOURCE%/*}/..
docker build -t cracofunew:latest -f containers/Dockerfile .
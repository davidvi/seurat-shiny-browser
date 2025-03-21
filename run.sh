#!/bin/bash
# Get data folder path from command line argument or use default
DATA_FOLDER=${1:-"/home/shiny-app/data"}

# Run the Docker container
docker run -p 3030:3030 \
  --platform linux/x86_64 \
  -v /Users/vanijzen/Developer/shiny-docker/data:${DATA_FOLDER} \
  shiny-docker ${DATA_FOLDER}
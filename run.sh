#!/bin/bash
# Get data folder path and port from command line arguments or use defaults
DATA_FOLDER=${1:-"/home/shiny-app/data"}
PORT=${2:-3030}

# Run the Docker container, mapping the specified port
docker run -p ${PORT}:${PORT} \
  --platform linux/x86_64 \
  -v /Users/vanijzen/Developer/shiny-docker/data:${DATA_FOLDER} \
  shiny-docker ${DATA_FOLDER} ${PORT}
#!/bin/bash
# Get data folder path and port from command line arguments or use defaults
DATA_FOLDER=${1:-"/home/shiny-app/data"}
PORT=${2:-3030}

# Get the local data directory path - default to ./data if not specified
LOCAL_DATA_DIR=${3:-"$(pwd)/data"}

echo "Using local data directory: $LOCAL_DATA_DIR"
echo "Mounting to Docker container at: $DATA_FOLDER"
echo "Using port: $PORT"

# Run the Docker container, mapping the specified port
docker run -p ${PORT}:${PORT} \
  --platform linux/x86_64 \
  -v ${LOCAL_DATA_DIR}:${DATA_FOLDER} \
  shiny-docker ${DATA_FOLDER} ${PORT}
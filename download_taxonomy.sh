#!/bin/bash

# Check if two arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <output_file> path/to/<extract_dir>"
    exit 1
fi

# Name of the output file (from argument)
OUTPUT_FILE="$1"

# Directory to extract files to (from argument)
EXTRACT_DIR="$2"

# URL of the file to download
URL="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"

# Remove existing file if it exists
if [ -f "$OUTPUT_FILE" ]; then
    echo "Removing existing $OUTPUT_FILE"
    rm "$OUTPUT_FILE"
fi

# Download the file using wget
echo "Downloading file..."
wget -O "$OUTPUT_FILE" "$URL"

# Check if the download was successful
if [ $? -eq 0 ]; then
    echo "File downloaded successfully as $OUTPUT_FILE"

    # Remove existing directory if it exists
    if [ -d "$EXTRACT_DIR" ]; then
        echo "Removing existing $EXTRACT_DIR directory"
        rm -rf "$EXTRACT_DIR"
    fi

    # Create extraction directory
    mkdir -p "$EXTRACT_DIR"

    # Extract the downloaded file
    echo "Extracting files..."
    tar -xzvf "$OUTPUT_FILE" -C "$EXTRACT_DIR"

    # Check if extraction was successful
    if [ $? -eq 0 ]; then
        echo "Files extracted successfully to $EXTRACT_DIR"
    else
        echo "Failed to extract files"
        exit 1
    fi
else
    echo "Failed to download the file"
    exit 1
fi

# Clean up: remove the downloaded tar.gz file
rm "$OUTPUT_FILE"
echo "Cleaned up: removed $OUTPUT_FILE"

echo "Process completed successfully"

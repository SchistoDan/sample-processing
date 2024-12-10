#!/bin/bash

#NCBI Taxonomy Downloader and Extractor
#This script downloads and extracts the NCBI taxonomy database from the NCBI FTP server.
#It handles the download, extraction, and cleanup process automatically, with error
#checking at each step.

#Features:
#- Downloads the latest new_taxdump.tar.gz from NCBI
#- Extracts the contents to a specified directory
#- Performs cleanup of temporary files
#- Includes error checking and status reporting
#- Handles existing files and directories safely
#- Asks for confirmation before deleting directories

#Usage:
#    ./download_taxonomy.sh <output_file> path/to/<extract_dir>

#Arguments:
#    output_file   : Temporary file path for the downloaded archive
#                   (will be deleted after successful extraction)
#    extract_dir   : Directory where the taxonomy files will be extracted
#                   (will be created if it doesn't exist)

#Example:
#    ./download_taxonomy.sh /tmp/taxdump.tar.gz /data/taxonomy


# Function to ask for confirmation before deleting directory
confirm_delete() {
    local dir="$1"
    read -p "Directory '$dir' exists. Do you want to delete it? (y/N): " response
    case "$response" in
        [yY][eE][sS]|[yY]) 
            return 0
            ;;
        *)
            return 1
            ;;
    esac
}

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

# Download the file using wget with specific options to handle encoding
echo "Downloading file..."
wget --no-check-certificate -O "${OUTPUT_FILE}" "$URL"

# Check if the download was successful and if the file exists
if [ $? -eq 0 ] && [ -f "$OUTPUT_FILE" ]; then
    echo "File downloaded successfully as $OUTPUT_FILE"
    echo "File size: $(ls -lh "$OUTPUT_FILE" | awk '{print $5}')"
    
    # Check if extraction directory exists and ask for confirmation before removing
    if [ -d "$EXTRACT_DIR" ]; then
        if ! confirm_delete "$EXTRACT_DIR"; then
            echo "Directory deletion cancelled by user"
            echo "Cleaning up downloaded file..."
            rm "$OUTPUT_FILE"
            exit 1
        fi
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
        echo "Extracted files:"
        ls -lh "$EXTRACT_DIR"
        
        # Only remove the downloaded file after successful extraction
        echo "Cleaning up..."
        rm "$OUTPUT_FILE"
        echo "Cleaned up: removed $OUTPUT_FILE"
        echo "Process completed successfully"
    else
        echo "Failed to extract files from $OUTPUT_FILE"
        echo "Tar command failed. Contents of output file:"
        file "$OUTPUT_FILE"
        exit 1
    fi
else
    echo "Failed to download the file or file not found after download"
    echo "Wget exit code: $?"
    echo "File exists check: [ -f $OUTPUT_FILE ] = $([ -f "$OUTPUT_FILE" ] && echo "true" || echo "false")"
    exit 1
fi

#!/bin/bash

# Clear the output file if it exists for M1
> kt1.M1.pro

# Loop through files 00000 to 00049
for i in {0..49}; do
    # Format number with leading zeros (5 digits)
    file_num=$(printf "%05d" $i)

    # Construct the S3 path
    s3_path="s3://serratus-rayan/beetles/logan_mar26_26_run/all-concat/kt1-diamond-mar26_26-${file_num}.pro.zst"

    echo "Processing ${s3_path}..."

    # Download file, decompress, grep, and append to output file
    aws s3 cp "$s3_path" - | zstd -d | grep "kt1.kt;M1;" >> kt1.M1.pro

    # Check if the download was successful
    if [ $? -ne 0 ]; then
        echo "Warning: Failed to process ${s3_path}"
    fi
done

echo "Done! Results written to kt1.M1.pro"
echo "Total matches found: $(wc -l < kt1.M1.pro)"
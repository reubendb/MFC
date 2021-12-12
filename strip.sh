#!/bin/bash

FILES="./src/pre_process_code/*.f90 ./src/simulation_code/*.f90 ./src/post_process_code/*.f90"

for f in $FILES; do
  echo "Processing $f file..."
  # cat "$f"
    # tail -n +27 "$f" > "$f.tmp" && mv "$f.tmp" "$f"
    sed '4,6d' "$f" > "$f.tmp" && mv "$f.tmp" "$f"
done

#!/usr/bin/env bash
# fetch_cmems.sh – Download a CMEMS subset via the copernicusmarine CLI
#
# Usage:
#   ./fetch_cmems.sh <dataset_id> <variable> <time_start> <time_end> \
#                    <min_lat> <max_lat> <min_lon> <max_lon> <output_file>
#
# Prerequisites:
#   pip install copernicusmarine
#   copernicusmarine login   (free CMEMS credentials – https://data.marine.copernicus.eu/)

set -euo pipefail

if [ "$#" -ne 9 ]; then
  echo "Usage: $0 <dataset_id> <variable> <time_start> <time_end> \\" >&2
  echo "          <min_lat> <max_lat> <min_lon> <max_lon> <output_file>" >&2
  exit 1
fi

dataset_id="$1"
variable="$2"
time_start="$3"
time_end="$4"
min_lat="$5"
max_lat="$6"
min_lon="$7"
max_lon="$8"
output_file="$9"

copernicusmarine subset \
  --dataset-id        "$dataset_id" \
  --variable          "$variable" \
  --start-datetime    "$time_start" \
  --end-datetime      "$time_end" \
  --minimum-latitude  "$min_lat" \
  --maximum-latitude  "$max_lat" \
  --minimum-longitude "$min_lon" \
  --maximum-longitude "$max_lon" \
  --force-download \
  --output-filename   "$output_file"

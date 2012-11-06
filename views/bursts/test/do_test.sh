#!/bin/sh

set -e

rm -rf tmp
mkdir -p tmp
for i in i*.jpg; do
	iion $i tmp/`basename $i .jpg`.png
done

cd tmp


burst_registration.sh i%02d.png reg_i%02d.png 0 5
burst_combination.sh reg_i%02d.png 0 5 denoised.png

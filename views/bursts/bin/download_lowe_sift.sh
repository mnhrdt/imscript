#!/bin/sh

wget http://www.cs.ubc.ca/~lowe/keypoints/siftDemoV4.zip
unzip siftDemoV4.zip
mv siftDemoV4/sift sift_executable_lowe
rm -rf siftDemoV4 siftDemoV4.zip

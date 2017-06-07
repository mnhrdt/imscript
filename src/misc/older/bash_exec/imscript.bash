#!/bin/bash

binpath="/path2bin"
exec="$binpath/$@"
$exec 0<&0

#!/usr/bin/env bash

export input1=$1

grep -v "GL\|KI\|Un" $input1 > tmp && mv tmp $input1


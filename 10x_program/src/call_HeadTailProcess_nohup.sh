#!/bin/bash
nohup python $1 -s $2 -f $3 -g $4 > $2.log 2> $2.err &
#!/bin/bash

for i in {0.02,0.04,0.05,0.06,0.08,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}
do
    python treeTrials.py -d -s $i -t 7000 -w
done

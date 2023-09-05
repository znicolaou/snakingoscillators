#!/usr/bin/bash

i=8
br='rv5'
filebase='data/janus'
mkdir -p $filebase/$i/invariant
cp {janus_rel.c,c.forward,invariant.auto} $filebase/$i/invariant
cp $filebase/$i/s.start_$br $filebase/$i/invariant/s.start
cp $filebase/$i/b.start_$br $filebase/$i/invariant/b.start
cp $filebase/$i/d.start_$br $filebase/$i/invariant/d.start
cd $filebase/$i/invariant && auto invariant.auto

#!/bin/bash
if [[ "$FC" == "pgf90" ]]
          then
              wget -q -O /dev/stdout \
              'https://raw.githubusercontent.com/nemequ/pgi-travis/master/install-pgi.sh' | \
              /bin/sh
fi
#!/bin/bash

./build/latticecut params/bcc.config 2> >(sed $'s/.*/\033[31m&\033[0m/' >&2)
./build/latticecut params/fcc.config 2> >(sed $'s/.*/\033[31m&\033[0m/' >&2)
./build/latticecut params/hc.config 2> >(sed $'s/.*/\033[31m&\033[0m/' >&2)
./build/latticecut params/sc.config 2> >(sed $'s/.*/\033[31m&\033[0m/' >&2)
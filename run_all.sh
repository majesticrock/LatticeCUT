#!/bin/bash

./build/default/latticecut params/bcc.config 2> >(sed $'s/.*/\033[31m&\033[0m/' >&2)
./build/default/latticecut params/fcc.config 2> >(sed $'s/.*/\033[31m&\033[0m/' >&2)
./build/default/latticecut params/hc.config 2> >(sed $'s/.*/\033[31m&\033[0m/' >&2)
./build/default/latticecut params/sc.config 2> >(sed $'s/.*/\033[31m&\033[0m/' >&2)
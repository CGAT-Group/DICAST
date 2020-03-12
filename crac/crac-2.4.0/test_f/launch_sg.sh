#!/bin/bash

code=0
for file in $(find . -name "*.should_get"); do
    ./should-get.sh $file
    if [ $? -ne 0 ]; then
        code=1
    fi
done

exit $code
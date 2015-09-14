#!/usr/bin/env bash

# Setting the PYTHONPATH
export PYTHONPATH=..:$PYTHONPATH

# Creating the HTML
make html latexpdf doctest

exit 0

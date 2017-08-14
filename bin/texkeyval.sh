#!/usr/bin/env bash

# modify file so it can be used as a latex key val pair

echo "\newcommand{\\$2}[1]{%"
echo "\IfEqCase{#1}{%"
cat $1
echo "}[\PackageError{tree}{Undefined option to tree: #1}{}]%"
echo "}%"



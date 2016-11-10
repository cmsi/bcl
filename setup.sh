#!/bin/sh

wget -O - https://github.com/todo-group/standards/archive/develop.tar.gz | tar zxf - standards-develop/config && mv standards-develop/config/* config/ && rm -rf standards-develop

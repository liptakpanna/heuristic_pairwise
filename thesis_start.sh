#!/bin/bash

monetdbd start dbfarm/
monetdb destroy thesis
monetdb create thesis
monetdb release thesis
monetdbd stop dbfarm/

mserver5 --dbpath=dbfarm/thesis --set embedded_py=true
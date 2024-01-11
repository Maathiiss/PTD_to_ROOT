#!/bin/bash

source /sps/nemo/sw/snswmgr/snswmgr.conf
snswmgr_load_setup falaise@5.1.0

export ROOT_INCLUDE_PATH=`bxquery --prefix`/include/bayeux:$ROOT_INCLUDE_PATH


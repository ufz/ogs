#!/usr/bin/env bash

set -ex
guix time-machine -C scripts/guix/channels.scm -- \
  shell --container --nesting --network --development $1
  openssl nss-certs coreutils bash git

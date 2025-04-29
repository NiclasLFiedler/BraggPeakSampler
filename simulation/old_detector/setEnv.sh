#!/bin/bash
read -p "How many detectors are used: " nbofdetectors
read -p "How many detectors are PWO: " nbofdetectorsPWO

export ENV_NBDETECTORS="$nbofdetectors"
export ENV_NBDETECTORSPWO="$nbofdetectorsPWO"
#!/usr/bin/env bash
set -e

version=`exec ls -v dist | tail -1 | sed 's/oncodriveclustl-//g' | sed 's/\.tar\.gz//g'`
singularity_file=$(mktemp)

if [ -z "${genome}" ]
then
    echo "genome is unset";
    image="oncodriveclustl_${version}.simg"
    sed "s@PKG_VERSION@${version}@g" Singularity > ${singularity_file}
else
    echo "genome is set to '${genome}'";
    image="oncodriveclustl_${version}_${genome}.simg"
    cat Singularity.genome | sed "s@PKG_VERSION@${version}@g" | sed "s@GENOME@${genome}@g"  > ${singularity_file}
fi

sudo singularity build ${image} ${singularity_file}
sudo chown ${USER}: ${image}

Bootstrap: docker
From: debian:buster-slim

%environment
    export LC_ALL=C.UTF-8

%runscript
    exec "/usr/local/bin/oncodriveclustl" "$@"

%files
    dist/oncodriveclustl-PKG_VERSION.tar.gz /oncodriveclustl.tar.gz

%post

    # Load environment
    export LC_ALL=C.UTF-8

    # Install dependencies
    apt update
    apt install -y python3 python3-pip python3-sklearn python3-statsmodels \
                   python3-matplotlib python3-scipy python3-pandas python3-numpy python3-pycurl \
                   python3-click python3-requests python3-tqdm python3-daiquiri python3-humanize \
                   python3-certifi python3-configobj python3-sortedcontainers

    # Install OncodriveCLUSTL
    pip3 --no-cache-dir install /oncodriveclustl.tar.gz

    # Clean unused things
    apt remove -y python3-pip
    apt autoremove -y
    rm -rf /var/lib/apt/lists/*

%test
    /usr/local/bin/oncodriveclustl --help


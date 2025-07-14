# hash:sha256:b38819aeda263b5d460c299dfb2b7bf9bdee8c432b87cccdad0c18596b762482
FROM registry.codeocean.com/codeocean/matlab:2020b-ubuntu20.04

ARG DEBIAN_FRONTEND=noninteractive
ARG MLM_LICENSE_FILE

RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys \
        E298A3A825C0D65DFD57CBB651716619E084DAB9 \
    && apt-get update \
    && apt-get install -y --no-install-recommends software-properties-common \
    && add-apt-repository -y 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' \
    && add-apt-repository -y 'deb https://mirror.linux.ncsu.edu/ubuntu/ focal-backports main restricted universe' \
    && apt-get purge -y --autoremove software-properties-common \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
        appstream=0.12.10-2 \
        build-essential=12.8ubuntu1.1 \
        cmake=3.16.3-1ubuntu1.20.04.1 \
        curl=7.68.0-1ubuntu2.24 \
        dmsetup=2:1.02.167-1ubuntu1 \
        elpa-ess=18.10.2-2focal0 \
        gcc=4:9.3.0-1ubuntu2 \
        gfortran=4:9.3.0-1ubuntu2 \
        gstreamer1.0-tools=1.16.3-0ubuntu1.1 \
        isoquery=3.2.3-1 \
        libblas-dev=3.9.0-1build1 \
        libcairo2-dev=1.16.0-4ubuntu1 \
        libcurl4-doc=7.68.0-1ubuntu2.24 \
        libcurl4-openssl-dev=7.68.0-1ubuntu2.24 \
        libfile-mimeinfo-perl=0.29-1 \
        libfontconfig1-dev=2.13.1-2ubuntu3 \
        libfreetype6-dev=2.10.1-2ubuntu0.3 \
        libfribidi-dev=1.0.8-2ubuntu0.1 \
        libgcrypt20-doc=1.8.5-5ubuntu1.1 \
        libgit2-dev=0.28.4+dfsg.1-2ubuntu0.1 \
        libharfbuzz-dev=2.6.4-1ubuntu4.2 \
        libidn11-dev=1.33-2.2ubuntu2 \
        libjpeg-dev=8c-2ubuntu8 \
        libkrb5-dev=1.17-6ubuntu4.7 \
        liblapack-dev=3.9.0-1build1 \
        libldap2-dev=2.4.49+dfsg-2ubuntu1.10 \
        libmbedtls-doc=2.16.4-1ubuntu2 \
        libnet-dbus-perl=1.2.0-1 \
        libnss-systemd=245.4-4ubuntu3.24 \
        libpam-cap=1:2.32-1ubuntu0.1 \
        libpng-dev=1.6.37-2 \
        librtmp-dev=2.4+20151223.gitfa8646d.1-2build1 \
        libssl-dev=1.1.1f-1ubuntu2.23 \
        libssl-doc=1.1.1f-1ubuntu2.23 \
        libtiff5-dev=4.1.0+git191117-2ubuntu0.20.04.14 \
        libx11-protocol-perl=0.56-7 \
        libxml2-dev=2.9.10+dfsg-5ubuntu0.20.04.7 \
        networkd-dispatcher=2.1-2~ubuntu20.04.3 \
        packagekit-tools=1.1.13-2ubuntu1.1 \
        pandoc=2.5-3build2 \
        pkg-config=0.29.1-0ubuntu4 \
        python-apt-doc=2.0.1ubuntu0.20.04.1 \
        python-dbus-doc=1.2.16-1build1 \
        python3-apt-dbg=2.0.1ubuntu0.20.04.1 \
        python3-cryptography=2.8-3ubuntu0.3 \
        python3-dbus-dbg=1.2.16-1build1 \
        python3-openssl=19.0.0-1build1 \
        python3-socks=1.6.8+dfsg-1ubuntu1 \
        r-base=4.4.1-1.2004.0 \
        r-base-core=4.4.1-1.2004.0 \
        r-base-html=4.4.1-1.2004.0 \
        r-doc-html=4.4.1-1.2004.0 \
        r-doc-info=4.4.1-1.2004.0 \
        r-doc-pdf=4.4.1-1.2004.0 \
        r-mathlib=4.4.1-1.2004.0 \
        r-recommended=4.4.1-1.2004.0 \
        rng-tools=5-1ubuntu2 \
        shared-mime-info=1.15-1 \
        systemd-container=245.4-4ubuntu3.24 \
        tcl8.6=8.6.10+dfsg-1 \
        tk8.6=8.6.10-1 \
        unattended-upgrades=2.3ubuntu0.3 \
        x11-utils=7.7+5 \
        x11-xserver-utils=7.7+8 \
        xdg-user-dirs=0.17-2ubuntu1 \
        xfonts-base=1:1.0.5 \
        xvfb=2:1.20.13-1ubuntu1~20.04.17 \
    && rm -rf /var/lib/apt/lists/*

RUN echo 'options(repos = c(CRAN = "https://cloud.r-project.org/"), download.file.method = "libcurl")' >> $(Rscript -e 'cat(R.home())')/etc/Rprofile.site \
    && echo 'options(Ncpus = parallel::detectCores())' >> $(Rscript -e 'cat(R.home())')/etc/Rprofile.site \
    && Rscript -e 'options(warn=2); install.packages("remotes")'
RUN Rscript -e 'remotes::install_version("BiocManager", "1.30.25")' \
    && Rscript -e 'remotes::install_version("MatrixModels", "0.5-3")' \
    && Rscript -e 'remotes::install_version("Rcpp", "1.0.13")' \
    && Rscript -e 'remotes::install_version("RcppEigen", "0.3.4.0.2")' \
    && Rscript -e 'remotes::install_version("carData", "3.0-5")' \
    && Rscript -e 'remotes::install_version("caret", "6.0-94")' \
    && Rscript -e 'remotes::install_version("conflicted", "1.2.0")' \
    && Rscript -e 'remotes::install_version("corrplot", "0.95")' \
    && Rscript -e 'remotes::install_version("data.table", "1.16.2")' \
    && Rscript -e 'remotes::install_version("dbplyr", "2.5.0")' \
    && Rscript -e 'remotes::install_version("devEMF", "4.5")' \
    && Rscript -e 'remotes::install_version("doBy", "4.6.24")' \
    && Rscript -e 'remotes::install_version("drda", "2.0.4")' \
    && Rscript -e 'remotes::install_version("dtplyr", "1.3.1")' \
    && Rscript -e 'remotes::install_version("export", "0.3.0")' \
    && Rscript -e 'remotes::install_version("fontquiver", "0.2.1")' \
    && Rscript -e 'remotes::install_version("ggpubr", "0.6.0")' \
    && Rscript -e 'remotes::install_version("ggrepel", "0.9.5")' \
    && Rscript -e 'remotes::install_version("ggsci", "3.2.0")' \
    && Rscript -e 'remotes::install_version("ggsignif", "0.6.4")' \
    && Rscript -e 'remotes::install_version("googlesheets4", "1.1.1")' \
    && Rscript -e 'remotes::install_version("gridExtra", "2.3")' \
    && Rscript -e 'remotes::install_version("haven", "2.5.4")' \
    && Rscript -e 'remotes::install_version("janitor", "2.2.0")' \
    && Rscript -e 'remotes::install_version("knitr", "1.48")' \
    && Rscript -e 'remotes::install_version("lubridate", "1.9.3")' \
    && Rscript -e 'remotes::install_version("numDeriv", "2016.8-1.1")' \
    && Rscript -e 'remotes::install_version("officer", "0.6.6")' \
    && Rscript -e 'remotes::install_version("openxlsx", "4.2.7")' \
    && Rscript -e 'remotes::install_version("pROC", "1.18.5")' \
    && Rscript -e 'remotes::install_version("polynom", "1.4-1")' \
    && Rscript -e 'remotes::install_version("readxl", "1.4.3")' \
    && Rscript -e 'remotes::install_version("remotes", "2.5.0")' \
    && Rscript -e 'remotes::install_version("reprex", "2.1.1")' \
    && Rscript -e 'remotes::install_version("rgl", "1.3.1")' \
    && Rscript -e 'remotes::install_version("rmarkdown", "2.28")' \
    && Rscript -e 'remotes::install_version("rvest", "1.0.4")' \
    && Rscript -e 'remotes::install_version("scales", "1.3.0")' \
    && Rscript -e 'remotes::install_version("selectr", "0.4-2")' \
    && Rscript -e 'remotes::install_version("stargazer", "5.2.3")' \
    && Rscript -e 'remotes::install_version("systemfonts", "1.1.0")' \
    && Rscript -e 'remotes::install_version("textshaping", "0.4.0")' \
    && Rscript -e 'remotes::install_version("tictoc", "1.2.1")' \
    && Rscript -e 'remotes::install_version("tidyverse", "2.0.0")' \
    && Rscript -e 'remotes::install_version("xtable", "1.8-4")'

RUN Rscript -e 'options(warn=2); install.packages("BiocManager")'
RUN Rscript -e 'options(warn=2); BiocManager::install(c( \
        "BiocVersion", \
        "cmapR" \
    ))' # Original versions: 3.19.1 1.16.0

COPY postInstall /
RUN /postInstall

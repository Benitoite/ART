FROM debian:testing

#   add the dependencies

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends build-essential locales cmake git curl libcanberra-gtk3-dev libexiv2-dev libexpat-dev libfftw3-dev libglibmm-2.4-dev libgtk-3-dev libgtkmm-3.0-dev libiptcdata0-dev libjpeg-dev liblcms2-dev liblensfun-dev liblensfun-bin liblensfun-data-v1 libpng-dev libsigc++-2.0-dev libtiff5-dev zlib1g-dev librsvg2-dev ca-certificates ssl-cert -y

#   prepare the environment

RUN locale-gen C.UTF-8
ENV LANG C.UTF-8
ENV LANGUAGE C.UTF-8:en
ENV LC_ALL C.UTF-8

#   clone source code, checkout dev branch

RUN mkdir -p ~/programs && git clone http://bitbucket.org/agriggio/ART.git ~/programs/code-ART && cd ~/programs/code-ART && git checkout master

#   update lensfun data
RUN cd ~/programs && curl -o lensfun-update-data https://raw.githubusercontent.com/Benitoite/lensfun/patch-1/apps/lensfun-update-data && chmod +x lensfun-update-data && ./lensfun-update-data

#   configure build system and compile

RUN cd ~/programs/code-ART && mkdir build && cd build && cmake \
    -DCMAKE_BUILD_TYPE="release"  \
    -DCACHE_NAME_SUFFIX="5-dev" \
    -DPROC_TARGET_NUMBER="1" \
    -DBUILD_BUNDLE="ON" \
    -DBUNDLE_BASE_INSTALL_DIR="$HOME/programs/ART" \
    -DOPTION_OMP="ON" \
    -DWITH_LTO="ON" \
    -DWITH_PROF="OFF" \
    -DWITH_SAN="OFF" \
    -DWITH_SYSTEM_KLT="OFF" \
    ..
RUN cd ~/programs/code-ART/build && make -j$(nproc --all) && make install

#   set the entrypoint command

LABEL maintainer="kd6kxr@gmail.com"
CMD echo "This is a test..." && ~/programs/ART/art && echo "THATS ALL FOLKS!!!"

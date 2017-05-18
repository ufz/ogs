wget http://downloads.conan.io/latest_debian -O conan.deb
sudo dpkg -i conan.deb
rm conan.deb

conan install -u -s compiler=gcc -s compiler.version=4.9

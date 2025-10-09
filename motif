#to install
1) Type the following commands and then follow the instructions
   printed by the "configure" command.

    $ tar zxf meme_VERSION.tar.gz
    $ cd meme_VERSION 
    $ ./configure --prefix=$HOME/meme --enable-build-libxml2 --enable-build-libxslt
# Compile (takes 5-10 minutes)
make
 
# Install
make install
 
# Now the bin directory should exist
ls ~/meme/bin/fimo
 
# Add to PATH
echo 'export PATH=$HOME/meme/bin:$PATH' >> ~/.zshrc
source ~/.zshrc
 
# Test
fimo --version

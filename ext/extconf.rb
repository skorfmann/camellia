require 'mkmf'
CAMELLIA_VERSION = "2.8.0"
puts "Unzipping Camellia library..."
`tar xvf ../CamelliaLib-#{CAMELLIA_VERSION}.tar.gz`
Dir.chdir "CamelliaLib-#{CAMELLIA_VERSION}"
puts "Configuring Camellia library..."
`./configure`
puts "Compiling and installing Camellia library..."
`make install`
Dir.chdir ".."
CONFIG['LDSHARED'] = "g++ -shared -lCamellia"
create_makefile('camellia')


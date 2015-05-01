require 'formula'

class Bissnp < Formula
  homepage 'http://sourceforge.net/projects/bissnp'
  version '0.82.2'
  url 'http://sourceforge.net/projects/bissnp/files/BisSNP-0.82.2.jar'
  sha1 'c3c0a48070675fda1d3dd4b7bc239876113cad77'

  def install
    java = share / 'java' / 'bissnp'
    java.install Dir['*.jar']
    bin.mkdir
    open(bin / 'bissnp', 'w') do |file|
      file.write <<-EOS.undent
        #!/bin/bash
        default_jvm_mem_opts="-Xms750m -Xmx2g"
        jvm_mem_opts=""
        jvm_prop_opts=""
        pass_args=""
        for arg in "$@"; do
            case $arg in
                '-D'*)
                    jvm_prop_opts="$jvm_prop_opts $arg"
                    ;;
                 '-Xm'*)
                    jvm_mem_opts="$jvm_mem_opts $arg"
                    ;;
                 '*8')
                    pass_args="$pass_args \\*8"
                    ;;
                 *)
                    pass_args="$pass_args $arg"
                    ;;
            esac
        done
        if [ "$jvm_mem_opts" == "" ]; then
            jvm_mem_opts="$default_jvm_mem_opts"
        fi
        eval java $jvm_mem_opts $jvm_prop_opts -jar #{java}/BisSNP-#{version}.jar $pass_args
      EOS
    end
  end

  test do
    system 'qsignature -h'
  end
end

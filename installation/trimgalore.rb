class Trimgalore < Formula
  homepage "http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/"
  # tag "bioinformatics"

  url "http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/trim_galore_v0.3.7.zip" 
  sha256 "f8adfab475452b9c1e9d0e94ffc7bfa2cd2ae6f6e0e7c0cf5c354d9a8fe66680"

  def install
    chmod 0755, "trim_galore"
    prefix.install Dir["*"]
    mkdir_p bin
    ln_s prefix/"trim_galore", bin/"trim_galore"
  end

  test do
    system "trim_galore", "-h"
  end
end

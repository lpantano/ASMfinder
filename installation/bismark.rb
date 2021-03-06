class Bismark < Formula
  homepage "http://www.bioinformatics.babraham.ac.uk/projects/bismark/"
  # tag "bioinformatics"

  url "http://www.bioinformatics.babraham.ac.uk/projects/bismark/bismark_v0.14.2.tar.gz"
  sha256 "19560979bf3d060dd48079696b595cedd5bd4e8685da95a5024365e071cee58b"

  def install
    bin.install Dir["bismark*"]
  end

  test do
    system "bismark", "-h"
  end
end

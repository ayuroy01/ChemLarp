# Homebrew formula for ChemLarp. To publish: create a tap repository (github.com/ayuroy01/homebrew-chemlarp),
# copy this file to Formula/chemlarp.rb there, and run scripts/release_manifests.py after tagging a release to
# fill in the sha256. Until then `brew install --HEAD chemlarp` builds from the main branch.
class Chemlarp < Formula
  desc "Chemistry engine: identifiers, thermodynamics, kinetics, process and molecular modelling"
  homepage "https://ayuroy01.github.io/ChemLarp/"
  url "https://github.com/ayuroy01/ChemLarp/archive/refs/tags/v@VERSION@.tar.gz"
  sha256 "@SHA256@"
  license "AGPL-3.0-only"
  head "https://github.com/ayuroy01/ChemLarp.git", branch: "main"

  def install
    # Only the CLI and the shared library; the test suite is not needed for an installation.
    system "make", "-j#{ENV.make_jobs}", "BUILD=build/brew", "install", "PREFIX=#{prefix}"
  end

  test do
    assert_equal "c1ccccc1", shell_output("#{bin}/chemlarp canon C1=CC=CC=C1").strip
    assert_match "LFQSCWFLJHTTHZ-UHFFFAOYSA-N", shell_output("#{bin}/chemlarp inchi CCO --key")
  end
end

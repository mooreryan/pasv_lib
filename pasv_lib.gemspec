# coding: utf-8
lib = File.expand_path("../lib", __FILE__)
$LOAD_PATH.unshift(lib) unless $LOAD_PATH.include?(lib)
require "pasv_lib/version"

Gem::Specification.new do |spec|
  spec.name          = "pasv_lib"
  spec.version       = PasvLib::VERSION
  spec.authors       = ["Ryan Moore"]
  spec.email         = ["moorer@udel.edu"]

  spec.summary       = %q{Library code for PASV}
  # spec.description   = %q{TODO: Write a longer description or delete this line.}
  spec.homepage      = "http://www.github.com/mooreryan/pasv_lib"

  spec.files         = `git ls-files -z`.split("\x0").reject do |f|
    f.match(%r{^(test|spec|features)/})
  end
  spec.bindir        = "exe"
  spec.executables   = spec.files.grep(%r{^exe/}) { |f| File.basename(f) }
  spec.require_paths = ["lib"]

  spec.add_development_dependency "bundler", "~> 1.15"
  spec.add_development_dependency "rake", "~> 10.0"
  spec.add_development_dependency "rspec", "~> 3.0"

  spec.add_runtime_dependency "blosum", "~> 0.1.0"
  spec.add_runtime_dependency "parse_fasta", "~> 2.5.2"
end

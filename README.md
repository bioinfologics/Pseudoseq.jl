# <img src="./sticker.svg" width="30%" align="right" /> Pseudoseq

_Fake genomes, fake sequencing, **real** insights._

[![Latest Release](https://img.shields.io/github/release/bioinfologics/Pseudoseq.jl.svg)](https://github.com/bioinfologics/Pseudoseq.jl/releases/latest)
[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/bioinfologics/Pseudoseq.jl/blob/master/LICENSE)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://bioinfologics.github.io/Pseudoseq.jl/stable)
[![Latest](https://img.shields.io/badge/docs-dev-blue.svg)](https://bioinfologics.github.io/Pseudoseq.jl/dev)
[![Pkg Status](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![DOI](https://zenodo.org/badge/160373866.svg)](https://zenodo.org/badge/latestdoi/160373866)


## Description

The `Pseudoseq` package allows you to build arbitrary genomes, and simulate DNA
sequencing experiments.
DNA sequencing experiments are modelled conceptually as a sampling process.


## Install

Pseudoseq is built with [BioJulia](https://biojulia.net), and is designed with
compatibility with the BioJulia ecosystem of tools in mind.
Pseudoseq is made available to install through BioJulia's package registry.

Julia by default only watches the "General" package registry, so before you
start, you should add the BioJulia package registry.

Start a julia terminal, hit the `]` key to enter pkg mode (you should see the
prompt change from `julia> ` to `pkg> `), then enter the following command:

```julia
registry add https://github.com/BioJulia/BioJuliaRegistry.git
```

After you've added the registry, you can install Pseudoseq from the julia REPL.
Press `]` to enter pkg mode again, and enter the following:

```julia
add Pseudoseq
```


## Testing

Pseudoseq is tested against Julia `1.X` on Linux, OS X, and Windows.



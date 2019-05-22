# API: MoleculePool

## Exported functions

### Making a pool of molecules

```@docs
makepool
```

### Molecule pool transformations

#### Fragment

```@docs
fragment(p::Pseudoseq.MoleculePool, meansize::Int)
```

#### Subsampling

```@docs
subsample(p::Pseudoseq.MoleculePool, n::Int)
```

#### Tagging

```@docs
tag(p::Pseudoseq.MoleculePool, ntags::Int)
```

#### Flipping

```@docs
flip(p::Pseudoseq.MoleculePool)
```
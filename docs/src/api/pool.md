```@meta
CurrentModule = Pseudoseq.Sequencing
```
# API: MoleculePool

## Exported functions

### Making a pool of molecules

```@docs
Molecules
```

### Molecule pool transformations

#### Fragment

```@docs
fragment(p::Molecules, meansize::Int)
```

#### Subsampling

```@docs
subsample(p::Molecules, n::Int)
```

#### Tagging

```@docs
tag(p::Molecules, ntags::Int)
```

#### Flipping

```@docs
flip(p::Molecules)
```
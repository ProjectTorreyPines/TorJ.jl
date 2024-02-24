## TorJ -- TOroilda Rays in Julia

#### Install

See https://pkgdocs.julialang.org/v1 to lear about Julia packages system

```
cd ~./julia/dev
git clone git@github.com:ProjectTorreyPines/TorJ.git TorJ
julia
]
dev TorJ
```

NOTE: TorJ uses ModelingToolkits.jl underneath, which can take some time (few minutes) to precompile.

#### Jupyter notebooks

Install jupyter-lab using conda. Then add the Julia kernel to Jupyter-lab with

julia
]
add IJulia WebIO

#### Usage

Open `TorJ/test/torJ.ipynb` in jupyter-lab (or VS-code with the Jupyter extension)

If you are not a fan of Jupyter notebooks, take a look and `TorJ/test/runtests.jl`
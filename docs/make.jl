using MQDT
using Documenter

DocMeta.setdocmeta!(MQDT, :DocTestSetup, :(using MQDT); recursive=true)

makedocs(;
    modules=[MQDT],
    authors="Acme Corp",
    sitename="MQDT.jl",
    format=Documenter.HTML(; canonical="https://pairinteraction.github.io/MQDT.jl", edit_link="main", assets=String[]),
    pages=["Home" => "index.md"],
    # doctest = :fix,  # uncomment to automatically fix doctests ouputs
)

deploydocs(; repo="github.com/pairinteraction/MQDT.jl", devbranch="main")

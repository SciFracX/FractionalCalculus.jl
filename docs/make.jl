using Documenter, FractionalCalculus

makedocs(
    sitename = "FractionalCalculus.jl",
    author = "Qingyu Qu",
    clean = true,
    doctest = false,
    modules = [FractionalCalculus],


    pages=[
        "Home" => "index.md",
    ]


    deploydocs(
        repo = "github.com/ErikQQY/FractionalCalculus.jl.git";
        push_preview = true
    )
)
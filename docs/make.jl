using Documenter, FractionalCalculus

makedocs(
    modules = [FractionalCalculus],
    sitename = "FractionalCalculus.jl",
    format = Documenter.HTML(
        assets = ["assets/favicon.ico"],
    ),
)

deploydocs(
    repo = "github.com/ErikQQY/FractionalCalculus.jl.git";
)

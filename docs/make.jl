using Documenter, FractionalCalculus

makedocs(
    modules = [FractionalCalculus],
    sitename = "FractionalCalculus.jl",
    format = Documenter.HTML(
        assets = ["assets/logo.png"],
    ),
)

deploydocs(
    repo = "github.com/ErikQQY/FractionalCalculus.jl.git";
)

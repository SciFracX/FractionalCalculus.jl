using Documenter, FractionalCalculus

makedocs(
    modules = [FractionalCalculus],
    sitename = "FractionalCalculus.jl",
    format = Documenter.HTML(
        assets = ["assets/favicon.ico"],
    ),

    pages = [
        "index.md",
        "Fractional Derivative" => [
            "Derivative/derivative.md",
            "Derivative/derivativeapi.md"
        ],
        "Fractional Integral" => [
            "Integral/integral.md",
            "Integral/integralapi.md"
        ]
    ]
)

deploydocs(
    repo = "github.com/ErikQQY/FractionalCalculus.jl.git";
)

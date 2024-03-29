using Documenter
using Demo

makedocs(
    sitename = "Demo",
    format = Documenter.HTML(),
    modules = [Demo]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/romanWSgit/Demo.jl.git"
)

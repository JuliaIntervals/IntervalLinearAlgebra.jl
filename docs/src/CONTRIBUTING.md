# IntervalLinearAlgebra.jl contribution guidelines

First of all, huge thanks for your interest in the package! ✨

This page has some hopefully useful guidelines. If this is your first time contributing, please read the [pull request-workflow](#Pull-request-workflow) section, mainly to make sure everything works smoothly and you don't get stuck with some nasty technicalities. 

You are also encouraged to read the coding and documentation guidelines, but you don't need to deeply study and memorize those. Core developers are here to help you. Most importantly, relax and have fun!

The core developers of the package can be found in the `#intervals` channel in the Julia slack or zulip, links to join the platforms can be found [here](https://julialang.org/community/)

## Opening issues

If you spot something strange in the software (something doesn't work or doesn't behave as expected) do not hesitate to open a [bug issue](https://github.com/JuliaIntervals/IntervalLinearAlgebra.jl/issues/new?assignees=&labels=bug&template=bug_report.md&title=%5BBUG%5D).

If have an idea of how to make the package better (a new feature, a new piece of documentation, an idea to improve some existing feature), you can open an [enhancement issue](https://github.com/JuliaIntervals/IntervalLinearAlgebra.jl/issues/new?assignees=&labels=enhancement&template=feature_request.md&title=%5Bfeature+request%5D%3A+). 

In both cases, try to follow the template, but do not worry if you don't know how to fill something. 

If you feel like your issue does not fit any of the above mentioned templates (e.g. you just want to ask something), you can also open a [blank issue](https://github.com/JuliaIntervals/IntervalLinearAlgebra.jl/issues/new).

## Pull request workflow

Pull requests are also warmly welcome. For small fixes/additions, feel free to directly open a PR. For bigger more ambitious PRs, it is preferable to open an issue first to discuss it. As a rule of thumb, every pull request should be as atomic as possible (fix one bug, add one feature, address one issue).

### Setup

!!! note
    This is just one way, you can do differently (e.g. clone your fork and add the original repo as `upstream`). In that case, make sure to use the correct remote names

This is something that needs to be done only once, the first time you start contributing

**1.** From the Julia REPL in package mode (you can enter package mode by typing `]`) do

```julia
pkg> dev IntervalLinearAlgebra
```

this will clone the repository into `.julia/dev/IntervalLinearAlgebra`. When you `dev` the package, Julia will use the code in the `dev` folder instead of the official released one. If you want to go back to use the released version, you can do `free IntervalLinearAlgebra`.

**2.** [Fork the repository](https://github.com/juliainterval/IntervalLinearAlgebra.jl).

**3.** Navigate to `.julia/dev/IntervalLinearAlgebra` where you cloned the original repository before. Now you need to add your fork as remote. This can be done with
```
git remote add $your_remote_name $your_fork_url
```

`your_remote_name` can be whatever you want. `your_fork_url` is the url you would use to clone your fork repository. For example if your github username is `lucaferranti` and you want to call the remote `lucaferranti` then the previous command would be

```
git remote add lucaferranti https://github.com/lucaferranti/IntervalLinearAlgebra.jl.git
```

you can verify that you have the correct remotes with `git remote -v` the output should be similar to

```
lucaferranti  https://github.com/lucaferranti/IntervalLinearAlgebra.jl.git (fetch)
lucaferranti  https://github.com/lucaferranti/IntervalLinearAlgebra.jl.git (push)
origin        https://github.com/JuliaIntervals/IntervalLinearAlgebra.jl.git (fetch)
origin        https://github.com/JuliaIntervals/IntervalLinearAlgebra.jl.git (push)
```

Now everything is set!

### Contribution workflow

**0.** Navigate to `.julia/dev/IntervalLinearAlgebra` and make sure you are on the main branch. You can check with `git branch` and if needed use `git switch main` to switch to the main branch. **The next steps assume you are in the `IntervalLinearAlgebra` folder**.

**1.** Before you start modifying, it's good to make sure that your local main branch is synchronized with the main branch in the package repo. To do so, run
```
git fetch origin
git merge origin/main
```
 Since you should **never** directly modify the main branch locally, this should not cause any conflicts. If you didn't follow the previous setup instructions, you may need to change `origin` with the appropriate remote name.

**2.** Now create a new branch for the new feature you want to develop. If possible, the branch should start with your name/initials and have a short but descriptive name of what you are doing (no strict rules). For example, if I (Luca Ferranti) want to fix the code that computes the eigenvalues of a symmetric matrix, I would call the branch `lf-symmetric-eigvals` or something like that. You can create a new branch and switch to it with

```
git switch -c lf-symmetric-eigvals
```
If you are targetting a specific issue, you can also name the branch after the issue number, e.g. `lf-42`.

**3.** Now let the fun begin! Fix bugs, add the new features, modify the docs, whatever you do, it's gonna be awesome! Check also the [coding guidelines](#Coding-guideline) and [documentation guidelines](#Documentation-guideline). Do not worry if it feels like a lot of rules, the core developers are here to help and guide.

**4.** It is important to run the tests of the package locally, to check that you haven't accidentally broken anything. You can run the tests with

```
julia --project test/runtests.jl
```

If you have changed the documentation, you can build it locally with

```
julia --project=docs docs/make.jl
```

This will build the docs in the `docs/build` folder, you can open `docs/build/index.html` and check that everything looks nice. Check also in the terminal that you don't have error messages (no broken links, doctests pass).

**5.** When you are ready, commit your changes. If example you want to commit src/file1.jl, src/file2.jl
```
git add src/file1.jl src/file2.jl
git commit -m "short description of what you did"
```

You can also add and commit all changes at once with

```
git commit -a -m "short description of what you did"
```

finally you are ready to push to your fork. If your fork remote is called `lucaferranti` and your branch is called `lf-symmetric-eigvals`, do

```
git push -u lucaferranti lf-symmetric-eigvals
```

The `-u` flag sets the upstream, so next time you want to push to the same branch you can just do `git push`.

**6.** Next, go to the [package repository](https://github.com/juliaintervals/IntervalLinearAlgebra.jl), you should see a message inviting you to open a pull request, do it! Make sure you are opening the PR to `origin/main`. Try to fill the blanks in the pull request template, but do not worry if you don't know anything. Also, your work needs not be polished and perfect to open the pull request! You are also very welcome to open it as a draft and request feedback, assistance, etc.

**7.** If nothing happens within 7 working days feel free to ping Luca Ferranti (@lucaferranti) every 1-2 days until you get his attention.

## Coding guideline

* Try to roughly follow the [bluestyle](https://github.com/invenia/BlueStyle) style guideline.
* If you add new functionalities, they should also be tested. Exported functions should also have a docstring.
* The test folder should roughly follow the structure of the src folder. That is if you create `src/file1.jl` there should also be `test/test_file1.jl`. There can be exceptions, the main point being that both `test` and `src` should have a logical structure and should be easy to find the tests for a given function.
* The `runtests.jl` should have only inlcude statements.

### Package version

Generally, if the pull request changes the source code, a new version of the package should be released. This means, that if you change the source code, you should also update the `version` entry in the `Project.toml`. Since the package is below version 1, the version update rules are
  * update minor version for breaking changes, e.g. `0.3.5`  => `0.4.0`
  * update patch version for non breaking changes, e.g. `0.3.5` => `0.3.6`
  * It is perfectly fine that you are not sure how to update the version. Just mention in the PR and you will receive guidance
* The person who merges the PR also register the new version.

### Add dependency
If the function you are adding needs an external package (say `Example.jl`), this should be added as dependency, to do so

1. Go to `IntervalLinearAlgebra.jl` and start a Julia session and activate the current environment with `julia --project`
2. Enter the package mode (press `]`) and add the package you want to add, e.g `]add Example`.
3. You can verify that the package was added by typing `st` while in package mode. You can exit the package mode by pressing backspace
4. Open the `Project.toml` file, your package should now be listed in the `[deps]` section.
5. In the `[compat]` section, specify the compatibility requirements. Packages are listed alphabetically. More details about specifying compatibility can be found [here](https://pkgdocs.julialang.org/v1/compatibility/)
6. In the `IntervalLinearAlgebra.jl` file, add the line `using Example` together with the other using statements, or `import Example: fun1, fun2` if you are planning to extend those functions.

If the dependency is quite heavy and used only by some functionalities, you may consider adding that as optional dependency. To do so,
1. Repeat the steps 1-5 above
2. In the `[deps]` section of `Project.toml` locate the package you want to make an optional dependency and move the corresponding line to `[extras]`, keep alphabetical ordering.
3. Add the dependency name to the `test` entry in the `[targets]` section
4. In the `IntervalLinearAlgebra.jl` file, locate the `__init__` function and add the line
```julia
@require """Example = "7876af07-990d-54b4-ab0e-23690620f79a" include("file.jl")"""
```
where `file.jl` is the file containing the functions needing `Example.jl`. The line `Example = "7876af07-990d-54b4-ab0e-23690620f79a"` is the same in the Project.toml
5. In `file.jl` the first line should be `using .Example` (or `import .Example: fun1, fun2`), note the dot before the package name. Then write the functions in the file normally

## Documentation guideline

* Documentation is written with [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl). Documentation files are in `docs/src`, generally as markdown file.
* If you want to include a Julia code example that is **not** executed in the markdown file, use ````` ```julia ````` blocks, e.g.

````
```julia
a = 1
b = 2
```
````

* Julia code that **is** executed should use ````` ```@example ````` blocks, e.g.

````
```@example
a = 1
b = 2
```
````

* If you want to reuse variables between `@example` blocks, they should be named, for example

````
```@example filename
a = 1
b = 2
```

... some text ...

```@example filename
c = a + b
```
````

* If you want to run a Julia code block but don't want the output to be displayed, add `nothing # hide` as last line of the code block.

* You can plot and include figures as follows

````
```@example
# code for plotting
savefig("figname.png") # hide
```

![][figname.png]
````

* Use single ticks for inline code ``` `A` ``` and double ticks for maths ``` ``A`` ```. For single line equations, use 
* For single-line equations, use ````` ```math ````` blocks, e.g.
````
```math
|A_cx-b_c| \le A_\Delta|x| + b_\Delta,
```
````

* You can refer to functions in the pacakge with ``` [`func_name`](@ref) ```
* You can quote references with `[[REF01]](@ref)`
* If you want to add references, you can use the following template
````
#### [REF01] 

```@raw html
<ul><li>
```
Author(s), [*Paper name in italic*](link_to_pdf_if_available), other infos (publisher, year, etc.)
```@raw html
<li style="list-style: none"><details>
<summary>bibtex</summary>
```
```
INSERT BIBTEX HERE
```
```@raw html
</details></li></ul>
```
---
````

* If the pdf of the paper is freely (and legally!) available online, make the title of the paper a link to it.
* The reference code should be first 3 letters of first author surname + last two digits of year, e.g `[FER87]`. To disambiguate duplicates, use letter, e.g. `[FER87a]`, `[FER87b]`.

### Docstrings

* Each exported function should have a docstring. The docstring should roughly follow the following structure
```julia
"""
    funname(param1, param2[, optional_param])

A short description (1-2 lines) of what the function does

### Input

Lis of inputs. Not needed if clear from description and signature.

### Output

List of outputs. Not needed if clear from description and signature.

### Notes

Anything else which is important.

### Algorithm

What algorithms the function uses, preferably with references.

### Example

At least one example, formatted as julia REPL, of what the function does.
Preferably, as a doctest.
"""
```

* Optional parameters in the function signature go around brackets.
* List of inputs and outputs can be omitted if the function has few parameters and they are already clearly explained by the function signature and description.
* Examples should be [doctests](https://juliadocs.github.io/Documenter.jl/stable/man/doctests/). Exceptions to this can occur if e.g. the function is not deterministic (random initialization) or requires a heavy optional dependency.

Here is an example

````julia
"""
    something(A::Matrix{T}, b::Vector{T}[, tol=1e-10]) where {T<:Interval}

this function computes the somethig product between the interval matrix ``A`` and 
interval vector ``b``.

### Input

`A`   -- interval matrix
`b`   -- interval vector
`tol` -- (optional), tolerance to compute the something product, default 1e-10

### Output

The interval vector representing the something product.

### Notes

If `A` and `b` are real, use the [`somethingelse`](@ref) function instead.

### Algorithm

The function uses the *something sometimes somewhere* algorithm proposed by Someone in [[SOM42]](@ref).

### Example

```jldoctest
julia> A = [1..2 3..4;5..6 7..8]
2×2 Matrix{Interval{Float64}}:
[1, 2]  [3, 4]
[5, 6]  [7, 8]

julia> b = [-2..2, -2..2]
2-element Vector{Interval{Float64}}:
[-2, 2]
[-2, 2]

julia> something(A, b)
2-element Vector{Interval{Float64}}:
[-1, 1]
[-7, 8]
```
"""
````

## Acknowledgments

Here is a list of useful resources from which this guideline was inspired

* [JuliaReach developers docs](https://github.com/JuliaReach/JuliaReachDevDocs)
* [Making a first Julia pull request](https://kshyatt.github.io/post/firstjuliapr/)
* [ColPrac](https://github.com/SciML/ColPrac)
* [Julia contributing guideline](https://github.com/JuliaLang/julia/blob/master/CONTRIBUTING.md)

#=
# This is the title

this is some text
=#

a = 1

b = 2

# here we talk more

c = a + b

# now we plot something

using Plots

f(x) = x^2

plot(f, -1, 1)
savefig("tmpfig.png"); nothing # hide

# ![](tmpfig.png)

using Plots

f(x, a) = begin
    sin(x+a*pi/2)
end
xs = collect(0.1:0.05:2*pi)
as = collect(0.2:0.1:2.0)
x_grid = [x for x = xs for y = as]
a_grid = [y for x = xs for y = as]
plot3d(x_grid, a_grid, f.(x_grid, a_grid), st = :surface)
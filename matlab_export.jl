using DelimitedFiles

function matlab_export(X,Y,Z,U, name::String)

    folder = "matlab_plots/"

    writedlm(folder*name*"_X.txt",X)
    writedlm(folder*name*"_Y.txt",Y)
    writedlm(folder*name*"_Z.txt",Z)

    # part = real
    writedlm(folder*name*"_ReUxy.txt",real.(U.XY))
    writedlm(folder*name*"_ReUxz.txt",real.(U.XZ))
    writedlm(folder*name*"_ReUyz.txt",real.(U.YZ))

    # part = real
    writedlm(folder*name*"_ImUxy.txt",imag.(U.XY))
    writedlm(folder*name*"_ImUxz.txt",imag.(U.XZ))
    writedlm(folder*name*"_ImUyz.txt",imag.(U.YZ))

end
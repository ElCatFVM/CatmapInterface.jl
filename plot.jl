using Pkg
Pkg.activate(@__DIR__)

using CairoMakie

    set_theme!(Theme(fontsize = 16, Axis = (titlesize =20,
                                            xlabelsize= 10,
                                            ylabelsize = 30, 
                                            xticklabelsize =15,
                                            yticklabelsize = 25,
                                           )
                    ))
    f = Figure()
    tbl1 = (cat = [1,2,3,4,5,6,7,8,9,10,11],
	       height= [1,2,3,4,5,6,7,8,9,10,11])
	
    barplot(f[1,1], tbl1.cat, tbl1.height,
	        axis = (xticks = (1:11, ["G1", "G2", "G3", "G4",
             "G5", "G6", "G7",
             "β1", "β2", "β3", "β4"]),
	                title = "First Order Indices"),
	        )

    save("figure.svg", f)

# Installing Eirene
If you haven't done so already, install Eirene with these steps: 0) ensure you are running the most recent version of Eirene 1) open a Julia REPL, 2) press `]` to switch to the Pkg REPL, and 3) type `add Eirene` and press return.

```
(v1.3) pkg> add Eirene
```
It's fairly common to get error messages when you do this.  Generally they will want you to install additional (supporting) packages.  Read these messages, and follow the instructions.  When in doubt, always double check that you are running the latest version of Eirene (it's been observed that installers such as Homebrew install dated versions of the language, so the only way to be sure you're up to date is to visit the Julia homepage), and use

```
(v1.3) pkg> update
```

to ensure that you are using the latest versions of the supporting packages.

You can now run any function in the Eirene package in the REPL.  If you want to run a fictional function called `fictionalfunction` that takes zero input arguments, you can enter

```
julia> Eirene.fictionalfunction()
```

Typing the prefix `Eirene.` quickly becomes tiresome. Entering command

```
julia> using Eirene
```

will allow you to call Eirene functions without the prefix.

# 1: The Noisy Circle

Let's load a point cloud sampled from a noisy circle, and compute its persistent homology.  The cloud we'll use is stored in a .csv file within the Eirene library, and its file path can be found via

```
julia> filepath = eirenefilepath("noisycircle")
```

The contents of this file can be loaded via

```
julia> using DelimitedFiles
julia> pointcloud = readdlm(filepath,',')
```

This stores the points in the cloud as the columns of matrix `pointcloud`.  To see this cloud for yourself, call

```
julia> ezplot_pjs(pointcloud)
```

![](/Users/gh10/a/c/j/eirene/code/Eirene.jl/images/sample_cloud_circle.png)

To compute persistent homology in dimensions 0 and 1, run

```
julia> C = eirene(pointcloud, model = "pc")
```

Anything you'd like to know about the persistent homology of the cloud can now be accessed via `C`.  For example, if you wanted to plot the persistence diagram or barcode in dimension 1, you could call

```
julia> plotbarcode_pjs(C,dim=1)

julia> plotpersistencediagram_pjs(C,dim=1)
```


![](/Users/gh10/a/c/j/eirene/code/Eirene.jl/images/sample_persistence_diagram.png)

If you hover the mouse over a point in the persistence diagram you'll see some information about the associated persistent homology class:

* `(x coordinate, y coordinate)`: birth and death times
* `class`: a unique id number for the class in question
* `size`: if a value has been given for `size`, then Eirene has computed a cycle representative for the homology class, and `size` refers to the number of cells in the support of that representative.

Suppose you are looking at the persistence diagram in dimension 1, and want to
know more about the unique point sitting high above the diagonal.  Hovering your
mouse over this point, you see that it comes from a persistent homology class
with unique id number 50.  The class was born at time t = 0.141 and died at time
t = 1.48.  Eirene computed a cycle representative for this class, and that
representative involved 153 cells (wow!).  To plot this class, call

```
julia> plotclassrep_pjs(C,class=50,dim=1)
```

![](/Users/gh10/a/c/j/eirene/code/Eirene.jl/images/sample_cycle_circle.png)

To plot the class without point labels, call

```
julia> plotclassrep_pjs(C,class=50,dim=1,showlabels=false)
```

![](/Users/gh10/a/c/j/eirene/code/Eirene.jl/images/sample_cycle_circle_nolabels.png)

and to plot only the points incident to the class representative, call

```
julia> plotclassrep_pjs(C,class=50,dim=1,showlabels=false,showcloud = false)
```

# 2: Noisy Torus

Let's try a larger point cloud.  Run

```
julia> filepath = eirenefilepath("noisytorus")

julia> pointcloud = readcsv(filepath)

julia> ezplot_pjs(pointcloud)

```

and checkout the new point cloud (it's a noisy torus).  Note that

```
julia> size(pointcloud)
(3, 1800)
```

so there are 1800 points in this cloud.  That's large for an H2 compuation, so
let's use a cutoff radius:

```
julia> C = eirene(pointcloud,maxdim=2,maxrad=0.3,model="pc")
```

To find and plot the two primary H1 generators for this cloud:

```
julia> plotbarcode_pjs(C,dim=0:1)

julia> plotpersistencediagram_pjs(C,dim=1)

julia> plotclassrep_pjs(C,class=768,dim=1)

julia> plotclassrep_pjs(C,class=769,dim=1)
```

For H2, use the following. Make sure to set showlabels=false, otherwise the
graphics will be expensive to render!

```
julia> plotbarcode_pjs(C,dim=2)

julia> plotpersistencediagram_pjs(C,dim=2)

julia> plotclassrep_pjs(C,class=12,dim=2,showlabels=false)
```

# 3: Clouds

You can generate several other point clouds and distance matrices via built-in
functions with Eirene.  Check out

```
julia> noisycircle()

julia> noisycircle3()

julia> torus(m = 100,n = 50,mrad=1,nrad = 0.2)

julia> noisytorus(m = 100,n=50,mrad=1,nrad = 0.2,noiserad= 0.5*nrad)

julia> sphere()

julia> matchingcomplex_symmat(2,3)

julia> chessboardcomplex_symmat(numrows=3,numcols=5)

julia> plane2torus(rand(1000,2)*pi/2)

julia> zerodrandmat(10)
```

# 4: Around the World

Suppose you want to explore world geography vis-a-vis networks of neighboring cities.  A wealth of data is available online, and for this exercise we will use a catalog of 7322 cities from simplemaps.com.  This data comes preloaded in the form of a .csv file with the Eirene library, and can be accessed via

```
julia> filepath = eirenefilepath("simplecity")
```

```
julia> a = ezread(filepath)
7323x9 Array{Any,2}:
 "city"           "city_ascii"       "lat"    "lng"  …  "country"      "iso2"  "iso3"  "province"
 "Qal eh-ye Now"  "Qal eh-ye"      34.983   63.1333     "Afghanistan"  "AF"    "AFG"   "Badghis"
 "Chaghcharan"    "Chaghcharan"    34.5167  65.25       "Afghanistan"  "AF"    "AFG"   "Ghor"
 "Lashkar Gah"    "Lashkar Gah"    31.583   64.36       "Afghanistan"  "AF"    "AFG"   "Hilmand"
 "Zaranj"         "Zaranj"         31.112   61.887      "Afghanistan"  "AF"    "AFG"   "Nimroz"
 "Tarin Kowt"     "Tarin Kowt"     32.6333  65.8667  …  "Afghanistan"  "AF"    "AFG"   "Uruzgan"
 "Zareh Sharan"   "Zareh Sharan"   32.85    68.4167     "Afghanistan"  "AF"    "AFG"   "Paktika"
 ⋮                                                   ⋱  ⋮
 "Gweru"          "Gweru"         -19.45    29.82       "Zimbabwe"     "ZW"    "ZWE"   "Midlands"
 "Mutare"         "Mutare"        -18.97    32.65       "Zimbabwe"     "ZW"    "ZWE"   "Manicaland"
 "Kadoma"         "Kadoma"        -18.33    29.9099     "Zimbabwe"     "ZW"    "ZWE"   "Mashonaland West"
 "Chitungwiza"    "Chitungwiza"   -18.0     31.1     …  "Zimbabwe"     "ZW"    "ZWE"   "Harare"
 "Harare"         "Harare"        -17.8178  31.0447     "Zimbabwe"     "ZW"    "ZWE"   "Harare"
 "Bulawayo"       "Bulawayo"      -20.17    28.58       "Zimbabwe"     "ZW"    "ZWE"   "Bulawayo"
 ```

A comment on the meaning of `Array{Any,2}`, which appears in the second line: Julia differentiates arrays based on the type (or class) of elements they can contain.  An n-dimensional array that is formally declared to have elements of type T is called an array of type `Array{T,n}`. A 5x5 matrix with 64-bit floating point entries can be realized either as an array of type `Array{Float64,2}` or, as the name suggests, as an array of type `Array{Any,2}`. Many functions in Julia are specialized to work with numerical arrays, so it's in our interest to extract the numeric part we are interested in, namely columns 3 and 4, sans headers, and convert to an array of type `Array{Float64,2}` before proceeding.

```
julia> b = a[2:end,3:4]
7322x2 Array{Any,2}:
  34.983   63.1333
  34.5167  65.25
  31.583   64.36
   ⋮
 -17.8178  31.0447
 -20.17    28.58
```

```
julia> b = convert(Array{Float64,2},b)
7322x2 Array{Float64,2}:
  34.983   63.1333
  34.5167  65.25
  31.583   64.36
   ⋮
 -17.8178  31.0447
 -20.17    28.58
 ```

Eirene has a built-in function to convert spherical coordinates (in degrees, with fixed radius 1) to 3D Euclidean coordinates.  The `rowsare` keyword argument determines wether rows are treated as points or dimensions.

```
julia> c = latlon2euc(b,model = "points")
3x7322 Array{Float64,2}:
 0.370265  0.344959  0.368622  0.403432  0.344318  0.309031  …   0.796253   0.822829   0.814358   0.81567    0.824296
 0.730885  0.748275  0.767998  0.755149  0.768533  0.781189      0.510205   0.473338   0.491252   0.490971   0.449048
 0.573333  0.566646  0.523733  0.516713  0.53926   0.542442     -0.325073  -0.31449   -0.309017  -0.305991  -0.344807
 ```

Now that we have a point cloud we can take a preliminary look at its shape and scale with the `ezplot_pjs` wrapper. Note you can always select "deny" when PlotlyJS asks for web access, without penalty.

```
julia> ezplot_pjs(c)
```

A cursory inspection shows that a number of interesting features appear at or below epsilon = 0.15, so we'll use that as our initial cutoff. To pass city names to the `eirene` function, extract column 2 of array a, sans header, and store in a 1-dimensional array of type `Array{Any}`. Due to potential errors resulting from nonstandard character strings, it's good practice to use the built-in label sanitzer ezlabel to clean this column before assigning to d.  The wrapper will replace any element of the column that cannot be expressed as an ACIIString with the number corresponding to its row.  Aside: n can be omitted from the expression `Array{T,n}` when n=1.

```
julia> d = ezlabel(a[2:end,2])
7322-element Array{Any,1}:
 "Qal eh-ye"
 "Chaghcharan"
 "Lashkar Gah"
 ⋮
 "Harare"
 "Bulawayo"
 ```

With our inputs in order, it's time to call eirene.  The calculation should take around 2.5GB of memory and 2 min to complete.

```
julia> C = eirene(c,model="pc",maxrad = 0.15,pointlabels=d)
elapsed time: 118.869523387 seconds
Dict{ASCIIString,Any} with 14 entries:
  "symmat"               => 7322x7322 Array{Int64,2}:…
  "filtration"           => Any[[515055,515055,515055,515055,515055,515055,515055,515055,515055,515055  …  515055,51505…
  "lowlab0"              => Any[Int64[],[1,2,3,4,5,6,7,8,9,10  …  7313,7314,7315,7316,7317,7318,7319,7320,7321,7322],[3…
  "firstv"               => Any[[1,2,3,4,5,6,7,8,9,10  …  7314,7315,7316,7317,7318,7319,7320,7321,7322,7323],[1,346,688…
  "filtrationtranslator" => [0.149999,0.149999,0.149999,0.149999,0.149999,0.149999,0.149999,0.149998,0.149998,0.149998 …
  ⋮                      => ⋮
  ```

The package JLD can be used to save this output.  Run `Pkg.add("JLD")` if you have not done so already, and enter

```
julia> JLD.save("demography.jld","C",C)
```

The file will be saved to your current working directory.  To recover what you have saved, run

```
julia> X = JLD.load("demography.jld")
```

This returns a dictionary object `X` with one key, ``"C"``.  Running

```
julia> C = X["C"]
```

will leave `C` unchanged.

We did not specify `bettimax`, so only the persistence modules in dimensions 0 and 1 are computed.  To view the persistence diagram in dimension 1, enter

```
julia> plotpersistencediagram_pjs(C)
```

Hovering over a point in the diagram will display the identification number of the corresponding persistent homology class, together with the size the cycle representative Eirene computed for it. Fix a feature of interest - say number 1757 - and plot the corresponding cycle, noting that only vertices will appear in the figure (cells of dimension 1 and higher are generally too numerous to plot efficiently).

Note: The cost of displaying text labels for individual points in the Plotly platform is higher than that of displaying points alone.  Consider using `showlabels=false`, if your graphics card has any difficulties.

```
julia> plotclassrep_pjs(C,class = 1757)
```

Class 1757, shown above, is among my favorites: the Himalayan branch of the silk road is clearly visible on its South West arc; to the North it follows the Trans-Siberian railroad from Moscow in the west to Vladivostok on the Sea of Japan, by way of Omsk, Irktsk, and Chita; from there, it follows the connecting route from Beijing to Hong Kong, passing Nanning, Hanoi, and close to Ho Chi Min on its way to Bangkok (not all of these appear in the cycle itself, but they are easily spotted nearby).  With a little exploration you can find large features formed by the Sahara Desert, Tapajos River, Falkland Islands, South China Sea, Guam, and Hudson Bay, and smaller ones shaped by the Andes mountains and Gulf of Mexico - all through the proxy of urban development.

To plot the same class-representative without text labels, use

```
julia> plotclassrep_pjs(C,class = 1757,showlabels=false)
```

To plot only the cities incident to the cycle representative, use

```
julia> plotclassrep_pjs(C,class = 1757,showcloud=false)
```

To plot these using MDS, use

```
julia> plotclassrep_pjs(C,class = 1757,showcloud=false,coords="mds",embeddingobj="hop")
```

To make any further changes to the appearance of your plot, share it with others, or export to still-frames, you can upload to the interactive, open-source Plotly web API.  This requires the setup of a personal account beforehand (see the github documentation for Plotly.jl for instructions), but once this has been done and the connection between your Julia REPL and Plotly account has been established,  uploading is quite simple: create an object p for the PlotlyJS figure, and call Plotly.post

```
julia> p = plotclassrep_pjs(C,class = 1757)

julia> Plotly.post(p)
Plotly.RemotePlot(URI(https://plot.ly/~henselmonster/78))
```

The figure will appear under the My Charts tab of your plotly library.

We hope to add more examples like this one in the future, but there is much left to learn from the Simplemaps data set.  It is instructive, for instance, to observe which features change, and how, under subsampling.  If you have a fun, interesting, or educational example you'd like to share with others, please let us know!

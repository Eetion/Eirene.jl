# Documentation

## One example says it all (almost)

Installing Julia is easy.  Download the latest version of the
[Julia Language](https://julialang.org/downloads/), open a
Julia REPL, and enter

```
julia> Pkg.clone("https://github.com/Eetion/Eirene.jl.git")
julia> using Eirene
```

It's fairly common to get error messages when you do this.  Generally they will want you to install additional (supporting) packages.  Read these messages, and follow the instructions.  When in doubt, *FIRST* check that you are running the latest version of Eirene, and use

```
julia> Pkg.update()
```

to ensure that you are using the latest versions of the supporting packages.  *SECOND*, close the Julia shell.  Some updates will only take effect in new windows.

Now let's plot the 1d persistence diagram for 50 points sampled from the uniform distribution on R<sup>20</sup>.  The first line below generates 20x50 iid matrix and stores it in a variable `x`.  The second asks Eirene to analyze `x` and store the results in a variable `C`.  The third plots the diagram.  Most of what you can do with Eirene is done with minor modifications to these commands.

```
julia> x = rand(20,50)
julia> C = eirene(x, model = "pc")
julia> plotpersistencediagram_pjs(C,dim=1)
```
**Note 1**  Check out `EXAMPLES.md` for more examples!

**Note 2**  If you encounter error messages like the following

```
ERROR: UndefVarError: eirene not defined
ERROR: UndefVarError: plotpersistencediagram_pjs not defined
```

then you will want to append `Eirene.` to the beginning of each function in the Eirene package.  For example, in pace of line 2, you would enter

```
julia> C = Eirene.eirene(x, model = "pc")
```

### Keywords
<br>
  In the example above, the expression `model = "pc"` that sits inside `C = eirene(x, model = "pc")` declares that the columns of `x` should be treated as points in a Euclidean point cloud. We call `model` a *keyword argument*, and `"pc"` its *value*.  Every keyword comes with a default value; you only have to declare a value if you want something other than a default.

`model = "pc", "vr", "complex"`   
* States format of the data.  
* Possible values: `"pc"` (point cloud), `"vr"` (vietoris-rips), `"complex"`.   Default value: `"vr"`.
* See 'Complexes' below for special

`maxdim = k`
* Compute persistent homology in dimensions 0, ..., k.
* Default value: `1`.

`minrad = t`
* Compute homology from time `t` onward.
* Default value: `-Inf`.

`maxrad = t`
* Stop computing homology after time `t`.
* Devalue value `Inf`.

`numrad = N`
* Divide the interval from `minrad` to `maxrad` into `N` equally spaced steps, and compute the homology of each step.  If the value of numrad is set to `Inf`, then homology will computed at every time point.
* Possible values: any positive integer, or `Inf`.  Default value: `Inf`.
* This keyword argument is currently only available for `model = "vr"`.

`record = "cyclerep", "none"`
* Determines wether or not to compute generators.  
* Possible values: `"cyclerep"`, `"none"`.  Default value: `"cyclerep"`.


### Inputs

Eirene computes persistent homology for three types of inputs: distance matrices, point clouds, and complexes.  Here are some formats these can come in, and how to analyze them.

##### Distance matrices

```
julia> C = eirene(x, <keyword arguments>)
```
Formats for `x`
* a symmetric matrix in Julia
* a file path to a symmetric matrix, recorded in a comma or space-delimited text file (.csv or .txt)

##### Point clouds

```
julia> C = eirene(x, model = "pc", <keyword arguments>)
```
Formats for `x`
* a numeric matrix in Julia (columns will be treated as points in Euclidean space)
* a file path to a numeric matrix, recorded in a comma or space-delimited text file (.csv or .txt)

##### Complexes

There are several ways to format the data of a cell complex `E`.  Here are the main ingredients.

* number the cells `1, ..., N` in ascending order, according to dimension
* let `D` be the `N x N` zero/one matrix with `D[i,j] = 1` iff `i` is a face of cell `j`
* store `D` as a sparse matrix in Julia (see Julia docs), and define
```
julia> rv = D.rowval
julia> cp = D.colptr
```

Define vectors `dv`, `ev`, and `dp` such that

* `dv[i]` is the dimension of cell `i`
* `ev[k]` is the total number of cells with dimension `k-1`
* `dp[k]` is *1 plus* the number of cells of dimension *strictly less than* `k-1`

If in addition we have a nested sequence of complexes `E_0 ≤ ... ≤ E_n = E`, then let `fv` be the vector such that

* `fv[i]` is the birthtime of cell `i`

###### simple format (from file)

The formats described below are great for efficiency, but they can be hard to read.  Alternatively, you can write a comma separated file with the following format

```
dv[1], fv[1], <first face of cell 1>, <second face of cell 1>, ...
dv[2], fv[2], <first face of cell 2>, <second face of cell 2>, ...
...
dv[N], fv[N], <first face of cell N>, <second face of cell N>, ...
```

If this file is saved as `Users/Adam/ez.csv`, then to compute PH call

```
julia> C = eirene("Users/Adam/complex.csv",model="complex",entryformat="sp")
```


###### sparse column format

Eirene has a keyword argument for every vector defined above.  To compute PH:

```
julia> C = eirene(rv=rv,cp=cp,dv=dv,fv=fv)
```
You can use either `ev = ev` or `dp=dp` instead of `dv=dv`.  Eirene only needs one of the three.


###### sparse column format (from file)

You can store `E` as a comma separated .txt or .csv file with four lines

```
dp_1, ..., dp_m
fv_1, ..., fv_n
rv_1, ..., rv_p
cp_1, ..., cp_q
```

where line one is `dp`, line 2 is `fv`, line 3 is `rv`, and line 4 is `cp`.
* If you know `D` and want to figure out what `rv` and `cp` should be, you can either view the Julia docs on sparse matrices or google "sparse column format".
* Instead of `dp`, you can put either `dv` or `ev` in the first line.

Say this file is saved as `Users/Adam/complex.csv`.  To compute PH, call

```
julia> C = eirene("Users/Adam/complex.csv",model="complex",entryformat="dp")
```

* Replace `"dp"` with either `"ev"` or `"dv"`, depending on what you've placed in the first line.


### Barcodes, Betti Curves, and Persistence Diagrams
<br>
Eirene stores a list of every persistent homology class (that is, every bar) in the output variable `C`.  To find the birth and death times of each class, run

```
julia> A = barcode(C, dim=k)
```

Variable `A` will be an `n x 2` matrix, where `n` is the number of persistent homology classes of dimension `k`.  The pth persistent homology class in Eirene's list is born at time `A[p,1]` and dies at time `A[p,2]`.  To obtain the betti curve, run

```
julia> B = betticurve(C, dim=k)
```

Variable `B` is an `n x 2` array for which `B[j,2]` is the `k`th betti number of the space at time `t = B[j,1]`.

To plot the persistence diagram, betti curve, and barcode in dimension `k`, run the following.  Note  that `_pjs` simply refers to PlotlyJS, with is the Julia package used for graphing.

```
julia> plotpersistencediagram_pjs(C, dim=k)`
julia> plotbetticurve_pjs(C, dim=k)
julia> plotbarcode_pjs(C, dim=k)
```

* Keyword `dim` defaults to `1`.
* A class that never dies will appear as red dot on the diagonal at the time of its birth.
* Hovering the cursor over a point in the persistence diagram will show a message with its precise Euclidean coordinates (birth and death times) and two additional pieces of information:  `class` and `size`.  Variable `class` is an integer.  If `class = p`, then the point where the cursor is hovering represents the birth/death time of `p`th persistent homology class in Eirene's list.  It is born at time `A[p,1]` and dies at time `A[p,2]`.  Variable `size` is also an integer.  It is the number of cells in the cycle representative Eirene computed for this persistent homology class (since Eirene computes homology over the two element field, we always refer to vectors by their support).
* `class` and `size` can be displayed permanently, without the need to hover a cursor, via the keyword argument `showlabels = true`.

### Representatives

To compute a cycle representative for the `p`th persistent homology class in Eirene's list for dimension `k`, run

```
julia> S = classrep(C, class=p, dim=k)
```

By default `class = 1` and `dim = 1`.  Recall that, since Eirene computes homology over the two element field, we can always specify a cycle representative by listing the cells that support it.  If `C` was generated from a point cloud or Vietoris-Rips complex, then the cycle representative `S` will be encoded as an array with `k+1` rows.  Each column of `S` represents the set of vertices that make up a `k`-dimensional simplex.  If `C` was generated from an `nxn` distance matrix `M`, then the entries of `S` correspond to the rows/columns of `M`.  If `C` is generated from a matrix `M` whose columns represent points in a point cloud, then the elements of `S` refer to columns of `M`.  If `C` was generated with `model = "complex"`, then `S` will be a vector of integers.  These numbers refer to cells in the complex; more precisely, the number `q` refers to the `q`th  dimension-`k` cell stored in memory.  

Plotting for class representatives is currently supported only for point clouds and Vietoris-Rips complexes.  If the input data is a point cloud, coordinates for each vertex will be inferred automatically.  If the input data is not a point cloud, or if you would like to overwrite the existing coordinates, set  `C["inputdata"]["pointlabels"] = N` where `N` is a numeric array whose columns correspond to sample points.  Alternatively, you can pass `N` directly to the plotting function using keyword `coords = N` (this will not overwrite the existing coordinate data in `C`).

There are several ways to assign text labels to the points in your figure.  If labels were extracted from the input data or file, they are made available automatically for plotting.  To create or overwrite existing text labels, set `C["inputdata"]["pointlables"] = U`, where `U` is an array of type `Any` with all string entries.  `U` may also be passed directly to a plotting function using keyword `textlabels = U`.  If no other information is available, Eirene will label data points 1 through n, according to their order in the input array.  As in the case of persistence diagrams, both wrappers accept the `showlabels` keyword argument.  Hover-tags in PlotlyJS will still show text labels when `showlabels = false`.  A third option, `showlabels = "cycle"`, will cause `plotclassrep_pjs` to display the labels of only those vertices incident to the cycle representative being plotted.  

### Visualizing non-Euclidean Data

Multidimensional scaling (MDS) is a general approach to generating Euclidean embeddings of non-Euclidean or
high-dimensional data which (to the extent possible) faithfully represent distances between points. Most implementations of
the MDS algorithm accept a symmetric `m x m` matrix `S` and return an ordered set of `m` points in R<sup>n</sup>, with point `i` corresponding
to the `i`th column of `S`. These points determine an `m x m` matrix `D`, with `D[i,j]` equal to the Euclidean distance between points `i`
and `j`. The algorithm aims to to select the points that will make `D` as similar as possible to `S`. To use multidmensional scaling with
Eirene, pass `coords = "mds"` to the `plotclassrep_pjs` function. To specify an embedding of dimension `n`, use `embeddingdim = n`.


By default, the input matrix `S` will be the distance matrix used to generate the filtered complex. However in the special case
where only the vertices of a specific class representative are to be shown, that is, when the user specifies `showcloud =
false`, there an additional option. The class to be represented consists of a finite family of simplices, or faces, and the
dimension-1 simplices incident to these faces form a graph. The hop distance on this graph yields a new distance matrix,
which the user may specify with `embeddingobj = "hop"`. This option often creates cleaner representations.

### Formatting

#### Point clouds and Vietoris Rips complexes

PlotlyJS plots may be exported to the Plotly web API for sharing, storage, and advanced editing. Every function that generates
a PlotlyJS visual in the Eirene library returns a PlotlyJS object which may be uploaded via the Plotly. post function. See the
documentation for Plotly.jl to learn more.

Arrays can be read from (and written to) text flies using the native `readdlm` and `writedlm` (see the Julia
wikibooks chapter for text flies). The Julia base includes a number of additional functions for importing data from text flies,
and for convenience several of these have been combined in the Eirene wrapper `ezread(<pathtofile>)`.

Import from and export to Matlab flies are available from the packages MATLAB and MAT. Import from and export to Numpy
flies is managed by NPZ.

To save an entire dictionary `C`, the most convenient option by far is the Julia JLD package.

### Performance Tips

The marginal time and memory cost of computing explicit cycle representatives is often small relative to that of computing
homological barcodes alone (usually around 10%, for Euclidean point cloud data). Setting the `record = "none"`
may yield significant savings for other spaces, however.

### Comments
It is common among persistence solvers to round the entries of a user-provided distance matrix to a fixed number of distinct
values or significant digits. This is not the default behavior with Eirene, so care must be taken when comparing outputs.


Multidimensional scaling is a blue-chip method in data visualization, but no embedding algorithm is perfect. Particular care must be taken when `embeddingobj = "hop"` , as performance often decreases appreciably with the number of connected
components in the underlying graph. Remember that you can always inspect class representatives directly with the
`classrep` function.

### Advanced Usage

#####  Keyword arguments

The following are optional keyword arguments for the main funciton, `eirene`.

`pointlabels = A`
* If `A` is a nonempty array, then assign the pth vertex in the complex a label of `A[p]`, for all p.
* If `A` is empty, assign the pth vertex a label of `p`.
* Possible values: `A` may be any array of numbers or strings.  Default value: `[]`.

`fastop = true, false`
* A value of `true` allows Eirene to stop computing homology at a given time point `t` if it can be determined that the homology of the complex does not change after time `t`.  There is no reason to change this setting under normal circumstances.
* Possible values: `true`, `false`.  Default value: `true`.

##### Complex formats

Let us say that for any vector `v`, the symbol `v{k}` denotes the set of all integers `j` such that `v[k] ≤ j < v[k+1]` (you might recall that we used this convention when we explained how to format data for an arbitrary complex).  Now let us suppose we have a matrix `M` with 0 - 1 entries.  To save `M` in computer memory, it suffices to record which rows of each column have nonzero entries.  In particular, suppose we are given two vectors, `rv` (row values) and `cp` (column pattern).  If for all `p` the set `rv[cp{p}]` equals the set of rows where column `p` is nonzero, then we can (essentially) record `M` in computer memory by saving `rv` and `cp`.  We call this the **column sparse format** of `M`.

When we run the command
```
C = eirene(rv=rv, cp=cp, fv=fv, dp=dp, <keyword arguments>)
```
we are passing `eirene` a pair of vectors (`rv` and `cp`) that record the *total boundary matrix*; this matrix has a row and column for every cell of every dimension.  The function of the vector `dp` is to record the dimension of each cell, and the function of `fv` is to specify the birth time of each cell.  

When Eirene processes these data, they are broken down by dimension.  Eirene creates *arrays* `RV`, `CP`, and `FV` such that (1) `RV[p]` and `CP[p]` are column sparse formats for the boundary matrix with columns indexed by the cells of dimension `p` and rows indexed by cells of dimension `p-1`, and (2) `FV[p][q]` is the birth time of the `q`th `p-1`-dimensional cell.  If you already have access to `RV`, `CP`, and `FV` you can pass them to Eirene directly with the command
```
C = eirene(rv=RV, cp=CP, fv=FV, dp=[], <keyword arguments>)
```



### Common Issues

The most common source of technical difficulties in working with Eirene is software update. The best way to check that you are running the most recent version of Julia is to log onto the Julia homepage. To ensure that you have the most recent
version of the the various support packages, enter `Pkg. update()` in a Julia REPL.

The second most common source of technical difficulties rregards operating systems. The fewest issues by far have been
reported for Mac OS X. If you have already visited the Julia homepage to make sure your Julia install is up-to-date, it may be worthwhile switching to Mac to see if your issue persists. Windows users have reported encountering an error message that includes the phrase "invalid escape sequence" in reference to file paths. We found that replacing the unary operator `\` with `\\` generally resolves this issue, e.g. `C:\\Users` versus `C:\Users`.

Calling the include function twice with the same input will frequently produce a large number of error messages. This is a well documented issue in the current Julia distribution, and does not impact performance.

### Documenting Use

If Eirene has helped your research, teaching, or artistic endeavors, please let us know! You can write us here and copy/paste the following bibtex entry to cite the original tech report. If you have ideas about how Eirene could help with education and
outreach, please let us know! We would love to hear from you.

@ARTICLE{henselmanghristl6,
	<br>
author= {{Henselman}, G. and {Ghrist}, R.},
	<br>
	title= "{Matroid Filtrations and Computational  Persistent Homology}",
	<br>
	journal= {ArXiv e-prints}, <br>
archivePrefix = "arXiv",
<br>
eprint = {1606.00199},
<br>
primaryClass = "math.AT",
<br>
keywords= {Mathematics -Algebraic Topology, Mathematics -Combinatorics},
<br>
year= 2016,
<br>
month= jun,
<br>
adsurl = {http://adsabs.harvard.edu/abs/2016arXiv160600199H},
<br>
 adsnote = {Provided by the SAO/NASA Astrophysics Data System}
 <br>
}






Eirene.jl package for Homological Algebra

Copyright (C) 2016, 2017, 2018  Gregory Henselman

www.gregoryhenselman.org

Eirene is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Eirene is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Eirene.  If not, see <http://www.gnu.org/licenses/>.

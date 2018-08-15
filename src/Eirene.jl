__precompile__()

module Eirene

#     You should have received a copy of the GNU General Public License
#     along with Eirene.  If not, see <http://www.gnu.org/licenses/>.

print("\n
Eirene Library for Homological Algebra
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

")
print_with_color(:blue,"
WELCOME TO EIRENE!
v$(Pkg.installed("Eirene"))

Please help us document Eirene's recent work! Bibtex entries and
contact information for teaching and outreach can be found at the
Eirene homepage, http://gregoryhenselman.org/eirene.\n\n
")

##########################################################################################

#### 	REQUIREMENTS

##########################################################################################

if typeof(Pkg.installed("Distances")) == Void
print_with_color(:green,"Please Note: Distances.jl may not be installed. This package is required
for use with Euclidean point cloud data. To install, enter the foll-
owing at the Julia prompt:

Pkg.add(\"Distances\")
using Distances \n\n
")
else
	using Distances
end

if typeof(Pkg.installed("JLD")) == Void
print_with_color(:green,"Please Note: JLD.jl may not be installed. This package is not required, but
it is the best means of saving Eirene output. To install, enter the
following at the Julia prompt:

Pkg.add(\"JLD\")
using JLD \n\n
")
else
	using JLD
end

if typeof(Pkg.installed("Blink")) == Void
print_with_color(:green,"Please Note: Blink.jl may not be installed. This package is required for
all Eirene functions ending in _pjs. To install, enter the following
at the Julia prompt:

Pkg.add(\"Blink\")
using Blink
Blink.AtomShell.install() \n\n
")
else
	using Blink
end

if typeof(Pkg.installed("PlotlyJS")) == Void
print_with_color(:green,"Please Note: PlotlyJS.jl may not be installed. This package is required
for all Eirene functions ending in _pjs. To install, enter the foll-
owing at the Julia prompt

Pkg.add(\"PlotlyJS\")
using PlotlyJS \n\n
")
else
	using PlotlyJS
end

if typeof(Pkg.installed("Plotly")) == Void
print_with_color(:green,"Please Note: Plotly.jl may not be installed. This package is required to
interface with the Plotly web API. To install, enter the foll-
owing at the Julia prompt

Pkg.add(\"Plotly\")
using Plotly \n\n
")
else
	using Plotly
end

if typeof(Pkg.installed("MultivariateStats")) == Void
print_with_color(:green,"Please Note: MultivariateStats.jl may not be installed. This package is required for
some operations pertaining to multidimensional scaling, but is not required.
To install, enter the following at the Julia prompt:

Pkg.add(\"MultivariateStats\")
using MultivariateStats \n\n
")
else
	using MultivariateStats
end

if typeof(Pkg.installed("Colors")) == Void
print_with_color(:green,"Please Note: Colors.jl may not be installed. This package is required for
some operations pertaining to plotting barcodes, but is not required.
To install, enter the following at the Julia prompt:

Pkg.add(\"Colors\")
using MultivariateStats \n\n
")
else
	using Colors
end

##########################################################################################

#### 	EXPORTS

##########################################################################################

export 	eirene,
		eirenefilepath,
		ezread,
		ezplot_pjs,
		plotpersistencediagram_pjs,
		plotclassrep_pjs,
		plotbarcode_pjs,
		plotbetticurve_pjs,
		ezplot_pjs,
		barcode,
		betticurve,
		classrep,
		latlon2euc,
		eirenefilepath,
		noisycircle,
		noisycircle3,
		torus,
		noisytorus,
		sphere,
		matchingcomplex_symmat,
		chessboardcomplex_symmat,
		plane2torus,
		zerodrandmat,
		ezlabel,
		unittest

##########################################################################################

#### 	SIMPLICIAL CONSTRUCTIONS

##########################################################################################

# NB: eirene permutes vertex labels prior to calculation;
# all versions of <vertexrealization> must be interpreted with
# respect to the PERMUTED labeling scheme
function vertexrealization(farfaces,firstv,facecardinality,facenames)
	numfaces::Int64 = length(facenames)
	fc::Int64 = facecardinality
	m::Int64 = length(firstv[2])-1
	preallocationspace = 0
	loci::Array{Int64,1} = copy(facenames)
	vrealization = Array{Int64}(facecardinality,numfaces)
	post0::Int64 = 1
	post1::Int64 = 1

	for sd = facecardinality:-1:1
		cp::Array{Int64,1} = firstv[sd]
		for i = 1:numfaces
			locus = loci[i]
			if cp[post0] > locus
				while cp[post0] > locus
					post0-=1
				end
				post1 = post0+1
			elseif cp[post1] <= locus
				while cp[post1] <= locus
					post1+=1
				end
				post0 = post1-1
			end
			loci[i] = farfaces[sd][locus]
			vrealization[sd,i]=post0
		end
	end
	return vrealization
end

# NB: eirene permutes vertex labels prior to calculation;
# all versions of <vertexrealization> must be interpreted with
# respect to the PERMUTED labeling scheme
function vertexrealization(D::Dict,facecardinality,facenames)
	return vertexrealization(D["farfaces"],D["firstv"],facecardinality,facenames)
end

# NB: eirene permutes vertex labels prior to calculation;
# all versions of <vertexrealization> must be interpreted with
# respect to the PERMUTED labeling scheme
function vertexrealization(D::Dict;dim = 1, class = 1)
	sd = dim+2
	facecard = dim+1

	if haskey(D,"cyclerep")
		rep = D["cyclerep"][sd][class]
	else
		cyclename = barname2cyclename(D,class;dim = dim)
		rep = getcycle(D,sd,cyclename)
	end

	vrealization = vertexrealization(D::Dict,facecard,rep)
end

# NB: eirene permutes vertex labels prior to calculation;
# the output of <incidentverts> must be interpreted with respect
# to the PERMUTED labeling scheme
function incidentverts(farfaces,firstv,facecardinality,facenames)
	numfaces::Int64 = length(facenames)
	fc::Int64 = facecardinality
	m::Int64 = length(firstv[2])-1
	preallocationspace = 0
	vsupp = falses(m)
	loci::Array{Int64,1} = copy(facenames)
	post0::Int64 = 1
	post1::Int64 = 1
	for sd = facecardinality:-1:1
		cp::Array{Int64,1} = firstv[sd]
		for i = 1:numfaces
			locus = loci[i]
			if cp[post0] > locus
				while cp[post0] > locus
					post0-=1
				end
				post1 = post0+1
			elseif cp[post1] <= locus
				while cp[post1] <= locus
					post1+=1
				end
				post0 = post1-1
			end
			loci[i] = farfaces[sd][locus]
			vsupp[post0]=true
		end
	end
	return find(vsupp)
end

# NB: eirene permutes vertex labels prior to calculation;
# the output of <incidentverts> must be interpreted with respect
# to the PERMUTED labeling scheme
function incidentverts(D::Dict,facecardinality,facenames)
	return incidentverts(D["farfaces"],D["firstv"],facecardinality,facenames)
end

# NB: eirene permutes vertex labels prior to calculation;
# the output of <incidentverts> must be interpreted with respect
# to the PERMUTED labeling scheme
function incidentverts(D::Dict;dim=1,class=1)
	facecardinality = dim+1

	if haskey(D,"cyclerep")
		rep = D["cyclerep"][dim+2][class]
	else
		cyclename = barname2cyclename(D,class;dim=dim)
		rep = getcycle(D,facecardinality,class)
	end
	return incidentverts(D,facecardinality,rep)
end

function buildclosefromclose(lrowval,lcolptr,lclosefaces,hrowval,hcolptr;facecard = size(lclosefaces,1)+1)
	m = length(hcolptr)-1
	n = length(hrowval)
	hclosefaces = Array{Int64}(facecard,n)
	if n == 0
		return hclosefaces
	else
		rowdepth = facecard-1
		rosettacol = Array{Int64}(maximum(lrowval))
		for i = 1:m
			rosettacol[lrowval[cran(lcolptr,i)]]=cran(lcolptr,i)
			for j = cran(hcolptr,i)
				farface = hrowval[j]
				for k = 1:rowdepth
					hclosefaces[k,j]=rosettacol[lclosefaces[k,farface]]
				end
				hclosefaces[facecard,j] = rosettacol[lrowval[farface]]
			end
		end
		return hclosefaces
	end
end

function buildclosefromclose_subr(rosettacol::Array{Int64,1},lrowval::Array{Int64,1},lcolptr::Array{Int64,1},hrowval::Array{Int64,1},hcolptr::Array{Int64,1},hclosefaces::Array{Int64,1},columnmarker::Int64,rowdepth::Integer)
	for i = 1:m
		rosettacol[lrowval[cran(lcolptr,i)]]=cran(lcolptr,i)
		for j = cran(hcolptr,i)
			if columnsupp[j]
				columnmarker+=1
				farface = hrowval[j]
				buildclosefromclose_subr_subr(rowdepth::Integer,hclosefaces::Array{Int64,1},columnmarker::Int64,rowettacol::Array{Int64,1},lclosefaces::Array{Int64,1},farface::Int64)
				hclosefaces[sd,columnmarker] = rosettacol[lrowval[farface]]
			end
		end
	end
end

function buildclosefromclose_subr_subr(rowdepth::Integer,hclosefaces::Array{Int64,1},columnmarker::Int64,rowettacol::Array{Int64,1},lclosefaces::Array{Int64,1},farface::Int64)
	for k = 1:rowdepth
		hclosefaces[k,columnmarker]=rosettacol[lclosefaces[k,farface]]
	end
end

function buildallfromclose(lrowval,lcolptr,lclosefaces,hrowval,hcolptr,selectedcolumnindices;verbose=false)
	if verbose
		println("PLEASE NOTE: COLUMNS MUST BE IN SORTED ORDER FOR THIS TO WORK PROPERLY")
	end
	m = length(hcolptr)-1
	numhigs = length(hrowval)
	numselected = length(selectedcolumnindices)
	rowdepth = size(lclosefaces,1)
	sd = rowdepth+1
	hclosefaces = Array{Int64}(sd+1,numselected)
	if numselected == 0
		return hclosefaces
	end
	rosettacol = Array{Int64}(maximum(lrowval))
	columnsupp = falses(numhigs)
	columnsupp[selectedcolumnindices]=true
	columnmarker = 0
	for i = 1:m
		rosettacol[lrowval[cran(lcolptr,i)]]=cran(lcolptr,i)
		for j = cran(hcolptr,i)
			if columnsupp[j]
				columnmarker+=1
				farface = hrowval[j]
				for k = 1:rowdepth
					hclosefaces[k,columnmarker]=rosettacol[lclosefaces[k,farface]]
				end
				hclosefaces[sd,columnmarker] = rosettacol[lrowval[farface]]
				hclosefaces[sd+1,columnmarker] = farface
			end
		end
	end
	return hclosefaces
end

function buildclosefaces!(lrowval,lcolptr,lclosefaces,lfarfaces,hrowval,hcolptr,destinationmatrix;verbose = false)
	m = length(hcolptr)-1
	n = length(hrowval)
	rowdepth = size(lclosefaces,1)
	sd = rowdepth+1
	rosettacol = Array{Int64}(maximum(lrowval))
	for i = 1:m
		rosettacol[lrowval[cran(lcolptr,i)]]=cran(lcolptr,i)
		for j = cran(hcolptr,i)
			farface = hrowval[j]
			for k = 1:rowdepth
				destinationmatrix[k,j]=rosettacol[lclosefaces[farface]]
			end
			destinationmatrix[sd,j] = rosettacol[lrowval[farface]]
		end
	end
	for j = 1:n
		for i = 1:sd
			lclosefaces[i,j]=destinationmatrix[i,j]
		end
	end
end

function buildclosefromfar(farfaces,firstv,sd)
	m = length(firstv[1])-1
	n = length(farfaces[sd])
	# destinationmatrix = Array{Int64}(sd,n)
	if sd == 1
		return Array{Int64}(0,m)
	end
	lclosefaces = Array{Int64}(1,firstv[2][end]-1)
	for i = 1:m
		lclosefaces[cran(firstv[2],i)]=i
	end
	if sd == 2
		return lclosefaces'
	end
	for i = 3:sd
		lclosefaces = buildclosefromclose(farfaces[i-1],firstv[i-1],lclosefaces,farfaces[i],firstv[i];facecard=i-1)
	end
	return lclosefaces
end

function buildclosefromfar(farfaces,firstv,sd,columnsinorder)
	m = length(firstv[1])-1
	n = length(farfaces[sd])
	if sd == 1
		return Array{Int64}(0,m)
	end
	lclosefaces = Array{Int64}(1,firstv[2][end]-1)
	for i = 1:m
		lclosefaces[cran(firstv[2],i)]=i
	end
	if sd == 2
		return lclosefaces[columnsinorder]'
	end
	for i = 3:(sd-1)
		lclosefaces = buildclosefromclose(farfaces[i-1],firstv[i-1],lclosefaces,farfaces[i],firstv[i];facecard=i-1)
	end
	lclosefaces = buildclosefromclose(farfaces[sd-1],firstv[sd-1],lclosefaces,farfaces[sd],firstv[sd],columnsinorder;facecard = sd-1)
	return lclosefaces
end

function buildallfromfar(farfaces,firstv,sd,columnsinorder;verbose = false)
	m = length(firstv[1])-1
	n = length(farfaces[sd])
	# destinationmatrix = Array{Int64}(sd,n)
	if sd == 1
		return Array{Int64}(0,m)
	end
	lclosefaces = Array{Int64}(1,firstv[2][end]-1)
	for i = 1:m
		lclosefaces[cran(firstv[2],i)]=i
	end
	if sd == 2
		return vcat(lclosefaces[columnsinorder]',farfaces[sd][columnsinorder]')
	end
	for i = 3:(sd-1)
		lclosefaces = buildclosefromclose(farfaces[i-1],firstv[i-1],lclosefaces,farfaces[i],firstv[i];facecard=i-1)
		#gc()
	end
	lclosefaces = buildallfromclose(farfaces[sd-1],firstv[sd-1],lclosefaces,farfaces[sd],firstv[sd],columnsinorder;verbose = verbose)
	#gc()
	return lclosefaces
end

function ff2boundary(farfaces,firstv;sd=1)
	rv = Array{Int64}(0)
	cp = [1]
	if sd == 1
		rv = Array{Int64}(0)
		cp = ones(Int64,length(farfaces[1])+1)
	else
		n = length(farfaces[sd])
		rv = ff2aflight(farfaces,firstv,sd,1:n)
		rv = vec(rv)
		cp = convert(Array{Int64,1},sd+1:sd:(1+n*sd))
		prepend!(cp,[1])
	end
	return rv,cp
end

function ff2complex(farfaces,firstv;maxsd = length(farfaces))
	Nrv 	= fill(Array{Int64}(0),maxsd)
	Ncp 	= fill(Array{Int64}(0),maxsd)
	Nrv		= convert(Array{Array{Int64,1}},Nrv)
	Ncp		= convert(Array{Array{Int64,1}},Ncp)
	Nrv[1] 	= Array{Int64}(0)
	Ncp[1]	= fill(1,length(farfaces[1])+1)
	for sd = 2:maxsd
		Nrv[sd],Ncp[sd] = ff2boundary(farfaces,firstv,sd=sd)
	end
	return Nrv,Ncp
end

function eirened2complex(C)
	if in(C["input"]["model"],["pc","vr"])
		rv,cp 	= 	boundarymatrices(C)
		fv 		= 	ocff2of(C["grain"],C["ocg2rad"])
	elseif in(C["input"]["model"],["complex"])
		rv 		= 	C["rv"]
		cp 		= 	C["cp"]
	else
		println("Error: the value of C[\"input\"][\"model\"] must be \"pc\", \"vr\", or \"complex\".")
	end
	fv 			= 	ocff2of(C["grain"],C["ocg2rad"])
	return 		rv,cp,fv
end

function ocff2of(grain::Array{Int64},ocg2rad::Array{Int64})
	m = length(grain)
	filt = Array{Int64}(m)
	for i = 1:m
		filt[i] = ocg2rad[grain[i]]
	end
	return filt
end

function ff2aflight_sc2(farfaces,firstv,columns)
	sd = 2
	if isempty(farfaces[sd])
		return Array{Int64}(2,0)
	end
	f0faces::Array{Int64,1} = farfaces[sd]
	colptr::Array{Int64,1} = firstv[2]
	columnpost::Int64   = 1
	columnpostp1::Int64 = 2
	faces::Array{Int64,2} = Array{Int64}(2,length(columns))

	for fp = 1:length(columns)
		f0 = columns[fp]
		if f0 >= colptr[columnpostp1]
			while f0 >= colptr[columnpostp1]
				columnpostp1+=1
			end
			columnpost = columnpostp1-1
		elseif f0 < colptr[columnpost]
			while f0 < colptr[columnpost]
				columnpost-=1
			end
			columnpostp1 = columnpost+1
		end
		faces[1,fp] = columnpost
		faces[2,fp] = f0faces[f0]
	end
	return faces
end

function ff2aflight_sc3(farfaces,firstv,columns)
	sd = 3

	if isempty(farfaces[sd])
		return Array{Int64}(3,0)
	end

	fcfaces::Array{Int64,2} = buildclosefromfar(farfaces,firstv,sd-1,1:length(farfaces[2]))

	f0faces::Array{Int64,1} = farfaces[sd]
	f1faces::Array{Int64,1} = farfaces[sd-1]

	fvscm0::Array{Int64,1}  = firstv[sd]
	fvscm1::Array{Int64,1}  = firstv[sd-1]
	fvscm2::Array{Int64,1}  = firstv[sd-2]

	holdi=[1];holdip1=[2]
	t1::Array{Int64,1} = Array{Int64}(fvscm2[end]-1);t1[crows(fvscm1,f1faces,1)]=cran(fvscm1,1)

	faces::Array{Int64,2} = Array{Int64}(3,length(columns))
	for fp = 1:length(columns)
		f0 = columns[fp]
		f1 = f0faces[f0]
		f2 = f1faces[f1]
		f3 = fcfaces[f1]
		updatetranslator!(f0::Int64,fvscm0::Array{Int64,1} ,holdi::Array{Int64,1},holdip1::Array{Int64,1},t1::Array{Int64,1},fvscm1::Array{Int64,1},f1faces::Array{Int64,1})
		faces[1,fp] = t1[f3]
		faces[2,fp] = t1[f2]
		faces[3,fp] = f1
	end
	return faces
end

function ff2aflight_scgt3(farfaces,firstv,sd,columns)

	if isempty(farfaces[sd])
		return Array{Int64}(sd,0)
	end

	f0faces::Array{Int64,1} = farfaces[sd]
	f1faces::Array{Int64,1} = farfaces[sd-1]
	f2faces::Array{Int64,1} = farfaces[sd-2]
	fcfaces::Array{Int64,2} = buildallfromfar(farfaces,firstv,sd-2,1:(firstv[sd-2][end]-1))

	fvscm0::Array{Int64,1}  = firstv[sd]
	fvscm1::Array{Int64,1}  = firstv[sd-1]
	fvscm2::Array{Int64,1}  = firstv[sd-2]
	fvscm3::Array{Int64,1}  = firstv[sd-3]

	holdi=[1];holdip1=[2];holdj=[1];holdjp1=[2]
	t1::Array{Int64,1} = Array{Int64}(fvscm2[end]-1);t1[crows(fvscm1,f1faces,1)]=cran(fvscm1,1)
	t2::Array{Int64,1} = Array{Int64}(fvscm3[end]-1);t2[crows(fvscm2,f2faces,1)]=cran(fvscm2,1)

	scm0::Int64 = sd; scm1::Int64 = sd-1; scm2::Int64 = sd-2
	faces::Array{Int64,2} = Array{Int64}(sd,length(columns))

	for fp = 1:length(columns)
		f0 = columns[fp]
		f1 = f0faces[f0]
		f2 = f1faces[f1]
		updatetranslator!(f0::Int64,fvscm0::Array{Int64,1},holdi::Array{Int64,1},holdip1::Array{Int64,1},t1::Array{Int64,1},fvscm1::Array{Int64,1},f1faces::Array{Int64,1})
		updatetranslator!(f1::Int64,fvscm1::Array{Int64,1},holdj::Array{Int64,1},holdjp1::Array{Int64,1},t2::Array{Int64,1},fvscm2::Array{Int64,1},f2faces::Array{Int64,1})
		for i = 1:scm2
			faces[i,fp] = t1[t2[fcfaces[i,f2]]]
		end
		faces[scm1,fp] = t1[f2]
		faces[scm0,fp] = f1
	end
	return faces
end

function updatetranslator!(f0,firstv0,holdi,holdip1,t,firstv1,farfaces1)
	if firstv0[holdip1[1]] <= f0
		while firstv0[holdip1[1]]<= f0
			holdip1[1]+=1
		end
		holdi[1] = holdip1[1]-1
		t[crows(firstv1,farfaces1,holdi[1])]=cran(firstv1,holdi[1])
	elseif firstv0[holdi[1]] > f0
		while firstv0[holdi[1]] > f0
			holdi[1]-=1
		end
		holdip1[1] = holdi[1]+1
		t[crows(firstv1,farfaces1,holdi[1])]=cran(firstv1,holdi[1])
	end
end

function ff2aflight_subr!(
	columns::UnitRange{Int64},f0faces::Array{Int64,1},f1faces::Array{Int64,1},f2faces::Array{Int64,1},fcfaces::Array{Int64,2},fvscm0::Array{Int64,1},
	fvscm1::Array{Int64,1},fvscm2::Array{Int64,1},holdi::Array{Int64,1},holdip1::Array{Int64,1},holdj::Array{Int64,1},
	holdjp1::Array{Int64,1},t1::Array{Int64,1},t2::Array{Int64,1},faces::Array{Int64,2},
	scm0::Int64,scm1::Int64,scm2::Int64)
	for fp = 1:length(columns)
		f0 = columns[fp]
		f1 = f0faces[f0]
		f2 = f1faces[f1]
		updatetranslator!(f0::Int64,fvscm0::Array{Int64,1} ,holdi::Array{Int64,1},holdip1::Array{Int64,1},t1::Array{Int64,1},fvscm1::Array{Int64,1},f1faces::Array{Int64,1})
		updatetranslator!(f1::Int64,fvscm1::Array{Int64,1},holdj::Array{Int64,1},holdjp1::Array{Int64,1},t2::Array{Int64,1},fvscm2::Array{Int64,1},f2faces::Array{Int64,1})
		for i = 1:scm2
			faces[i,fp] = t1[t2[fcfaces[i,f2]]]
		end
		faces[scm1,fp] = t1[f2]
		faces[scm0,fp] = f1
	end
end

function ff2aflight(farfaces,firstv,sd,columns)
	if sd == 1
		return Array{Int64}(0,length(columns))
	elseif sd == 2
		return ff2aflight_sc2(farfaces,firstv,columns)
	elseif sd == 3
		return ff2aflight_sc3(farfaces,firstv,columns)
	else
		return ff2aflight_scgt3(farfaces,firstv,sd,columns)
	end
end

function ff2aflight(D::Dict,sd,columns)
	farfaces = D["farfaces"]; firstv = D["firstv"]
	faces = ff2aflight(farfaces,firstv,sd,columns)
	return faces
end

#=

NB
- Input argument <grain> must be arranged least to greatest

OUTPUTS
- higlab
The concatenated vector [pphigs,nphigs)]
- lowlab
The concatenated vector [pplows,nplows[perm]], where perm is a permutation such
that the entries of lowgrain[nplows[perm]] appear in ascending order, numer-
ically.
- Mrv, Mcp, Mm
Sparse matrix representation of transpose(D[lowlab,higlab]), where D is
submatrix of the total boundary operator indexed by cells of dimension sd-1
(along the columns) and sd-2 (along the rows).

=#
function filteredmatrixfromfarfaces{Tv}(
	farfaces,
	firstv,
	prepairs,
	grain,
	sd::Integer,
	lowbasisnames::Array{Tv,1};
	verbose = false)

	numhigs = length(farfaces[sd])
	numlows = length(farfaces[sd-1])
	numppair= length(prepairs[sd])

	pphigs = prepairs[sd]
	pplows = farfaces[sd][pphigs]
  	lpls = lowbasisnames
	hphs = farfaces[sd+1][prepairs[sd+1]]
	nplows = intervalcomplementuniqueunsortedinput(vcat(lpls,pplows),numlows)
	nphigs = intervalcomplementuniqueunsortedinput(vcat(hphs,pphigs),numhigs)

	numnhph = numhigs-length(hphs)
	Ml = numlows - length(lpls)
	Mh = numhigs - length(hphs)

	higtranslator = zeros(Tv,numnhph)
	lowtranslator = zeros(Tv,numlows)
	lowtranslator[pplows] = 1:numppair

	# if !isempty(nplows) && sd > 2
	# 	npfilt = grain[sd-1][nplows]
	# 	nporder = integersinsameorder(npfilt)
	# 	addinteger!(nporder,numppair)
	# else
	# 	nporder = (numppair+1):(numppair+length(nplows))
	# end

	if sd > 1
		npfilt 	= 	grain[sd-1][nplows]
		nporder = 	integersinsameorder(npfilt)
		addinteger!(nporder,numppair)
	else
		npfilt 	= 	zeros(Int64,0)
		nporder =	zeros(Int64,0)
	end

	lowtranslator[nplows] = nporder
	higsinpointorder = intervalcomplementuniqueunsortedinput(hphs,numhigs)
	lowlab = Array{Int64}(Ml)
	lowlab[1:numppair]=pplows
	lowlab[nporder]=nplows
	higlab = vcat(pphigs,nphigs)

	if verbose
		comparisonsuppvec = trues(numhigs)
		comparisonsuppvec[hphs]=false
		comparisonvec=find(comparisonsuppvec)
		differencecounter = 0
		for i = 1:length(higsinpointorder)
			if higsinpointorder[i]!=comparisonvec[i]
				differencecounter+=1
			end
		end
		if differencecounter>0
			print(["hi ho comparison vec" differencecounter])
			print(length(higsinpointorder))
			print(length(comparisonvec))
			print(comparisonvec[1:20])
			print(higsinpointorder[1:20])
			sleep(5)
		end
	end
	ppsupp = falses(numhigs)
	ppsupp[pphigs]=true
	ppmarker = 0
	nppmarker = numppair
	for i = 1:numnhph
		hig = higsinpointorder[i]
		if ppsupp[hig]
			ppmarker+=1
			higtranslator[i]=ppmarker
		else
			nppmarker+=1
			higtranslator[i]=nppmarker
		end
	end
	allfaces = buildallfromfar(farfaces,firstv,sd,higsinpointorder;verbose = verbose)
	if verbose
		print("done building allfromfar")
	end
	Mrv,Mcp,Mm = presparsefull2unsortedsparsetranspose(allfaces,lowtranslator,higtranslator;verbose=verbose)
	higtranslator = [];npfilt = [];ppsupp = [];allfaces = []
	#gc()
	if verbose && length(Mrv)>(Mcp[end]-1)
		print("There was the thought that Mrv should have no extra elements")
		sleep(3)
	end
	return Mrv,Mcp,lowlab,higlab,Mm
end

function getmaxdim(farfaces)
	l = length(farfaces)
	maxdim = l
	for i = 1:l
		if length(farfaces[i]) == 0
			maxdim = i
			break
		end
	end
	return maxdim
end

function grain2maxsd(grain)
	c = 0
	for i = 1:length(grain)
		if !isempty(grain[i])
			c = i
		end
	end
	return c
end

function skelcount(numvertices,maxsdinality)
	c = 0
	for i = 1:maxsdinality
		c += binom(numvertices,i)
	end
	return c
end

##########################################################################################

####	CELL OPERATIONS

##########################################################################################

function cellcount(C)
	c = 0
	for i = 1:length(C["grain"])
		c+= length(C["grain"][i])
	end
	return c
end

function ocff2of(grain::Array{Int64},ocg2rad::Array{Float64})
	m = length(grain)
	filt = Array{Float64}(m)
	for i = 1:m
		filt[i] = ocg2rad[grain[i]]
	end
	return filt
end

function ocff2of(grain::Array{Array{Int64,1},1},ocg2rad::Array{Float64})
	n = length(grain)
	filt = Array{Array{Float64}}(n)
	for i = 1:n
		filt[i] = ocff2of(grain[i],ocg2rad)
	end
	return filt
end

function floatgrain(C)
	grain = convert(Array{Array{Int64}},C["grain"])
	ocg2rad = C["ocg2rad"]
	return ocff2of(grain,ocg2rad)
end

##########################################################################################

####	SCHUR COMPLEMENTS

##########################################################################################

function 	schurit4!(	Mrv,Mcp,Mm,Mn,Mn0,
						rowlab,collab,
						Jprows,Jpcols,numjunpairs,
						Sprows,Spcols,numsenpairs,
						comprows,compcols,
						Trv,Tcp,Srv,Scp;
						updatetransform = true,
						verbose = false,
						diagonstic = false
					)

	Mm0 = copy(Mm[1])
	Mm[1] = length(comprows)
	Mn[1] = length(compcols)

	copycolind2colind!(Trv,Tcp,Jpcols[numjunpairs[1]:-1:1],Srv,Scp,numsenpairs[1]+1,0)

	topspot = numsenpairs[1]+numjunpairs[1]
	for i = 1:numjunpairs[1]
		placementlocation = numsenpairs[1]+i
		extractionlocation = numjunpairs[1]-i+1
		Sprows[placementlocation] = rowlab[Jprows[extractionlocation]]
		Spcols[placementlocation] = collab[Jpcols[extractionlocation]]
	end

	numsenpairs[1]+=numjunpairs[1]

	rowsum = falses(Mm0)#zeros(Int64,Mm0)
	for jp = 1:Mn[1]
		j = compcols[jp]
		for ip = cran(Mcp,j)
			rowsum[Mrv[ip]]=true
		end
	end

	keptlist = finddownstreamelements_embeddedupperunitriangularmatrix(
					Mrv,Mcp,Mm0,find(rowsum),Jprows[1:numjunpairs[1]],Jpcols[1:numjunpairs[1]];verbose=verbose
					)

	if verbose
		println()
		println([length(keptlist) "=numkept" numjunpairs[1] "=numinputp" Mn[1] "=Mn" Mm[1] "=Mm" (Mcp[Mn[1]+1]-1) "=nnz(M)"])
	end

	keptmarker = length(keptlist)
	prows = Array{Int64}(keptmarker)
	pcols = Array{Int64}(keptmarker)
	for i = 1:keptmarker
		keptindex = keptlist[i]
		prows[i] = Jprows[keptindex]
		pcols[i] = Jpcols[keptindex]
	end

	Arv,Crv,Acp,Ccp = stackedsubmatrices(Mrv,Mcp,prows,comprows,pcols,Mm0)
	Brv,Drv,Bcp,Dcp = stackedsubmatrices(Mrv,Mcp,prows,comprows,compcols,Mm0)
	Lrv,Lcp = copycolumnsubmatrix(Trv,Tcp,pcols)
	Rrv,Rcp = copycolumnsubmatrix(Trv,Tcp,compcols)

	translator = Array{Int64}(Mm0)
	translator[prows]=1:keptmarker
	yafterx!(translator,Arv)
	yafterx!(translator,Brv)

	translator[comprows]=1:length(comprows)
	yafterx!(translator,Crv)
	yafterx!(translator,Drv)

	Airv,Aicp = morseInverseF2orderedColsUnsortedRowsInSilentOut(Arv,Acp)
	Brv,Bcp = spmmF2silentLeft(Airv,Aicp,Brv,Bcp,keptmarker)

	for j = 1:keptmarker
		translator[j]=collab[pcols[j]] #repurposing rowsum as col-to-row translator
	end
	collabcopy = copy(collab)
	for i = 1:Mm[1]
		rowlab[i] = rowlab[comprows[i]]
	end
	for j = 1:Mn[1]
		collab[j] = collab[compcols[j]]
	end

	blockprodsumWrite2!(Crv,Ccp,Brv,Bcp,Drv,Dcp,Mm[1],Mn[1],Mrv,Mcp,0)

	if updatetransform
		blockprodsumsilenticolsleftWrite2!(Lrv,Lcp,Brv,Bcp,Rrv,Rcp,Mn0,Mn[1],Trv,Tcp,0,translator)
	end
end

##########################################################################################

####	CHAIN OPERATIONS

##########################################################################################

#=
notes on morselu!
- 	the output array tid has size equal to the number of columns of the
	input array, NOT the column-rank of the input array
-	the columns of M must be ordered by grain (ascending)
-	the first (rank of M) elements of tlab index the complete set
	of nonzero columns in the reduced matrix
=#
function morselu!{Tv<:Integer}(
	Mrv::Array{Tv,1},
	Mrowgrain::Array{Tv,1},
	Mcp::Array{Tv,1},
	Mcolgrain::Array{Tv,1},
	lowlab::Array{Tv,1},
	higlab::Array{Tv,1},
	pplow::Array{Tv,1},
	pphig::Array{Tv,1},
	Mm::Integer;
	storetransform = true,
	verbose = false,
	diagnostic = false)

 	rowlab = higlab;collab = lowlab

	Mm = [length(higlab)]
	Mn = [length(lowlab)]
	Mn0 = Mn[1]
	maxnz = Mcp[Mn[1]+1]

	maxnumpairs = min(Mm[1],Mn[1]); numjunpairs = [length(pplow)]; numsenpairs = [0]
	Sprows=Array{Tv}(maxnumpairs);Spcols=Array{Tv}(maxnumpairs);
	Jprows=Array{Tv}(maxnumpairs);Jpcols=Array{Tv}(maxnumpairs);
	Jprows[1:numjunpairs[1]]=pphig;Jpcols[1:numjunpairs[1]]=pplow
	comprows = convert(Array{Tv,1},(numjunpairs[1]+1):Mm[1])
	compcols = convert(Array{Tv,1},(numjunpairs[1]+1):Mn[1])

	Trv=Array{Tv}(0);Srv = Array{Tv}(0)
	Tcp=ones(Tv,Mn[1]+1);Scp=ones(Tv,Mn[1]+1)

	if diagnostic
		numsenpairsOLD = numsenpairs[1]
	end
	schurit4!(		Mrv,Mcp,Mm,Mn,Mn0,
					rowlab,collab,
					Jprows,Jpcols,numjunpairs,
					Sprows,Spcols,numsenpairs,
					comprows,compcols,
					Trv,Tcp,Srv,Scp;
					updatetransform = storetransform,
					verbose = verbose
				)
	maxnz = max(maxnz,Mcp[Mn[1]+1])

	if verbose
		println("first shurr finished")
		if Mn[1]>0
			println([Mcp[Mn[1]] "nnz(M)" Mm[1] "Mm" Mn[1] "Mn" length(Mrv) "length(Mrv)"])
		else
			println("Mn = 0")
		end
	end
	#gc()
	rowfilt = Array{Tv}(length(comprows)); colfilt = Array{Tv}(length(compcols))
	counter = 0
	while Mcp[Mn[1]+1]>1
		if verbose
			println("starting block operation $(counter)")
		end
		counter+=1
		for i = 1:Mm[1]
			rowfilt[i] = Mrowgrain[rowlab[i]]
		end
		for j = 1:Mn[1]
			colfilt[j] = Mcolgrain[collab[j]]
		end
		getPairsLightWrite2!(Mrv,Mcp,rowfilt,colfilt,Mm[1],Mn[1],Jprows,Jpcols,numjunpairs,verbose=verbose)
		comprows = intervalcomplementuniqueunsortedinput(Jprows[1:numjunpairs[1]],Mm[1])
		compcols = intervalcomplementuniqueunsortedinput(Jpcols[1:numjunpairs[1]],Mn[1])
		if diagnostic
			numsenpairsOLD = numsenpairs[1]
		end

		schurit4!(		Mrv,Mcp,Mm,Mn,Mn0,
						rowlab,collab,
						Jprows,Jpcols,numjunpairs,
						Sprows,Spcols,numsenpairs,
						comprows,compcols,
						Trv,Tcp,Srv,Scp;
						updatetransform = storetransform,
						verbose = verbose
					)
		maxnz = max(maxnz,Mcp[Mn[1]+1])
	end
	lastSrowmarker = Scp[numsenpairs[1]+1]-1
	lastTrowmarker = Tcp[Mn[1]+1]-1
	deleteat!(Srv,(lastSrowmarker+1):length(Srv))
	deleteat!(Trv,(lastTrowmarker+1):length(Trv))
	deleteat!(Scp,(numsenpairs[1]+1):length(Scp))
	deleteat!(Tcp,(Mn[1]+2):length(Tcp))
	deleteat!(Sprows,(numsenpairs[1]+1):maxnumpairs)
	deleteat!(Spcols,(numsenpairs[1]+1):maxnumpairs)
	Tcp+=lastSrowmarker
	append!(Scp,Tcp)
	append!(Srv,Trv[1:lastTrowmarker])
  	tlab = Spcols[1:numsenpairs[1]]
	append!(tlab,collab[1:Mn[1]])
	return Srv,Scp,Sprows,Spcols,tlab,maxnz
end

function persistf2_core_cell(
	Nrv,
	Ncp,
	grain;
	maxsd = length(Nrv),
	record="cyclerep",
	verbose=false,
	prepairs = fill(Array{Int64}(0),maxsd+1)
	)
	if record == "all" || record == "cyclerep"
		storetransform = true
	else
		storetransform = false
	end
	tcp::Array{Array{Int64,1},1}			=Array{Array{Int64,1},1}(maxsd+1)
	trv::Array{Array{Int64,1},1}			=Array{Array{Int64,1},1}(maxsd+1)
	phi::Array{Array{Int64,1},1}			=Array{Array{Int64,1},1}(maxsd+1)
	plo::Array{Array{Int64,1},1}			=Array{Array{Int64,1},1}(maxsd+1)
	tid::Array{Array{Int64,1},1}			=Array{Array{Int64,1},1}(maxsd+1)
	for i in [1,maxsd+1]
		tcp[i] 				=[1]
		trv[i]				=Array{Int64}(0)
		phi[i]				=Array{Int64}(0)
		plo[i]				=Array{Int64}(0)
		tid[i]				=Array{Int64}(0)
	end
	maxnzs 									=zeros(Int64,maxsd+1)
	m 										=length(Ncp[1])-1
	for sd = 2:maxsd
		if sd > length(Nrv)
			trv[sd] = Array{Int64}(0)
			tcp[sd] = ones(Int64,1)
			tid[sd] = Array{Int64}(0)
			plo[sd] = Array{Int64}(0)
			phi[sd] = Array{Int64}(0)
			continue
		elseif sd>2
			lowbasisnames = phi[sd-1]
		else
			lowbasisnames = Array{Int64}(0)
		end
		Mrv 			= Nrv[sd] # temporary
		Mcp 			= Ncp[sd] # temporary
		if isempty(Mrv)
			plo[sd] 	= Array{Int64}(0)
			phi[sd] 	= Array{Int64}(0)
			tid[sd] 	= Array{Int64}(1:numcols(Ncp[sd-1]))
						  deleteat!(tid[sd],sort(phi[sd-1]))
		    perm 		= sortperm(grain[sd-1][tid[sd]],alg=MergeSort)
			tid[sd] 	= tid[sd][perm]
			trv[sd] 	= Array{Int64}(0)
			tcp[sd] 	= ones(Int64,1+length(tid[sd]))
			continue
		end
		Mm0				= maximum(Mrv) # only used for the definition of the low-translator, which I believe is only there for debugging purposes
		Mm				= [maximum(Mrv)]  # temporary
		Mn 				= [length(Mcp)-1] # temporary
		higlab			= convert(Array{Int64},1:Mn[1])	# we'll assume prepairs to be empty, for now
		lowlab	 		= intervalcomplementuniqueunsortedinput(lowbasisnames,Mm[1])	# temporary, and we'll assume prepairs is empty for now
		nporder			= sortperm(grain[sd-1][lowlab],alg=MergeSort)
		lowlab			= lowlab[nporder]
		Mrv,Mcp			= transposeLighter_submatrix(
							Mrv, 					# the Arv argument
							Mcp, 					# the Acp argument
							maximum(Mrv), 			# the Am argument
							rows = lowlab,			# the rows selected
							cols = higlab)			# the columns selected
		Mn[1]			= length(Mcp)-1
		Mm[1] 			= maximum(Mrv)
 		lowlabtemp 		= convert(Array{Int64,1},1:length(lowlab))
 		higlabtemp 		= convert(Array{Int64,1},1:length(higlab))
 		higfilttemp 	= grain[sd][higlab]
 		lowfilttemp 	= grain[sd-1][lowlab]
		pplow 			= convert(Array,length(prepairs[sd]):-1:1) # we'll assume prepairs is empty, for now
		pphig 			= convert(Array,length(prepairs[sd]):-1:1) # we'll assume prepairs is empty, for now
		if verbose
			println("Constructed Morse boundary operator, columns indexed by cells of dimension $(sd-1)")
		end
		# NB: It is critical that the columns of the input array should
		# be ordered according to filtration; in particular, the entries of
		# lowfiltemp should increase monotonically
		Srv,Scp,Sphigs,Splows,tlab,maxnz =
		morselu!(
			Mrv,
			higfilttemp,
			Mcp,
			lowfilttemp,
			lowlabtemp,
			higlabtemp,
			pplow,
			pphig,
			Mm[1],
			storetransform = storetransform,
			verbose = verbose)
		trv[sd] 		= lowlab[Srv]
		tcp[sd] 		= Scp
		tid[sd] 		= lowlab[tlab]
		plo[sd] 		= lowlab[Splows]
		phi[sd] 		= higlab[Sphigs]
		maxnzs[sd] 	= maxnz
	end
	return trv,tcp,plo,phi,tid,maxnzs
end

function persistf2_core_vr(
	farfaces::Array{Array{Int64,1},1},
	firstv::Array{Array{Int64,1},1},
	prepairs::Array{Array{Int64,1},1},
	grain::Array{Array{Int64,1},1},
	maxsd::Integer;
	record="cyclerep",
	verbose=false)

	farfaces::Array{Array{Int64,1},1}
	firstv::Array{Array{Int64,1},1}
	prepairs::Array{Array{Int64,1},1}
	grain::Array{Array{Int64,1},1}

	if record == "all" || record == "cyclerep"
		storetransform = true
	else
		storetransform = false
	end

	m = length(firstv[1])-1

	trv::Array{Array{Int64,1},1}			=Array{Array{Int64,1},1}(maxsd+1);
	tcp::Array{Array{Int64,1},1}			=Array{Array{Int64,1},1}(maxsd+1);
	phi::Array{Array{Int64,1},1}			=Array{Array{Int64,1},1}(maxsd+1);
	plo::Array{Array{Int64,1},1}			=Array{Array{Int64,1},1}(maxsd+1);
	tid::Array{Array{Int64,1},1}			=Array{Array{Int64,1},1}(maxsd+1);
	for i in [1,maxsd+1]
		tcp[i]  			=[1];
		trv[i]				=Array{Int64}(0)
		tid[i]				=Array{Int64}(0)
		phi[i]				=Array{Int64}(0)
		plo[i]				=Array{Int64}(0)
	end

	maxnzs 									=zeros(Int64,maxsd+1);

	for sd = 2:maxsd
		if sd > length(farfaces)
			continue
		elseif sd>2
			lowbasisnames = phi[sd-1]
		else
			lowbasisnames = Array{Int64}(0)
		end
		Mrv::Array{Int64,1},
		Mcp::Array{Int64,1},
		lowlab::Array{Int64,1},
		higlab::Array{Int64,1},
		Mm =
		filteredmatrixfromfarfaces(farfaces,firstv,prepairs,grain,sd,lowbasisnames;verbose=verbose)
 		lowlabtemp = convert(Array{Int64,1},1:length(lowlab))
 		higlabtemp = convert(Array{Int64,1},1:length(higlab))
 		higfilttemp = grain[sd][higlab]
 		lowfilttemp = grain[sd-1][lowlab]
		if verbose
			println("Constructed Morse boundary operator, columns indexed by cells of dimension $(sd-1)")
		end
		pplow = convert(Array,length(prepairs[sd]):-1:1)
		pphig = convert(Array,length(prepairs[sd]):-1:1)

		# NB: It is critical that the columns of the input array
		# be ordered according to filtration; in particular, the entries of
		# lowfiltemp should increase monotonically
		Srv,Scp,Sphigs,Splows,tlab,maxnz =
		morselu!(
			Mrv,
			higfilttemp,
			Mcp,
			lowfilttemp,
			lowlabtemp,
			higlabtemp,
			pplow,
			pphig,
			Mm,
			storetransform = storetransform,
			verbose = verbose)
		trv[sd] = lowlab[Srv]
		tcp[sd] = Scp
		tid[sd] = lowlab[tlab]
		plo[sd] = lowlab[Splows]
		phi[sd] = higlab[Sphigs]
		maxnzs[sd]= maxnz
	end
	return trv,tcp,plo,phi,tid,maxnzs
end

function persistf2!(
	D::Dict;maxsd=0,
	dictionaryoutput::Bool = true,
	verbose::Bool = false,
	record = "cyclerep")

	farfaces = D["farfaces"]
	firstv = D["firstv"]
	prepairs = D["prepairs"]
	grain = D["grain"]
	if maxsd == 0
		maxsd = length(farfaces)-1
	end

	trv,tcp,plo,phi,tid,maxnzs =
	persistf2_core_vr(farfaces,firstv,prepairs,grain,maxsd::Integer;record=record,verbose = verbose)
	if dictionaryoutput == true
		D["trv"] = trv
		D["tcp"] = tcp
		D["tid"] = tid
		D["plo"] = plo
		D["phi"] = phi
		D["maxnz"] = maxnzs
		return D
	else
		return trv,tcp,plo,phi,tid
	end
end

function persistf2vr(
	s,
	maxsd;
	model 			= "vr",
	entryformat 		= "textfile",
	minrad			= -Inf,
	maxrad			= Inf,
	numrad			= Inf,
	nodrad  		= [],
	filfun 			= "n/a", # stands for filtration function; only to be used with fastop=false; function must take finite values; not compatible with minrad, maxrad, or numrad
	fastop			= true,
	vscale			= "diagonal",
	pointlabels 	= [],
	verbose 		= false,
	record 			= "cyclerep")

	#### Start timer
# 	tic()

	#### Must not stop early if an explicit filtration function is passed; also do no rounding
	if filfun 	   != 	"n/a"
		fastop 		= 	false
		numrad 		= 	1
		println("must define <clearprepairs>")
		return
	end

	#### Extract data as necessary
	inputisfile = false
	if typeof(s) == String
		inputisfile = true
		filename = s #modified 12/29/2017; note that strings are immutable
		if entryformat == "textfile"
			if typeof(readdlm(filename,','))<:Array{Float64}
				s = readdlm(s,',')
			elseif typeof(readdlm(filename,' '))<:Array{Float64}
				s = readdlm(s,' ')
			else
				print("Error reading text file.  Input files of this format must be either comma or space delimited.")
				return
			end
		elseif entryformat == "perseus_distmat"
			s,minrad,maxrad,maxsd = parse_perseusdistmat(filename)
			model = "clique_perseus"
			numrad = Inf
		elseif entryformat == "perseus_brips"
			s,minrad,maxrad = parse_perseusbrips(filename)
			s = s'  # the perseus brips format stores points as rows
			model = "pointcloud_perseusbrips"
			numrad = Inf
		end
	else
		filename = "user input julia array"
	end
	if model == "pc" || model == "perseus_brips"
		pc 	= 	"genera"
	else
		pc 	= 	"na"
	end

	#### Store the input
	input = Dict(
		"model"			=> model,
		"genera"		=> copy(s),
		"pc"			=> pc,
		"source"		=> filename,
		"maxdim" 		=> maxsd-2,
		"maxrad"		=> maxrad,
		"minrad"		=> minrad,
		"numrad"		=> numrad,
		"nodrad"		=> nodrad,
		"fastop"		=> fastop,
		"record"	 	=> record,
		"filfun" 		=> filfun,
		"version" 		=> Pkg.installed("Eirene"),
		"date"			=> string(Dates.Date(now())),
		"time"			=> string(Dates.Time(now()))
		)

	#### Determine the number of points
	if model == "pc"
		numpoints = size(s,2)
	elseif model == "vr"
		numpoints = size(s,1)
		if !issymmetric(s)
			print("It appears the input matrix is not symmetric.  Only symmetric distance matrices are accepted when the <model> keyword argument has value \"clique\".")
			return
		end
	else
		print("Keyword argument <model> must be \"vr\", \"pc\", or \"cell\".  Please see documentaiton for further details.")
	end

	#### Extract labels
	if pointlabels == "none" || pointlabels == []
		pointlabels = 1:numpoints
	elseif pointlabels in ["left","right","top","bottom"]
		(s,pointlabels) = separatelabels(s,pointlabels)
	elseif typeof(pointlabels) == String
		input["source_pointlabels"] = pointlabels
		pointlabels = ezread(pointlabels)
	end
	if length(pointlabels) != numpoints
		warn("It appears the number of vertex labels does not match the number of vertices.")
	end
	pointlabels = ezlabel(pointlabels)
	for i = 1:length(pointlabels)
		pointlabels[i] = "$(pointlabels[i])"
	end
	input["pointlabels"] = pointlabels

	#### type the matrix
	s = convert(Array{Float64,2},s)
	if model == "pc"
		d = Distances.pairwise(Euclidean(),s)
		if !isempty(nodrad)
			for 	i 	= 	1:numpoints
				d[i,i] 	= 	nodrad[i]
			end
		end
	elseif model == "vr"
		d = convert(Array{Float64,2},s)
	end

	################################################################################

	if fastop
		maxrad_alt 	= 	minimum(maximum(d,1))
		maxrad_alt  = 	min(maxrad_alt,maxrad)
	else
		maxrad_alt 	= 	maxrad
	end

	d 			= 	minmaxceil(		d,
									minrad 	= 	minrad,
									maxrad 	= 	maxrad_alt,
									numrad 	= 	numrad)

	# vfilt 		= 	diag(t)
	# recall that t will have the same order (NOT inverted) as d
	# <trueordercanonicalform> is a bit like <integersinsameorder>, just valid for floating point inputs, and with a bit more data in the output
	t,ocg2rad 	= 	trueordercanonicalform(d,factor=true)

	t 			= 	1+maximum(t)-t
	ocg2rad 	= 	flipdim(ocg2rad,1)

	if 	any(d.>maxrad)
		t 		= 	t-1
		deleteat!(ocg2rad,1)
	end

	vertices2keep 	= 	find(diag(t).!=0)  # this step is necessary in order to cover the case where some vertices never enter the filtration
	t 				= 	t[vertices2keep,vertices2keep]

	#### Build the complex
	D = buildcomplex3(t,maxsd;verbose = verbose)
	D["ocg2rad"]=ocg2rad

	################################################################################

	################################################################################
	# arbitrary function values
	# NB assumes the function takes only finite values
	if filfun != "n/a"

		fv 	= 	Array{Array{Float64,1}}(maxsd) # fv stands for filtration values
		for sd = 1:maxsd
			fv[sd] 	= 	Array{Float64}(length(D["farfaces"][sd]))
			for p 	= 	1:length(D["grain"][sd])
				# syntax reminder: vertexrealization(D::Dict,facecardinality,facenames)
				fv[sd][p] 	= 	filfun(vertexrealization(D,sd,[p]))
			end
		end

		ocg,ocg2rad 	= 	trueordercanonicalform(cat(1,fv...),factor=true) # ocg stands for order canonical grain
		ocg 			= 	1+maximum(ocg)-ocg
		ocg2rad 		= 	flipdim(ocg2rad,1)


		D["ocg2rad"]	=	ocg2rad
		D["grain"] 		= 	ocg
	end

	################################################################################

	#### Compute persistence
	persistf2!(D;verbose = verbose,record = record)

	#### Store input data
	D["input"] 		= 	input
	D["nvl2ovl"]	= 	vertices2keep[D["nvl2ovl"]]  # this covers the case where some vertices never enter the filtration

	#### Store generators
	#gc()
	if record == "all" || record == "cyclerep"
 		unpack!(D)
	end
	#gc()
	if record == "cyclerep"
		delete!(D,"trv")
		delete!(D,"tcp")
		delete!(D,"Lrv")
		delete!(D,"Lcp")
		delete!(D,"Rrv")
		delete!(D,"Rcp")
		delete!(D,"Lirv")
		delete!(D,"Licp")
		delete!(D,"prepairs")
	end

	#### Record time
# 	D["input"]["computationtime"] = toc()

	return D
end

#=
- 	drafted 12/27/2017
- 	version with 1 (string) non-keyword argument
- 	formatting guidelines:
	- the input file should be a csv in which
	- line 1: cell dimensions
	- line 2: cell filtrations
	- line 3: boundary matrix row values
	- line 4: boundary matrix column pattern

NB: the default value for maxdim has not been tested, and may cause errors;
the -3 accounts for (1) julia uses 1-indexed arrays, (2) to calculate
homology in dimension p, one must inspect the (p+1)-dimensional boundary
operator, (3) this operator should be givent the same treatment as those
that precede it ... generally this assumes that the next one up is at least
defined, even if it is trivial.
=#
function persistf2complex(filepath::String;
						maxdim=Inf,
						entryformat = "sp",
						record = "cyclerep",
						verbose=false)

	if entryformat 		== 	"sp"
		dp,fv,rv,cp 	= 	humanreadablefilepath2unsegmentedfilteredcomplex(filepath)
	elseif in(entryformat,["dp","dv","ev"])
		dp,fv,rv,cp 	= 	cscfilepath2unsegmentedfilteredcomplex(filepath,toprow = entryformat)
	end

	if 	maxdim 		== 	Inf
		maxdim 		= 	length(dp)-3
	end

	C  				= 	persistf2complex(
						rv = rv,
						cp = cp,
						fv = fv,
						dp = dp,
						maxdim=maxdim,
						record=record,
						verbose=verbose)
	C["input"]["source"]		= 	filepath;
	C["input"]["entryformat"]	= 	entryformat;
	return C
end

function emptyunsegmentedfilteredcomplex_dp()
	rv 	= zeros(Int64,0)
	cp 	= ones(Int64,1)
	fv	= zeros(Float64,0)
	dp  = ones(Int64,1)
	return dp,fv,rv,cp #rv,cp,fv,dp
end

# possible values for toprow: dp, dv, ev
function cscfilepath2unsegmentedfilteredcomplex(filepath;toprow="dp")
	M = readcsv(filepath)
	endpoints = csvimport2linends(M)

	if size(M,1) 	== 1
		if in(toprow,["dv","ev"]) && M == ones(Float64,1,1)
			return emptyunsegmentedfilteredcomplex_dp() # this returns values for rv,cp,fv,dp
		else
			print("Error: please check formatting of input file.")
			return
		end
	end

	if size(M,1) 	== 2
		if in(toprow,["dp"]) && M == ones(Float64,2,1)
			return emptyunsegmentedfilteredcomplex_dp() # this returns values for rv,cp,fv,dp
		else
			print("Error: input file should have .csv format with four lines.")
			return
		end
	end

	if size(M,1) 	== 3
		if any(M[3,:].!=1)
			print("Error: please check formatting of input file.")
			return
		else
			xx 	    = M[1,1:endpoints[1]]  # stands for dimension pattern
			fv      = M[2,1:endpoints[2]]  # stands for filtration values
			rv 	    = zeros(Int64,0)       # stands for row values
			cp      = M[4,1:endpoints[4]]  # stands for column pattern
		end
	end

	if size(M,1)	== 4
		xx 	    	= M[1,1:endpoints[1]]  # stands for dimension pattern
		fv      	= M[2,1:endpoints[2]]  # stands for filtration values
		rv 	    	= M[3,1:endpoints[3]]  # stands for row values
		cp      	= M[4,1:endpoints[4]]  # stands for column pattern
	end

	if toprow 			== 	"dp"
		dp 				= 	xx
	elseif toprow	 	== 	"dv"
		dp 				= 	dimensionvalues2dimensionpattern(xx)
	elseif toprow	 	== 	"ev"
		dp 				= 	eulervector2dimensionpattern(xx)
	end

	dp 		= convert(Array{Int64,1},dp)
	fv      = convert(Array{Float64,1},fv)
	rv 		= convert(Array{Int64,1},rv)
	cp 		= convert(Array{Int64,1},cp)

	return dp,fv,rv,cp
end

function humanreadablefilepath2unsegmentedfilteredcomplex(filepath)
	M 					= 	readcsv(filepath)
	m 					= 	size(M,1)
	endpoints 			= 	csvimport2linends(M)

	dv 					= 	Array{Int64}(M[:,1])
	fv 					= 	Array{Float64}(M[:,2])
	dp 					= 	dimensionvalues2dimensionpattern(dv)

	nrv 				= 	sum(endpoints)-2*m
	rv 					= 	zeros(Int64,nrv)
	cp 					= 	zeros(Int64,m+1)
	cp[1]				= 	1

	rvpost 				= 	0
	cppost 				= 	1
	for p 				= 	1:m
		numval 			= 	endpoints[p]-2
		readran 		= 	3:endpoints[p]
		fillran_a 		= 	rvpost+1
		fillran_b 		= 	rvpost+numval
		fillran 		= 	fillran_a:fillran_b

		vals 			= 	Array{Int64,1}(M[p,readran])
		rv[fillran] 	= 	vals
		cp[cppost+1] 	= 	cp[cppost]+numval

		rvpost 			= 	fillran_b
		cppost 			= 	cppost+1
	end

	return dp,fv,rv,cp
end

function segmentarray{Tv}(vr::Tv,vp)
	m 	= 	length(vp)-1
	u 	= 	Array{Tv}(m)
	for 	p 	=	1:m
		u[p] 	= 	crows(vp,vr,p)
	end
	return u
end

function checksegmentarray(numits)
	for p 			= 	1:numits
		cp 			= 	rand(1:100,50)
		dp 			= 	eulervector2dimensionpattern(cp)
		rv 			= 	rand(dp[end]-1)
		A 			= 	segmentarray(rv,dp)
		for q 		= 	1:50
			if 	A[q]	!=	crows(dp,rv,q)
				return rv,cp
			end
		end
	end
	return []
end

function unsegmentedfilteredcomplex2segmentedfilteredcomplex(rv,cp,fv,dp;ncd=Inf)
	# ncd stands for number of chain dimensions
	# nsd stands for number of stored dimensions
	nsd 	= 	length(dp)-1
	if ncd == Inf
		ncd  	= nsd
	end
	m 		= 	min(nsd,ncd)

	fvc     = Array{Array{Float64,1}}(ncd)
	rvc     = Array{Array{Int64,1}}(ncd)
	cpc     = Array{Array{Int64,1}}(ncd)

	# remarks:
	# (a) #{cells of dimension ≤ (p   = k+1)} 				= dp[k+2]-1		= 	dp[p+1]-1
	# (b) #{cells of dimension ≤ (p-2 = (k-2)+1 = k-1)} 	= dp[k]-1 		=	dp[p-1]-1
	for p = 1:m
		rvc[p],cpc[p]   = 	copycolumnsubmatrix(rv,cp,cran(dp,p))
		rvc[p]			= 	rvc[p] - ec(dp,p-1,1) + 1  # we subtract off the number of cells of dimension 2 less than those of interest, since starts the _faces_ of the cells of interest off at index 1; we get away with this in dimension 0 bc rvc[1] is the empty set
		fvc[p] 			=   convert(Array{Float64,1},fv[cran(dp,p)])
	end

	for p = (m+1):(ncd)
		rvc[p]		=	zeros(Int64,0)
		cpc[p]		=	ones(Int64,1)
		fvc[p] 		=	zeros(Int64,0)
	end

	dpc 			= vcat(dp,fill(dp[end],ncd+1-length(dp))) # extend dp to the proper length
	return rvc,cpc,fvc,dpc
end

function segmentedfilteredcomplex2unsegmentedfilteredcomplex(rv,cp,fv)
	numsd 					 	= 	length(fv)
	eulervec 					= 	zeros(Int64,numsd)
	for p 						= 	1:numsd
		eulervec[p]				= 	length(cp[p])-1
	end
	dp 							= 	eulervector2dimensionpattern(eulervec)

	cp 							= 	copy(cp)
	for 	p 					= 	2:numsd
		cp[p]					= 	cp[p-1][end]-1+cp[p]
		deleteat!(cp[p],1)
	end
	cp 							= 	cat(1,cp...)

	rv 							= 	copy(rv)
	for 	p 					= 	2:numsd
		rv[p]					= 	rv[p] + ec(dp,p-1,1) - 1
	end
	rv 							= 	cat(1,rv...)

	fv 							= 	cat(1,fv...)
	return rv,cp,fv,dp
end

function checksegVdesegcomplex(numits)
	for p 			= 	1:numits
		x 			= 	rand(30,30)
		x 			= 	x+x'
		C 			= 	eirene(x,model="vr",maxdim=2)
		rv,cp 		= 	boundarymatrices(C)
		fv 			= 	ocff2of(C["grain"],C["ocg2rad"])
		rvold 		= 	copy(rv)


		rv1,cp1,fv1,dp1 	= 	segmentedfilteredcomplex2unsegmentedfilteredcomplex(rv,cp,fv)
		rv2,cp2,fv2,dp2 	= 	unsegmentedfilteredcomplex2segmentedfilteredcomplex(rv1,cp1,fv1,dp1)

		check0 	= 	dp1 == dp2
		check1 	= 	rv2 == rv
		check2 	= 	cp2 == cp
		check3 	= 	fv2 == fv
		check4 	= 	true
		for q 	= 	1:length(dp1)-1
			if 	length(fv[q]) != length(cran(dp1,q))
				check4 	= 	false
				break
			end
		end

		if !all([check0, check1, check2,check3,check4])
			println([check0, check1,check2,check3,check4])
			return x,rv1,cp1,fv1,dp1,rv2,cp2,fv2,rv,cp,fv,rvold
		end
	end
	return []
end

#=
version with 4 (integer array) non-keyword arguments
- 	rv:	row values
- 	cp: column pattern
- 	fv: filtration values
- 	dp: dimension pattern (satisfies cran(dp,j) = {cells of dimension j-1})

NB: the default value for maxdim has not been tested, and may cause errors;
the -3 accounts for (1) julia uses 1-indexed arrays, (2) to calculate
homology in dimension p, one must inspect the (p+1)-dimensional boundary
operator, (3) this operator should be givent the same treatment as those
that precede it ... generally this assumes that the next one up is at least
defined, even if it is trivial.
=#
function persistf2complex(	;
							rv 			= zeros(Int64,0),
							cp 			= zeros(Int64,0),
							fv 			= zeros(Int64,0),
							dp 			= zeros(Int64,0),
							dv 			= zeros(Int64,0),
							ev 			= zeros(Int64,0),
							maxdim		= [],
							numrad 		= Inf,
							maxrad 		= Inf,
							minrad 		= -Inf,
							record 		= "cyclerep",
							verbose		= false)

	if isempty(maxdim)
		if isempty(rv)
			ncd 			= 	1;
		elseif typeof(rv)  != 	Array{Int64,1}
			ncd 			= 	length(rv)-1 # we subtract one b/c it's usually convenient to have an extra empty array at the end of rv
		elseif !isempty(ev)
			ncd 			= 	length(ev)
		elseif !isempty(dv)
			ncd 			= 	dv[end]+1
			if !isempty(checkdv(rv,cp,dv))
				println("error: please check that the input operator is graded of degree 1")
				return
			# else
			# 	println("input operator is graded of degree 1")
			end
		elseif !isempty(dp)
			ncd 			= 	length(dp)-1
		else
			println("error: keyword input <rv> is a vector, and none of <dv>, <dp>, and <ev> is nonempty.")
		end
	else
		ncd 				= 	maxdim+2
	end


	if typeof(rv) == Array{Int64,1}
		if isempty(dp)
			if !isempty(ev)
				dp 				= 	eulervector2dimensionpattern(ev)
			end
			if !isempty(dv)
				dp 				= 	dimensionvalues2dimensionpattern(dv)
			end
		end
	end

	if !isempty(dp) # segment the complex if necessary
		rv,cp,fv,dp 	= 	unsegmentedfilteredcomplex2segmentedfilteredcomplex(rv,cp,fv,dp;ncd=ncd)
		# #############################################################################
		# if !pairwiseisequal([rv,cp,fv],under=length)
		# 	println("Error: please check unsegmentedfilteredcomplex2segmentedfilteredcomplex")
		# end
		# #############################################################################
	else # if it has not been defined by this point, define the dimension pattern
		ev 			= 	zeros(Int64,ncd)
		for p 			= 	1:ncd
			ev[p] 	= 	length(fv[p])
		end
		dp 				= 	eulervector2dimensionpattern(ev)
	end

	### Record the input parameters
	input = Dict(
		"model"			=> "complex",
		"version" 		=> Pkg.installed("Eirene"),
		"date"			=> string(Dates.Date(now())),
		"time"			=> string(Dates.Time(now())),
		"maxdim"		=> ncd-2,
		"record" 		=> record
		)

	if typeof(fv) == Array{Float64,1}
		# println("nonerror message: TYPE IS FLOATING ARRAY!")
		ocg,ocg2rad			=   trueordercanonicalform(
								fv,
								firstval=1,
								factor=true,
								rev=true)
	else
	 	ocg,ocg2rad			=   trueordercanonicalform(
								cat(1,fv...),
								firstval=1,
								factor=true,
								rev=true)
	end

	ocg 					= 	segmentarray(ocg,dp)

	### Perform the persistence computation
	trv,tcp,plo,phi,tid =
	persistf2_core_cell(
		rv,
		cp,
		ocg;
		maxsd = ncd,
		record=record,
		verbose=false,
		prepairs = convert(Array{Array{Int64,1},1},fill(Array{Int64}(0),ncd))
		)

	### Create the dictionary that stores all relevant data
	D = Dict(
		"rv" 		=> rv,
		"cp" 		=> cp,
		"grain"		=> ocg,
		"trv" 		=> trv,
		"tcp" 		=> tcp,
		"tid" 		=> tid,
		"plo" 		=> plo,
		"phi" 		=> phi,
		"ocg2rad" 	=> ocg2rad,
		"input"		=> input
		)

	#### Store generators
	#gc()
	if record == "all" || record == "cyclerep"
 		unpack!(D)
	end
	#gc()
	if record == "cyclerep"
		delete!(D,"trv")
		delete!(D,"tcp")
		delete!(D,"Lrv")
		delete!(D,"Lcp")
		delete!(D,"Rrv")
		delete!(D,"Rcp")
		delete!(D,"Lirv")
		delete!(D,"Licp")
		delete!(D,"prepairs")
	end

	return D
end


function parse_perseusdistmat(filename)
	A = readdlm(filename)
	s = convert(Array{Float64},A[3:end,:])
	minrad = convert(Float64,A[2,1])
	maxrad = minrad+convert(Float64,A[2,2])*convert(Float64,A[2,3])
	maxsd = convert(Int64,A[2,4])+1
	return s,minrad,maxrad,maxsd
end

function parse_perseusbrips(filename)
	A = readdlm(filename)
	s = convert(Array{Float64},A[3:end,1:end-1])
	minrad = 0
	maxrad = convert(Array{Float64},A[2,1:3])
	maxrad = prod(maxrad)
	return s,minrad,maxrad
end

function boundarylimit_Lionly(trv,tcp,tid,sd;verbose=false)
	if isempty(tid[sd])
		Lirv = zeros(Int64,0)
		Licp = ones(Int64,1)
	else
		tsize			= length(tcp[sd])-1
		Lrv 			= copy(trv[sd])
		lowtranslator 	= zeros(Int64,maximum(tid[sd]))
		lowtranslator[tid[sd]] = 1:length(tid[sd])
		yafterx!(lowtranslator,Lrv)
		Lirv,Licp   	= morseInverseF2orderedColsUnsortedRowsSilentInSilentOut(Lrv,tcp[sd])
		Lirv,Licp   	= transposeLighter(Lirv,Licp,tsize)
	end
	return Lirv,Licp
end

#=
-	L, Li, and R are all indexed by tid, just like (trv,tcp)
- 	L (not Li) is the transpose of (trv,tcp)
- 	up to permutation and deletion of some zero rows, RLM is the fundamental
 	cycle matrix of the d-dimensional boundary operator with respect
	to the doubly-minimal basis we have constructed, where M the submatrix of
	the boundary with rows indexed by tid.
=#
function boundarylimit_core(brv,bcp,trv,tcp,tid,numl,nump,numnlpl)
	lowtranslator = zeros(Int64,numl)
	lowtranslator[tid] = 1:numnlpl
	trv = copy(trv)
	yafterx!(lowtranslator,trv)
	#################################################################################
	if !issorted(tcp)
		println("error: tcp is not sorted")
		sleep(2)
	end
	if length(trv) != tcp[end]-1
		println("error: tcp appears to end in the wrong place")
		sleep(2)
	end
	if tcp[1] != 1
		println("error: tcp does not begin with 1")
		sleep(2)
	end
	#################################################################################
	Lirv,Licp	= morseInverseF2orderedColsUnsortedRowsSilentInSilentOut(trv,tcp)
	Lirv,Licp   = transposeLighter(Lirv,Licp,numnlpl)
	Lrv,Lcp 	= transposeLighter(trv,tcp,numnlpl)
	brv,bcp 	= spmmF2silentLeft(Lrv,Lcp,brv,bcp,numnlpl)
	Rrv,Rcp 	= morseInverseF2orderedColsUnsortedRowsInSilentOut(brv,bcp)
	return Lrv,Lcp,Lirv,Licp,Rrv,Rcp
end

function boundarylimit_simplex(farfaces,firstv,trv,tcp,plo,phi,tid,sd;verbose=false)
	numl 		= length(farfaces[sd-1])
	nump 		= length(phi[sd])
	numnlpl 	= numl-length(plo[sd-1])
	brv,bcp 	= maxnonsingblock_simplex(farfaces,firstv,plo,phi,tid,sd;verbose=false)
	return 	boundarylimit_core(brv,bcp,trv[sd],tcp[sd],tid[sd],numl,nump,numnlpl)
end

function boundarylimit_cell(rv,cp,trv,tcp,plo,phi,tid,sd;verbose=false)
	numl 		= length(cp[sd-1])-1
	nump 		= length(phi[sd])
	numnlpl 	= numl-length(phi[sd-1])
	brv,bcp		= maxnonsingblock_cell(rv,cp,plo,phi,tid,sd;verbose=false)
	Lrv,Lcp,Lirv,Licp,Rrv,Rcp = boundarylimit_core(brv,bcp,trv[sd],tcp[sd],tid[sd],numl,nump,numnlpl)
	return 	boundarylimit_core(brv,bcp,trv[sd],tcp[sd],tid[sd],numl,nump,numnlpl)
end

function boundarylimit(D::Dict,sd)
	# special case sd = 1 may possibly be degenerate
	trv = D["trv"];tcp=D["tcp"];plo=D["plo"];phi=D["phi"];tid=D["tid"]
	if haskey(D,"farfaces")
		farfaces = D["farfaces"];firstv = D["firstv"]
		return boundarylimit_simplex(farfaces,firstv,trv,tcp,plo,phi,tid,sd;verbose=false)
	else
		rv = D["rv"];cp = D["cp"]
		# Lrv,Lcp,Lirv,Licp,Rrv,Rcp = boundarylimit_cell(rv,cp,trv,tcp,plo,phi,tid,sd;verbose=false)
		return boundarylimit_cell(rv,cp,trv,tcp,plo,phi,tid,sd;verbose=false)
	end
end

function maxnonsingblock_simplex(farfaces,firstv,plo,phi,tid,sd;verbose=false)
	numl = length(farfaces[sd-1])
	nump = length(phi[sd])
	numnlpl = numl-length(plo[sd-1])

	lowtranslator = zeros(Int64,numl)

	lowtranslator[tid[sd]] = 1:numnlpl
	brv = ff2aflight(farfaces,firstv,sd,phi[sd])
	brv = reshape(brv,length(brv))
	bcp = convert(Array{Int64,1},1:sd:(nump*sd+1))
	supportedmatrix!(brv,bcp,tid[sd],1:nump,numl)
	yafterx!(lowtranslator,brv)
	append!(bcp,fill(bcp[end],numnlpl-nump)) # <-- note there's no +1 b/c bcp is already 1 elt. too long
	return brv,bcp
end

function maxnonsingblock_cell(rv,cp,plo,phi,tid,sd;verbose=false)
	numl = length(cp[sd-1])-1
	nump = length(phi[sd])
	numnlpl = numl-length(plo[sd-1])

	if isempty(phi[sd])
		brv 	= 	zeros(Int64,0)
		bcp 	=   ones(Int64,numnlpl+1)
		return 		brv,bcp
	end

	lowtranslator = zeros(Int64,numl)
	lowtranslator[tid[sd]] = 1:numnlpl
	dummy0 = zeros(Int64,0)

	brv,dummy1,bcp,dummy2 = stackedsubmatrices(rv[sd],cp[sd],tid[sd],dummy0,phi[sd],numl)
	yafterx!(lowtranslator,brv)
	append!(bcp,fill(bcp[end],numnlpl-nump)) # <-- note there's no +1 b/c bcp is already 1 elt. longer than the # of columns
	return brv,bcp
end

function unpack!(D::Dict)
	# l = length(D["grain"])
	maxsd = D["input"]["maxdim"]+2 # grain2maxsd(D["grain"])

	Lirv = Array{Array{Int64,1},1}(maxsd);  Licp = Array{Array{Int64,1},1}(maxsd)
	Lrv  = Array{Array{Int64,1},1}(maxsd);  Lcp  = Array{Array{Int64,1},1}(maxsd)
	Rrv  = Array{Array{Int64,1},1}(maxsd);  Rcp  = Array{Array{Int64,1},1}(maxsd)

	if 	D["input"]["record"] == "all"
		N 	= 	maxsd
	else
		N 	= 	maxsd-1
		if maxsd >= 2
			Lirv[maxsd],Licp[maxsd] = boundarylimit_Lionly(D["trv"],D["tcp"],D["tid"],maxsd)
		elseif maxsd == 1
			Lirv[maxsd] = Array{Int64,1}(0)
			Licp[maxsd] = ones(Int64,1)
		end
		Lrv[maxsd]=Array{Int64,1}(0)
		Lcp[maxsd]=Array{Int64,1}(0)
		Rrv[maxsd]=Array{Int64,1}(0)
		Rcp[maxsd]=Array{Int64,1}(0)
	end

	for i = 2:N
		Lrv[i],Lcp[i],Lirv[i],Licp[i],Rrv[i],Rcp[i] = boundarylimit(D,i)
		if isempty(Lcp[i])
			println()
			println("ERROR MESSAGE IN unpack!: Lcp[i] = 0 and i = $(i)")
		end
	end

	Lirv[1]		=Array{Int64}(0)
	Lrv[1]		=Array{Int64}(0);
	Rrv[1]		=Array{Int64}(0)
	Licp[1]		=ones(Int64,1)
	Lcp[1]		=ones(Int64,1)
	Rcp[1]		=ones(Int64,1)

	D["Lrv"] = Lrv
	D["Lcp"] = Lcp
	D["Lirv"]= Lirv
	D["Licp"]= Licp
	D["Rrv"] = Rrv
	D["Rcp"] = Rcp

	D["cyclerep"] = fill(Array{Array{Int64,1},1}(0),maxsd)  # for each dimension one stores an array of arrays

	for i = 2:maxsd
		dim = i-2
		m = nnzbars(D,dim=dim)
		cyclenames = barname2cyclename(D,1:m,dim=dim)
		D["cyclerep"][i] = getcycle(D,cyclenames,dim=dim)
	end
	return
end

##########################################################################################

####	INVERSION

##########################################################################################

function morseInverseF2orderedColsUnsortedRowsInSilentOut{Tv<:Integer}(Arowval::Array{Tv,1},Acolptr::Array{Tv,1})
	mA = length(Acolptr)-1
	const colptrA = Acolptr
	const rowvalA = Arowval
    const preallocationIncrement = colptrA[end]

    colptrC = Array{Tv}(mA+1); colptrC[1]=1
	rowSupp = zeros(Tv, mA)
	rowList = Array{Tv}( mA)
	rowvalCj = Array{Bool}( mA)
	rowvalC = Array{Tv}( mA)
    totalrowscounter = 0
    onepast = 0
	for i in 1:mA
		if colptrC[i]+mA > length(rowvalC)+1
			append!(rowvalC, Array{Int64}(preallocationIncrement))
		end
		if colptrA[i]+1 == colptrA[i+1]
			colptrC[i+1] = colptrC[i]
		elseif colptrA[i]+2==colptrA[i+1]
			if Arowval[colptrA[i]]<i
				k = Arowval[colptrA[i]]
			else
				k = Arowval[colptrA[i]+1]
			end
			ccap = colptrC[i]+colptrC[k+1]-colptrC[k]-1
			rowvalC[colptrC[i]:ccap]= rowvalC[colptrC[k]:(colptrC[k+1]-1)]
			rowvalC[ccap+1]=k
			colptrC[i+1]=ccap+2
		else
			eyerange = cran(Acolptr,i)
			newRowsCounter = length(eyerange)
			for l=1:newRowsCounter
				ll = Arowval[eyerange[l]]
				rowList[l] = ll
				rowSupp[ll] = i
				rowvalCj[ll]=true
			end
			rowvalCj[i] = false ## note we have to make this correction
			for jp in eyerange
				j = rowvalA[jp]
				if j < i
					for kp in colptrC[j]:(colptrC[j+1] - 1)
						k = rowvalC[kp]
						if rowSupp[k] != i
							rowSupp[k] = i
							newRowsCounter +=1
							rowList[newRowsCounter] = k
							rowvalCj[k] = true
						else
							rowvalCj[k] = !rowvalCj[k]
						end
					end
				end
			end
			marker = colptrC[i]
			for l = 1:newRowsCounter
				if rowvalCj[rowList[l]]
					rowvalC[marker]=rowList[l]
					marker+=1
				end
			end
			colptrC[i+1]=marker
		end

	end
	deleteat!(rowvalC,colptrC[end]:length(rowvalC))
	return rowvalC, colptrC
end

function morseInverseF2orderedColsUnsortedRowsSilentInSilentOut{Tv<:Integer}(Arowval::Array{Tv,1},Acolptr::Array{Tv,1})
	mA = length(Acolptr)-1
	const colptrA = Acolptr
	const rowvalA = Arowval
    const preallocationIncrement = colptrA[end]

    colptrC = Array{Tv}(mA+1); colptrC[1]=1
	rowSupp = zeros(Tv, mA)
	rowList = Array{Tv}( mA)
	rowvalCj = Array{Bool}( mA)
	rowvalC = Array{Tv}( mA)
    totalrowscounter = 0
    onepast = 0
	for i in 1:mA
		if colptrC[i]+mA > length(rowvalC)+1
			append!(rowvalC, Array{Int64}(preallocationIncrement))
		end
		if colptrA[i] == colptrA[i+1]
			colptrC[i+1] = colptrC[i]
		elseif colptrA[i]+1==colptrA[i+1]
			k = Arowval[colptrA[i]]
			ccap = colptrC[i]+colptrC[k+1]-colptrC[k]
			rowvalC[colptrC[i]:(ccap-1)]= crows(colptrC,rowvalC,k)
			rowvalC[ccap]=k
			colptrC[i+1]=ccap+1
		else
			eyerange = cran(Acolptr,i)
			newRowsCounter = length(eyerange)
			for l=1:newRowsCounter
				ll = Arowval[eyerange[l]]
				rowList[l] = ll
				rowSupp[ll] = i
				rowvalCj[ll]=true
			end
			for jp in eyerange
				for kp in cran(colptrC,rowvalA[jp])
					k = rowvalC[kp]
					if rowSupp[k] != i
						rowSupp[k] = i
						newRowsCounter +=1
						rowList[newRowsCounter] = k
						rowvalCj[k] = true
					else
						rowvalCj[k] = !rowvalCj[k]
					end
				end
			end
			marker = colptrC[i]
			for l = 1:newRowsCounter
				if rowvalCj[rowList[l]]
					rowvalC[marker]=rowList[l]
					marker+=1
				end
			end
			colptrC[i+1]=marker
		end
	end
	deleteat!(rowvalC,colptrC[end]:length(rowvalC))
	return rowvalC, colptrC
end

##########################################################################################

####	MULTIPLICATION

##########################################################################################

function addcol!(
	oddfloods::BitArray{1},shoreline::Array{Int64,1},watermark::Int64,
	peakcounter::Array{Int64,1},Acp::Array{Int64,1},Arv::Array{Int64,1},
	flippedlist::Array{Int64,1},j::Int64)
	for i = cran(Acp,j)
		ii = Arv[i]
		if shoreline[ii] != watermark
			peakcounter[1]+=1
			flippedlist[peakcounter]=ii
			shoreline[ii] = watermark
			oddfloods[ii] = true
		else
			oddfloods[ii] = !oddfloods[ii]
		end
	end
end

function spmmF2{Tv<:Integer}(Arowval::Array{Tv,1},Acolptr::Array{Tv,1},Browval::Array{Tv,1},Bcolptr::Array{Tv,1},Am)
    const mA = Am
    const nB = length(Bcolptr)-1
    const rowvalA = Arowval; colptrA = Acolptr
    const rowvalB = Browval; colptrB = Bcolptr
    const preallocationIncrement = colptrA[end]+colptrB[end]

	colptrC = Array{Tv}( nB+1)
    colptrC[1] = 1
	rowSupp = zeros(Tv, mA)
	rowList = Array{Tv}( mA)
	rowvalCj = Array{Bool}( mA)
	rowvalC = Array{Tv}( preallocationIncrement)
	for i in 1:nB
		newrowscounter = 0
		for jp in colptrB[i]:(colptrB[i+1] - 1)
			j = rowvalB[jp]
			for kp in colptrA[j]:(colptrA[j+1] - 1)
				k = rowvalA[kp]
				if rowSupp[k] != i
					rowSupp[k] = i
					newrowscounter +=1
					rowList[newrowscounter] = k
					rowvalCj[k] = true
				else
					rowvalCj[k] = !rowvalCj[k]
				end
			end
		end
		nzRows = find(rowvalCj[rowList[1:newrowscounter]])
		colptrC[i+1] = colptrC[i]+length(nzRows)

		if colptrC[i+1] > length(rowvalC)+1
			append!(rowvalC, Array{Int}(preallocationIncrement))
		end
		rowvalC[colptrC[i]:(colptrC[i+1]-1)] = sort(rowList[nzRows])
	end
	deleteat!(rowvalC,colptrC[end]:length(rowvalC))
	return rowvalC, colptrC
end

function spmmF2_testfun()
	N 			=	1000; 	# matrix dimension
	n 			= 	100;  	# number of samples
	for p = 1:n
		m1 	=	sprand(N,N,0.01); 	rv1 =	m1.rowval;  	cp1 =	m1.colptr;
		m2 	=	sprand(N,N,0.01);	rv2 =	m2.rowval;		cp2 =	m2.colptr;
		if !isempty(rv1)
			m1m2 		= ceil.(Int64,m1)*ceil.(Int64,m2)
			m1m2 		= mod.(m1m2,2)
			qrv 		= m1m2.rowval;
			qcp 		= m1m2.colptr;
			prv,pcp 	= spmmF2(rv1,cp1,rv2,cp2,maximum(rv1))
			for q 	= 	1:N
				if crows(qcp,qrv,q) != sort(crows(pcp,prv,q))
					print("error 1: please see spmmF2_testfun")
					return
				end
			end
		end
		prv,pcp 	= 	spmmF2(zeros(Int64,0),ones(Int64,N+1),rv2,cp2,1)
		if (prv 	!=	zeros(Int64,0)) || (pcp != ones(Int64,N+1))
			print("error 2: please see spmmF2_testfun")
			return
		end
	end
	print("test complete - no errors detected")
end

function spmmF2silentLeft{Tv<:Integer}(Arowval::Array{Tv,1},Acolptr::Array{Tv,1},Browval::Array{Tv,1},Bcolptr::Array{Tv,1},Am)
    const mA = Am
    const nB = length(Bcolptr)-1
    const rowvalA = Arowval; colptrA = Acolptr
    const rowvalB = Browval; colptrB = Bcolptr
    const preallocationIncrement = colptrA[end]+colptrB[end]

	colptrC = Array{Tv}( nB+1)
    colptrC[1] = 1
	rowSupp = zeros(Tv, mA)
	rowList = Array{Tv}( mA)
	rowvalCj = Array{Bool}( mA)
	rowvalC = Array{Tv}( preallocationIncrement)
	for i in 1:nB
		newrowscounter = 0
		eyerange = cran(Bcolptr,i)
		newrowscounter = length(eyerange)
		for l=1:newrowscounter
			ll = rowvalB[eyerange[l]]
			rowList[l] = ll
			rowSupp[ll] = i
			rowvalCj[ll]=true
		end
		for jp in colptrB[i]:(colptrB[i+1] - 1)
			j = rowvalB[jp]
			for kp in colptrA[j]:(colptrA[j+1] - 1)
				k = rowvalA[kp]
				if rowSupp[k] != i
					rowSupp[k] = i
					newrowscounter +=1
					rowList[newrowscounter] = k
					rowvalCj[k] = true
				else
					rowvalCj[k] = !rowvalCj[k]
				end
			end
		end
		nzRows = find(rowvalCj[rowList[1:newrowscounter]])
		colptrC[i+1] = colptrC[i]+length(nzRows)

		if colptrC[i+1] > length(rowvalC)+1
			append!(rowvalC, Array{Int}(preallocationIncrement))
		end
		rowvalC[colptrC[i]:(colptrC[i+1]-1)] = sort(rowList[nzRows])
	end
	deleteat!(rowvalC,colptrC[end]:length(rowvalC))
	return rowvalC, colptrC
end

function blockprodsumWrite2!(Crv,Ccp,Brv,Bcp,Drv,Dcp,Dm,Dn,Mrv,Mcp,preallocationIncrement)
	extend!(Mcp,Dn+1)
	deleteat!(Mcp,(Dn+2):length(Mcp))
	Mcp[1]=1

	rowSupp = zeros(Int64, Dm)
	rowList = Array{Int64}( Dm)
	rowvalCj = BitArray(Dm)
	for i in 1:Dn
		if length(Mrv)<Mcp[i]+Dm
			extend!(Mrv,length(Mrv)+Dm+preallocationIncrement)
		end
		eyerange = cran(Dcp,i)
		newrowscounter = 0
		for kp = eyerange
			row = Drv[kp]
			rowSupp[row]=i
			newrowscounter+=1
			rowList[newrowscounter]=row
			rowvalCj[row] = true
		end
		if Bcp[i]<Bcp[Dn+1]
			for jp in Bcp[i]:(Bcp[i+1] - 1)
				j = Brv[jp]
				for kp in Ccp[j]:(Ccp[j+1] - 1)
					row = Crv[kp]
					if rowSupp[row] != i
						rowSupp[row] = i
						newrowscounter +=1
						rowList[newrowscounter] = row
						rowvalCj[row] = true
					else
						rowvalCj[row] = !rowvalCj[row]
					end
				end
			end
		end
		ii = i+1
		Mcp[ii]=Mcp[i]
		for k = 1:newrowscounter
			row = rowList[k]
			if rowvalCj[row]
				Mrv[Mcp[ii]]=row
				Mcp[ii]+=1
			end
		end
	end
end

function blockprodsumsilenticolsleftWrite2!(Crv,Ccp,Brv,Bcp,Drv,Dcp,Dm,Dn,Mrv,Mcp,preallocationIncrement,col2silenti)
	extend!(Mcp,Dn+1)
	deleteat!(Mcp,(Dn+2):length(Mcp))
	Mcp[1]=1

	rowSupp = zeros(Int64, Dm)
	rowList = Array{Int64}( Dm)
	rowvalCj = BitArray(Dm)
	for i in 1:Dn
		if length(Mrv)<Mcp[i]+Dm
			extend!(Mrv,length(Mrv)+Dm+preallocationIncrement)
		end
		eyerange = cran(Dcp,i)
		newrowscounter = 0
		for kp = eyerange
			row = Drv[kp]
			rowSupp[row]=i
			newrowscounter+=1
			rowList[newrowscounter]=row
			rowvalCj[row] = true
		end
		for jp in Bcp[i]:(Bcp[i+1] - 1)
			j = Brv[jp]
			row = col2silenti[j]
			if rowSupp[row] != i
				rowSupp[row] = i
				newrowscounter +=1
				rowList[newrowscounter] = row
				rowvalCj[row] = true
			else
				rowvalCj[row] = !rowvalCj[row]
			end
		end
		for jp in Bcp[i]:(Bcp[i+1] - 1)
			j = Brv[jp]
			for kp in Ccp[j]:(Ccp[j+1] - 1)
				row = Crv[kp]
				if rowSupp[row] != i
					rowSupp[row] = i
					newrowscounter +=1
					rowList[newrowscounter] = row
					rowvalCj[row] = true
				else
					rowvalCj[row] = !rowvalCj[row]
				end
			end
		end
		ii = i+1
		Mcp[ii]=Mcp[i]
		for k = 1:newrowscounter
			row = rowList[k]
			if rowvalCj[row]
				Mrv[Mcp[ii]]=row
				Mcp[ii]+=1
			end
		end
	end
end

##########################################################################################

####	COPY, SHIFT, INDEX, AND SLICE

##########################################################################################

function cran(A::SparseMatrixCSC,j)
	return A.colptr[j]:(A.colptr[j+1]-1)
end

function cran(colptr::Array,j::Int64)
	return colptr[j]:(colptr[j+1]-1)
end

#	added 12/27/2017
# 	may be tested via function <testcran>
function cran(colptr::Array,J::Array{Int64,1})
	m = nval(colptr,J)
	v = zeros(Int64,m)
	c = 0
	for p=1:length(J)
		k = nval(colptr,J[p])
		v[c+1:c+k]=cran(colptr,J[p])
		c += k
	end
	return v
end

#	added 12/27/2017
# 	may be tested via function <testcran>
function cran(colptr::Array,J::UnitRange{Int64})
	m = nval(colptr,J)
	v = zeros(Int64,m)
	c = 0
	for p=1:length(J)
		k = nval(colptr,J[p])
		v[c+1:c+k]=cran(colptr,J[p])
		c += k
	end
	return v
end

#	added 12/27/2017
function testcran(m,n)
	for p = 1:m
		a = rand(n,n).<0.1
		a = convert(Array{Int64},a)
		J1= rand(1:n,5)
		J2= minimum(J1):maximum(J1)
		v1= []
		v2= []
		for q in J1
			k = countnz(a[:,1:(q-1)])
			d = countnz(a[:,q])
			append!(v1,k+1:k+d)
		end
		for q in J2
			k = countnz(a[:,1:(q-1)])
			d = countnz(a[:,q])
			append!(v2,k+1:k+d)
		end
		rv,cp = full2ss(a)
		if v1 != cran(cp,J1)
			print("error at J1"); return a,v1,J1
		elseif v2 != cran(cp,J2)
			print("error at J2"); return a,v1,J2
		end
	end
	print("test successful")
end

function cran(colptr::UnitRange,j)
	return colptr[j]
end

#	added 12/27/2017
"""
	nval(colptr,j::Int64)

For a column sparse matrix with column pattern `colptr`, counts the number of values stored for column `j`.
"""
function nval(colptr,j::Int64)
	return colptr[j+1]-colptr[j]
end

#	added 12/27/2017
"""
	nval(colptr,J)

For a column sparse matrix with column pattern `colptr`, counts the number of values stored for column j, for each j in `J`.
"""
function nval(colptr,J)
	c = 0
	for p = 1:length(J)
		c += nval(colptr,J[p])
	end
	return c
end

#	added 12/27/2017
function nvaltest(m,n)
	for p = 1:m
		a = rand(n,n).<0.1
		a = convert(Array{Int64},a)
		rv,cp = full2ss(a)
		J1 = rand(1:n,5)
		J2 = minimum(J1):maximum(J1)
		if nval(cp,J1)!=countnz(a[:,J1])
			print("error"); return a,cp,J1
		elseif nval(cp,J2)!=countnz(a[:,J2])
			print("error"); return a,cp,J2
		end
	end
	print("test successful")
end

function crows(A::SparseMatrixCSC,j)
	return A.rowval[cran(A,j)]
end

function crows(colptr::Array,rowval::Array,j)
	return rowval[cran(colptr,j)]
end

function findcol(cp,k)
	i = 1
	while cp[i]<=k
		i+=1
	end
	return i-1
end

function extend!{Tv}(x::Array{Tv,1},n::Integer)
	if length(x)<n
		append!(x,Array{Tv}(n-length(x)))
	end
end

function copycolumnsubmatrix{Tv<:Integer}(Arv::Array{Tv,1},Acp,columnindices)
	allocationspace = 0
	for j in columnindices
		allocationspace+= Acp[j+1]-Acp[j]
	end
	Brv = Array{Tv}(allocationspace)
	Bcp = Array{Int64}(length(columnindices)+1)
	Bcp[1]=1
	for jp = 1:length(columnindices)
		j = columnindices[jp]
		Bcp[jp+1] = Bcp[jp]+Acp[j+1]-Acp[j]
		Brv[Bcp[jp]:(Bcp[jp+1]-1)]=Arv[Acp[j]:(Acp[j+1]-1)]
	end
	return Brv,Bcp
end

# NB: the following was a facsimile of the preceding funciton, written at a time
# when it was regarded as a good style in Julia to be extremely specific with
# regard to input type.
# function copycolumnsubmatrix{Tv<:Integer}(Arv::Array{Tv,1},Acp::Array{Tv,1},columnindices::UnitRange{Int64})
# 	allocationspace = 0
# 	for j in columnindices
# 		allocationspace+= Acp[j+1]-Acp[j]
# 	end
# 	Brv = Array{Tv}(allocationspace)
# 	Bcp = Array{Tv}(length(columnindices)+1)
# 	Bcp[1]=1
# 	for jp = 1:length(columnindices)
# 		j = columnindices[jp]
# 		Bcp[jp+1] = Bcp[jp]+Acp[j+1]-Acp[j]
# 		Brv[Bcp[jp]:(Bcp[jp+1]-1)]=Arv[Acp[j]:(Acp[j+1]-1)]
# 	end
# 	return Brv,Bcp
# end

function copycolind2colind!{Tv<:Integer}(rowvalA::Array{Tv,1},colptrA::Array{Tv,1},columnindices,rowvalB::Array{Tv,1},colptrB::Array{Tv,1},startingDestination::Integer,growthIncrement::Integer)
	numnewrows = 0
	for col in columnindices
		numnewrows+=colptrA[col+1]-colptrA[col]
	end
	if length(rowvalB)<colptrB[startingDestination]+numnewrows-1
		append!(rowvalB,Array{Tv}(numnewrows+growthIncrement))
	end
	if length(colptrB)<startingDestination+length(columnindices)
		append!(colptrB,Array{Tv}(startingDestination+length(columnindices)))
	end
	colptrB[1]=1
	for i = 1:length(columnindices)
		k = startingDestination+i #the index of colptr pointing to end of this new column
		col = columnindices[i]
		colptrB[k]=colptrB[k-1]+colptrA[col+1]-colptrA[col]
		rowvalB[colptrB[k-1]:(colptrB[k]-1)]=rowvalA[colptrA[col]:(colptrA[col+1]-1)]
	end
end

function extendcolumnlight!{Ti}(rowval::Array{Ti,1},colptr::Array{Ti,1},v::Array{Ti},k::Ti,growthincrement::Ti)
	r = rowval
	c = colptr
	startpoint = copy(c[k+1])
	c[k+1]=c[k]+length(v)
	if length(r)<c[k+1]-1
		append!(r,Array{Int64}(max(growthincrement,length(v))))
	end
	r[startpoint:(c[k+1]-1)]=v
end

# colsinorder must be in sorted order
function supportedmatrix!{Tv<:Integer}(Mrowval::Array{Tv},Mcolptr::Array{Tv,1},rows1,colsinorder,Mm::Tv)
	n = length(colsinorder)
	suppcol1 = falses(Mm)
	suppcol1[rows1]=true
	cpHolder = 1
	nz1 = 0
	for jp = 1:n
		for ip in cran(Mcolptr,colsinorder[jp])
			i = Mrowval[ip]
			if suppcol1[i]
				nz1+=1
				Mrowval[nz1]=i
			end
		end
		Mcolptr[jp]=cpHolder
		cpHolder = nz1+1
	end
	Mcolptr[n+1] = cpHolder
	deleteat!(Mcolptr,(n+2):length(Mcolptr))
	deleteat!(Mrowval,Mcolptr[end]:length(Mrowval))
end

function stackedsubmatrices{Tv<:Integer}(
	Mrowval,#::Array{Tv,1},
	Mcolptr,#::Array{Tv,1},
	rows1,#::Array{Tv,1},
	rows2,#::Array{Tv,1},
	cols,#::Array{Tv,1},
	Mm::Tv)
	n = length(cols)
	suppcol1 = falses(Mm)
	suppcol2 = falses(Mm)
	suppcol1[rows1]=true
	suppcol2[rows2]=true
	nz1 = 0; nz2 = 0
	for jp = 1:n
		for ip in cran(Mcolptr,cols[jp])
			i = Mrowval[ip]
			if suppcol1[i]
				nz1+=1
			elseif suppcol2[i]
				nz2+=1
			end
		end
	end
	rv1 = Array{Tv}(nz1)
	rv2 = Array{Tv}(nz2)
	cp1 = Array{Tv}(n+1); cp1[1]=1
	cp2 = Array{Tv}(n+1); cp2[1]=1
	nz1 = 0; nz2 = 0
	for jp = 1:n
		for ip in cran(Mcolptr,cols[jp])
			i = Mrowval[ip]
			if suppcol1[i]
				nz1+=1
				rv1[nz1]=i
			elseif suppcol2[i]
				nz2+=1
				rv2[nz2]=i
			end
		end
		cp1[jp+1] = nz1+1
		cp2[jp+1] = nz2+1
	end
	return rv1,rv2,cp1,cp2
end

##########################################################################################

####	TRANSPOSITION OPERATORS

##########################################################################################

function transposeLighter{Ti}(Arowval::Array{Ti},Acolptr::Array{Ti},Am)
    Annz = Acolptr[end]-1
    An = length(Acolptr)-1
    # Attach destination matrix
    Cm = An
    Cn = Am
    Ccolptr = Array{Ti}(Am+1)
    Crowval = Array{Ti}(Annz)
    # Compute the column counts of C and store them shifted forward by one in
	# Ccolptr
    Ccolptr[1:end] = 0
    @inbounds for k in 1:Annz
        Ccolptr[Arowval[k]+1] += 1
    end
    # From these column counts, compute C's column pointers
    # and store them shifted forward by one in Ccolptr
    countsum = 1
    @inbounds for k in 2:(Cn+1)
        overwritten = Ccolptr[k]
        Ccolptr[k] = countsum
        countsum += overwritten
    end
    # Distribution-sort the row indices and nonzero values into Crowval and
	# Cnzval, tracking write positions in Ccolptr
    for Aj in 1:An
        for Ak in Acolptr[Aj]:(Acolptr[Aj+1]-1)
            Ai = Arowval[Ak]
            Ck = Ccolptr[Ai+1]
            Crowval[Ck] = Aj
            Ccolptr[Ai+1] += 1
        end
    end
    # Tracking write positions in Ccolptr as in the last block fixes the colptr
	# shift, but the first colptr remains incorrect
    Ccolptr[1] = 1

	return Crowval, Ccolptr
end

function transposeLighter{Tv,Ti}(Arowval::Array{Ti},Acolptr::Array{Ti},Anzval::Array{Tv},Am::Integer)
    Annz = Acolptr[end]-1
    An = length(Acolptr)-1
    Cm = An
    Cn = Am
    Ccolptr = Array{Ti}(Am+1)
    Crowval = Array{Ti}(Annz)
    Cnzval = Array{Tv}(Annz)
    # Compute the column counts of C and store them shifted forward by one in Ccolptr
    Ccolptr[1:end] = 0
    @inbounds for k in 1:Annz
        Ccolptr[Arowval[k]+1] += 1
    end
    # From these column counts, compute C's column pointers
    # and store them shifted forward by one in Ccolptr
    countsum = 1
    @inbounds for k in 2:(Cn+1)
        overwritten = Ccolptr[k]
        Ccolptr[k] = countsum
        countsum += overwritten
    end
    # Distribution-sort the row indices and nonzero values into Crowval and Cnzval,
    # tracking write positions in Ccolptr
    @inbounds for Aj in 1:An
        for Ak in Acolptr[Aj]:(Acolptr[Aj+1]-1)
            Ai = Arowval[Ak]
            Ck = Ccolptr[Ai+1]
            Crowval[Ck] = Aj
            Cnzval[Ck] = Anzval[Ak]
            Ccolptr[Ai+1] += 1
        end
    end
    # Tracking write positions in Ccolptr as in the last block fixes the colptr shift,
    # but the first colptr remains incorrect
    Ccolptr[1] = 1
	return Crowval, Ccolptr, Cnzval
end

#=
-	Returns the sparse equivalent of A[rows,cols]'.
- 	The rows of the output array will be listed in the order that they appear in
	input vector <cols>.
- 	Does not assume that the entries of <rows> appear in sorted order.
- 	I *believe* that repeated rows and columns are allowed.
=#
function transposeLighter_submatrix{Ti}(Arowval::Array{Ti},Acolptr::Array{Ti},Am;rows = 1:Am,cols = 1:length(Acolptr)-1)
	if rows == 1:Am && cols == 1:(length(Acolptr)-1)
		Crowval, Ccolptr = transposeLighter(Arowval::Array{Ti},Acolptr::Array{Ti},Am)
		return Crowval, Ccolptr
	end
    # Attach destination matrix
    Cm = length(cols)
    Cn = length(rows)
    Ccolptr = Array{Ti}(Cn+1)
    # Compute the column counts of C and store them shifted forward by one in Ccolptr
    Ccolptr[1:end] = 0
	rs = rowsupportsum(Arowval,Acolptr,Am,cols)
	for i = 1:Cn
	    Ccolptr[i+1] = rs[rows[i]]
	end
	Cnnz = sum(Ccolptr)
    Crowval = Array{Ti}(Cnnz)
    # From these column counts, compute C's column pointers
    # and store them shifted forward by one in Ccolptr
    countsum = 1
    @inbounds for k in 2:(Cn+1)
        overwritten = Ccolptr[k]
        Ccolptr[k] = countsum
        countsum += overwritten
    end
    # Distribution-sort the row indices and nonzero values into Crowval and Cnzval,
    # tracking write positions in Ccolptr
    rowtranslator = zeros(Int64,Am)
    rowtranslator[rows] = 1:Cn
    @inbounds for Cj in cols
        for Ak in cran(Acolptr,Cj)
            Ai = rowtranslator[Arowval[Ak]]
            if Ai > 0
	            Ck = Ccolptr[Ai+1]
	            Crowval[Ck] = Cj
    	        Ccolptr[Ai+1] += 1
    	    end
        end
    end
    # Tracking write positions in Ccolptr as in the last block fixes the colptr shift,
    # but the first colptr remains incorrect
    Ccolptr[1] = 1
	return Crowval, Ccolptr
end

#=
Accepts an mxn integer array, M.  To M we implicitly associate an array N,
as follows. If M has no nonzero entries, then N is the zeron-one array such that
supp(N[:,i]) = M[:,i].  Here we regard M[:,i] as a set.  I *believe* one ignores
duplicate entries, but have not checked. A similar interpretation holds when M
has zero entries - one simply discards the zeros.
Let S denote the support of N, r denote row02row1translator, and c denote
col02col1translator.  This function returns the data specifying a sparse
zero-one matrix whose support is
{ [c[j],r[i]] : [i,j] \in N and neither c[j] nor r[i] are zero }.
=#
function presparsefull2unsortedsparsetranspose{Tv<:Integer}(
	M::Array{Tv,2},
	row02row1translator,
	col02col1translator;
	verbose::Bool=false)
	Mm,Mn = size(M)

	if Mn == 0
		rowval1 = Array{Int64}(0)
		if isempty(row02row1translator)
			colptr1 = ones(Int64,1)
		else
			colptr1 = ones(Int64,1+maximum(row02row1translator))
		end
		return rowval1,colptr1,Mn
	end

	for i = 1:(Mm*Mn)
		if M[i]>0
			M[i] = row02row1translator[M[i]]
		end
	end
	m0 = maximum(row02row1translator)
	rowcounter = zeros(Int64,m0)
	for k in M  #allow for zero values
		if k>0
			rowcounter[k]+=1
		end
	end
	colptr1 = Array{Int64}(m0+1)
	colptr1[1]=1
	for i = 1:m0
		colptr1[i+1]=colptr1[i]+rowcounter[i]
	end
	rowval1 = Array{Int64}(colptr1[end]-1)
	placer = copy(colptr1)
	if verbose
		coverageTestVec = trues(colptr1[end]-1)
	end
	for j = 1:Mn
		jval = col02col1translator[j]
		for i = 1:Mm
			row = M[i,j]
			if row > 0
				if verbose
					coverageTestVec[placer[row]]=false
				end
				rowval1[placer[row]]=jval
				placer[row]+=1
			end
		end
	end
	if verbose
		if any(coverageTestVec)
			print("please refer to the coverageTestVec")
			sleep(4)
		end
		println([length(rowval1) "length(output of presparsefull2unsortedsparsetranspose)"])
	end
	#gc()
	return rowval1,colptr1,Mn
end

##########################################################################################

####	SHAPE GENERATORS

##########################################################################################

function noisycircle()
	theta = 1:100
	theta = theta*2*pi/100
	x = cos.(theta)
	y = sin.(theta)
	pcloud = hcat(x,y)'
	pcloud = hcat(pcloud,pcloud,pcloud)
	pcloud = pcloud + 0.3*rand(2,300)
	return pcloud
end

function noisycircle3()
	theta = 1:100
	theta = theta*2*pi/100
	x = cos.(theta)
	y = sin.(theta)
	z = zeros(100)
	pcloud = hcat(x,y,z)'
	pcloud = hcat(pcloud,pcloud,pcloud)
	pcloud = pcloud + 0.3*rand(3,300)
	return pcloud
end

function torus(;m = 100,n = 50,mrad=1,nrad = 0.2)
	theta = (1:m)*2*pi/m;
	torus1 = mrad*repmat(cos.(theta),1,n)
	torus2 = mrad*repmat(sin.(theta),1,n)
	torus3 = nrad*repmat(sin.(2*pi*(1:n)./n),1,m)'
	torus4 = repmat(cos.(2*pi*(1:n)./n),1,m)'
	torus1 = torus1+nrad*torus1.*torus4;
	torus2 = torus2+nrad*torus2.*torus4;
	return hcat(torus1[:],torus2[:],torus3[:])'
end

function noisytorus(;m = 100,n=50,mrad=1,nrad = 0.2,noiserad= 0.5*nrad)
	pcloud = torus(m=m, n=n,nrad = nrad,mrad=mrad)
	pcloud = pcloud + noiserad*rand(size(pcloud,1),size(pcloud,2))
end

function sphere()
	halfnumrings = 40
	latitd = pi*(-halfnumrings:halfnumrings)/(2*halfnumrings)

	x = Array{Float64}(0)
	y = Array{Float64}(0)
	z = Array{Float64}(0)

	for t in latitd
		numpts = round.(Int64,40*cos.(t))
		theta = 2*pi*(1:numpts)/numpts
		append!(x,cos.(t)*cos.(theta))
		append!(y,cos.(t)*sin.(theta))
		append!(z,fill(sin.(t),numpts))
	end
	return x,y,z
end

function matchingcomplex_symmat(m,n)
	#### 1-A-I, where A is the adjacency matrix of the matching complex M(m,n)
	#### and I is identity
	combs = collect(combinations(Array(1:n),m))  #an array of arrays
	numcombs = binom(n,m)
	numpointedcombs = binom(n-1,m-1)
	symmat = zeros(Int64,numcombs,numcombs)
	indexvec = zeros(Int64,numpointedcombs)
	marker = 0

	for i = 1:n
		marker = 0
		for j = 1:numcombs
			if i in combs[j]
				marker+=1
				indexvec[marker]=j
			end
		end
		if marker!=length(indexvec) || marker!=numpointedcombs
			print("please check construction algorithm")
			sleep(10)
		end
		symmat[indexvec,indexvec] = 1
	end
	for i = 1:numcombs
		symmat[i,i]=0
	end
	return symmat
end

function chessboardcomplex_symmat(;numrows=3,numcols=4)
	#### 1-A-I, where A is the adjacency matrix of the chessboard complex C(m,n)
	#### and I is identity
	m 		= numrows
	n 		= numcols
	numrooks = m*n
	symmat = ones(Int64,numrooks,numrooks)
	for i = 1:m
		for j = 1:n
			rooknum1 = (i-1)*m+j
			for ii = 1:m
				for jj = 1:n
					if ii != i && jj !=j
						rooknum2 = (ii-1)*m+jj
						symmat[rooknum1,rooknum2] = 0
						symmat[rooknum2,rooknum1] = 0
					end
				end
			end
		end
	end
	for k = 1:numrooks
		symmat[k,k] = 0
	end
	return symmat
end

function plane2torus(A)
	theta1 = A[:,1]
	theta2 = A[:,2]
	x = cos.(theta1)
	y = sin.(theta1)
	z = 0.25*sin.(theta2)
	alpha = 1-0.25*cos.(theta2)
	x = x[:]'
	y = y[:]'
	z = z[:]'
	alpha = alpha[:]'
	x = alpha.*x
	y = alpha.*y
	return vcat(x,y,z)
end

function latlon2euc(A;model = "pc")
	if model == "pc"
		theta1 = A[1,:]
		theta2 = A[2,:]
	elseif model == "points"
		theta1 = A[:,1]'
		theta2 = A[:,2]'
	end
	x = cosd.(theta2)
	y = sind.(theta2)
	z = sind.(theta1)
	w = cosd.(theta1)
	x = w.*x
	y = w.*y
	return vcat(x,y,z)
end

function worldexample()
	a = readdlm("/Users/greghenselman/Desktop/citylonglat.csv",',',Float64,'\r')
	d = latlon2sphere(a')
	x = JLD.load("/Users/greghenselman/JuliaFiles/citysample3_copy.jld")
	supp = x["supp"]
	C = eirene(d[:,supp];model = "pc",maxdim=1,upperlim = 0.25,record="cyclerep")
	return a,d,supp,C
end

function zerodrandmat(n)
	# input: 	an integer n
	# output: 	symmetric matrix with zeros on the diagonal and iid uniform entries elsewhere
	x = rand(n,n)
	for p = 1:n
		x[p,p]=0
	end
	x = (x+x')/2
	return x
end

function eirenefilepath(filedescription)
	if 		filedescription 	== 	"simplecity"
			return joinpath(@__DIR__,"examples/simplemapscitydata.csv")
	elseif 	filedescription 	== 	"noisycircle"
			return joinpath(@__DIR__,"examples/noisycircle.csv")
	elseif 	filedescription 	== 	"noisytorus"
			return joinpath(@__DIR__,"examples/noisytorus.csv")
	end
end

##########################################################################################

####	MATRIX WEIGHTS AND FORMATTING

##########################################################################################

#=
The output of this function is obtatined by *first* rounding all values below
minrad up to minrad, *then* rounding all values above maxrad up to Inf, *then*
rounding the entries valued in [minrad,maxrad] to the nearest greater element of
[minrad:numrad:maxrad], where by definition [minrad:Inf:maxrad] is the closed
inteveral containing all reals between minrad and maxrad.
If minrad == -Inf, then we set it to minimum(N) before performing any operations.
If maxrad == Inf, then we set it to maximum(N) before performing any operations.
=#
function minmaxceil(N;minrad=minimum(N),maxrad=maximum(N),numrad=Inf)
	S 	= 	copy(N) # NB: this is probably an important step; several pernicious problems in the development of rounding procedures turned out to be rooted in accidental rewriting of matrix entries
	return minmaxceil!(S;minrad=minrad,maxrad=maxrad,numrad=numrad)
end

#=
NB THIS FUNCTION MODIFIES ITS INPUT

The result of this function is obtatined by *first* rounding all values below
minrad up to minrad, *then* rounding all values above maxrad up to Inf, *then*
rounding the entries valued in [minrad,maxrad] to the nearest greater element of
[minrad:numrad:maxrad], where by definition [minrad:Inf:maxrad] is the closed
inteveral containing all reals between minrad and maxrad.
If minrad == -Inf, then we set it to minimum(N) before performing any operations.
If maxrad == Inf, then we set it to maximum(N) before performing any operations.
=#
function minmaxceil!(S;minrad=minimum(N),maxrad=maximum(N),numrad=Inf)

	if 	(minrad == Inf) || (maxrad == -Inf)
		return fill(Inf,size(S)...)
	end

	if minrad == "minedge"
		minrad = minimum(offdiagmin(S))
	end

	S[S.<minrad] 	= 	minrad
	S[S.>maxrad] 	= 	Inf

	fi 	= 	find(isfinite,S) # stands for finite indices
	fv 	= 	S[fi]

	if isempty(fv)
		return S
	end

	if minrad 	== 	-Inf
		minrad 	=  	minimum(fv)
	end
	if maxrad 	== 	Inf
		maxrad 	= 	maximum(fv)
	end
	if numrad 	== 	1
		S[fi]	= 	maxrad
		return 	S
	end
	if numrad 	== 	Inf
		return S
	end

	ran 		= 	linspace(minrad,maxrad,numrad)
	ran 		= 	Array{Float64}(ran)
	ran[end] 	= 	maxrad

	fvr 		= 	ceil2grid(fv,ran)

	S[fi]		= 	fvr

	return 		S
end

# ran should be an array in sorted order, with ran[end] >= maximum(A)
function ceil2grid(A,ran)
	if 	ran 			== 		"all"
		return A
	end
	B 					= 		copy(A)
	for j 				= 		1:length(A)
		post 			= 		1
		while ran[post] < A[j]
			post+=1
		end
		B[j] 			= 		ran[post]
	end
	return B
end

function checkceil2grid(numits)
	for 	p 	=	1:numits
		A 		= 	rand(70)
		ran 	= 	sort(rand(20))
		ran 	= 	ran*(1/maximum(ran)) # this guarantees that maximum(ran) > maximum(A)
		crct 	= 	crosscheckceil2grid(A,ran)
		if !crct
			println()
			println("error: please check checkceil2grid")
		end
	end
	return []
end

function crosscheckceil2grid(A,ran)
	B 		= 	ceil2grid(A,ran)
	for j 	= 	1:length(A)
		q 	= 	findfirst(ran.>A[j])
		val = 	ran[q]
		if B[j] 	!= 	val
			println()
			println("error: please check <checkceil2grid>")
			return false
		end
	end
	return true
end

function ocfcheckfun3()
	n = 10
	m = 1000
	for q = 1:n
		M = rand(m,m)
		M = M+M'
		if q < n-5
			numrad = rand(10:90)
			upperlim = rand()*2
		elseif q == n-5
			numrad = 10
			upperlim = rand()*2
		elseif q == n-4
			numrad = 1
			upperlim = -1
		elseif q == n-3
			numrad = 1
			upperlim = rand()*2
		elseif q == n-2
			numrad = 1
			upperlim = 0
		elseif q == n-1
			numrad = 1
			upperlim = 3
		elseif q == n
			numrad = 1
			upperlim = -1
		end

		ocf1 = Array{Any,1}(8)
		ocf2 = Array{Any,1}(8)
		j = 0
		for a in [Inf, upperlim]
			for b in [Inf,numrad]
				for c in [true, false]
					j+=1
					print("RRR1[[j=[$(j)]]]")
					ocf1[j] = ordercanonicalform_3(M,maxrad=a,numrad=b,fastop=c)
				end
			end
		end
		print("ocf3-done-")

		if numrad > 1
			N = copy(M)
			for i = 1:m
				N[i,i] = Inf
			end
			minentry = minimum(N)
			for i = 1:m
				N[i,i] = -Inf
			end
			maxentry = maximum(N)

			Q = minmaxceil(copy(M),minentry,upperlim,numrad)
			N = minmaxceil(copy(M),minentry,maxentry,numrad)
			for i = 1:m
				N[i,i] = Inf
				Q[i,i] = Inf
			end
		end
		if numrad == 1
			N = copy(M)
			if upperlim == Inf
				for i = 1:m
					N[i,i] = -Inf
				end
				maxrad = maximum(N)
			else
				maxrad = upperlim
			end
			N[N.<=maxrad]=maxrad
			N[N.>maxrad]= maxrad+1
			for i = 1:m
				N[i,i] = Inf
			end
			Q = copy(N)
		end

		j = 0
		for a in [Inf, upperlim]
			for b in [Inf,numrad]
				if b == numrad
					if a == upperlim
						L = copy(Q)
					else
						L = copy(N)
					end
				else
					L = copy(M)
				end

				for c in [true, false]
					j+=1
					LL = copy(M)

					for i = 1:m
						LL[i,i] = Inf
					end
					publicmin = minimum(LL)

					for i = 1:m
						LL[i,i]=-Inf
					end
					if a == Inf
						publicmax = maximum(LL)
					else
						publicmax = a
					end

					if publicmax < publicmin
						ocf2[j] = zeros(Int64,m,m),Array{Float64}(0)
					elseif b == 1 && c
						privatemax = minimum(maximum(LL,1))
						if publicmax >= privatemax
							ocf = zeros(Int64,m,m)
							index = findfirst(maximum(LL,1),privatemax)
							ocf[index,:]=1
							ocf[:,index]=1
							ocf[index,index]=0
							ocf2[j]= (ocf,[publicmax],M)
						else
							ocf = zeros(Int64,m,m)
							ocf[LL.<=publicmax]=1
							for i=1:m
								ocf[i,i]=0
							end
							ocf2[j] = (ocf,[publicmax],M)
						end
					elseif b == 1 && a == Inf
						ocf = ones(Int64,m,m)
						for i=1:m
							ocf[i,i]=0
						end
						ocf2[j] = (ocf,[maximum(LL)],M)
					else
						ocf2[j] = ordercanonicalform(L,maxrad=a,fastop=c)
					end
				end
			end
		end
		for k = 1:8
			if ocf1[k][1]!=ocf2[k][1] || sum(abs(ocf1[k][2]-ocf2[k][2])) > 0.00000001
				println("locus = [$(k)]")
				println("q = [$(q)]")
				println("numrad = $(numrad)")
				println("upperlim = $(upperlim)")
				return ocf1,ocf2,M,N
			end
		end
	end
	return "passedtest"
end

#=
minimum(deleteat(S[:,i],i))
=#
function offdiagmin(S,i)
	if i == 1
		return(minimum(S[2:end,i]))
	elseif i == size(S,1)
		return minimum(S[1:end-1,i])
	else
		return min(minimum(S[1:i-1,i]),minimum(S[i+1:end,i]))
	end
end

function checkoffdiagmin(numits)
	for p 					= 	1:numits
		S 					= 	rand(50,50)
		for q 				= 	1:50
			checkval 		= 	minimum(deleteat!(S[:,q],q))
			if checkval 	!= 	offdiagmin(S,q)
				return S
			end
		end
	end
	return zeros(Int64,0)
end

function offdiagmean(S;defaultvalue=[])
	m,n 	= 	size(S)
	if m 	!= 	n
		println()
		println("error in <offdiagmean>: input matrix should be square")
	end
	if isempty(defaultvalue)
		println()
		println("warning: no defaulvalue was set for <offdiagmin>; this parameter has been set to 0")
		defaultvalue 	= 	0
	end
	if m 		== 	1
		return 	defaultvalue
	end
	mu 			= 	zeros(m)
	for j 		= 	1:m
		v 		= 	S[1:(j-1),j]
		u 		= 	S[(j+1):m,j]
		mu[j]  	= 	mean(vcat(v[:],u[:]))
	end
	return 		mu
end

# NB: Assumes that the input array S has only finite entries.
# NB: The value for keyword <numrad> must be either a positive integer or Inf
function ordercanonicalform_3{Tv}(
	S::Array{Tv};
	maxrad=Inf,
	minrad=-Inf,
	numrad=Inf,
	vscale="default",
	fastop::Bool=true,
	verbose::Bool=false)

    if size(S,1) == 0
    	ocf 		= 	zeros(Int64,0,0)
    	ocg2rad		= 	zeros(Float64,0)
    	return ocf,ocg2rad
    end

    # Format input
    symmat		= convert(Array{Float64},copy(S))
	m 			= size(symmat,1)

	if vscale == "default"
		for i = 1:m
			symmat[i,i] = minimum(symmat[:,i])
		end
	elseif typeof(vscale) <: Array
		vscale = convert(Array{Float64},copy(vscale))
		if length(vscale) != m
			print("Error: keyword <vscale> must take a value of <defualt>, <diagonal>, or <v>, where v is a vector of length equal to the number of vertices in the complex.")
			return
		end
		for i=1:m
			if offdiagmin(symmat,i) < vscale[i]
				print("Error: the \"birth time\" assigned a vertex by keyword argument <vscale> may be no greater than that of any incident edge.  The value assigned to vertex $(i) by <vscale> is $(vscale[i]), and i is incident to an edge with birth time $(offdiagmin(symmat,i)).")
				return
			else
				symmat[i,i] = vscale[i]
			end
		end
	elseif vscale 	== 	"diagonal"
		vscale		=	Array{Float64,1}(m)
		for i=1:m
			vscale[i]	=	symmat[i,i]
		end
		# the following is in prnciple unnecessary, but it simplifies rounding
		for i=1:m
			if offdiagmin(symmat,i) < vscale[i]
				print("Error: the \"birth time\" assigned a vertex by keyword argument <vscale> may be no greater than that of any incident edge.  The value assigned to vertex $(i) by <vscale> is $(vscale[i]), and i is incident to an edge with birth time $(offdiagmin(symmat,i)).")
				return
			end
		end
	end

	# Deterime the public maxrad
	if maxrad == Inf
		publicmax = maximum(symmat)
	else
		publicmax 	= copy(maxrad)
	end

	# Deterime the public minrad
	publicmin	= minimum(symmat)
	publicmin	= max(publicmin,minrad)

	# It's important that this precede the other cases
	if publicmax < publicmin
		return zeros(Int64,m,m),Array{Float64}(0)
	end

	# This covers all remaining cases where numrad ==1
	if numrad == 1
		privatemax = minimum(maximum(symmat,1))
		ocf = zeros(Int64,m,m)
		if fastop && publicmax >= privatemax
			index = findfirst(maximum(symmat,1),privatemax)
			ocf[index,:]=1
			ocf[:,index]=1
			for i = 1:m
				ocf[i,i]=1
			end
		else
			ocf[symmat.<=publicmax]=1
		end
		return ocf,[publicmax]
	end

	# If necessary, determine step size.  Recall we have already treated every case where numrad == 1.
	if numrad < Inf
		alpha 		= (publicmax-publicmin)/(numrad-1)
	elseif numrad == Inf
		alpha 		= Inf
	end

	# If one stops early, determine when to stop
	if fastop
		privatemax = minimum(maximum(symmat,1))
		privatemax = min(privatemax,publicmax)
		if numrad < Inf
			post = publicmin
			stepcounter = 1
			while post < privatemax
				stepcounter+=1
				if stepcounter == numrad
					post = publicmax # must take this rather cumbersome step on account of numerical error.
				else
					post+=alpha
				end
			end
			privatemax = post
		end
	else
		privatemax = publicmax
	end

	# Extract sortperm
	p 						= sortperm(vec(symmat),alg=MergeSort)

	# Compute the ocf
	val						= publicmin
	ocg2rad 				= Array{Float64}(binom(m,2)+m) #the plus 1 covers values taken from the diagonal
	ocg2rad[1]				= val
	post 					= 1
	exceededmax 			= false
	ocf						= fill(-1,m,m) #will reverse the order on this after it's been filled
	stepcounter = 1
	for i = 1:length(p)
		if symmat[p[i]] <= val
			ocf[p[i]] = post
		else
			if numrad == Inf
				val = symmat[p[i]]
			else
				if symmat[p[i]] == Inf
					val = Inf
				else
					while symmat[p[i]] > val
						stepcounter+=1
						if stepcounter == numrad
							val = publicmax # must take this rather cumbersome final step b/c of numerical error, since o/w it can and does happen that the (numrad)th grain value fails to equal publicmax
						else
							val+=alpha
						end
					end
				end
			end
			post+=1
			ocf[p[i]] 		= post
			ocg2rad[post]	= val
		end
		if val > privatemax
			ocf[p[i:end]] 	= post
			exceededmax		= true
			break
		end
	end
	if exceededmax
		cutoff = post
	else
		cutoff = post+1
	end
	deleteat!(ocg2rad,cutoff:length(ocg2rad))
	ocg2rad = flipdim(ocg2rad,1)
	ocf = cutoff - ocf
	return ocf,ocg2rad # additional outputs for diagnostic purposes -> #,privatemax,S,maxrad,publicmax,publicmin,maxrad
end

function ordercanonicalform_4(
	d;
	minrad=-Inf,
	maxrad=Inf,
	numrad=Inf, # note we have already performed the necessary rounding by this point
	fastop=true,
	vscale="diagonal",
	verbose = verbose)

	# round as necessary
	if 	minrad 		== 	"minedge"
		minrad 		= 	minimum(offdiagmin(d))
	end
	d 				= 	minmaxceil(d,minrad=minrad,maxrad=maxrad,numrad=numrad)
	d[d.>maxrad] 	= 	maxrad+1; # we make this reassignment b/c ordercanonicalform_3 takes only arguments with finite entries

	(t,ocg2rad) = ordercanonicalform_3(
		d;
		minrad=minrad,
		maxrad=maxrad,
		numrad=Inf, # note we have already performed the necessary rounding by this point
		fastop=fastop,
		vscale=vscale,
		verbose = verbose)
	return t, ocg2rad
end

# Under development as of 12/30/2017
#
# function 	graduate(
# 			A;
# 			minval		=	minimum(A),
# 			maxval		= 	Inf,
# 			numval		= 	Inf,
# 			stepsize	= 	[],
# 			privatemax	= 	Inf)
#
# 	if numval == 1
# 		return ones(Int64,size(A)...),minimum(A)
# 	end
#
# 	#	Note that maxval and numval cannot both be specified by the user
# 	if stepsize 	= 	[]
# 			alpha	=	(maxval-minval)/(numval-1)
# 		end
# 	else
# 		alpha 		= 	stepsize
# 		if  (maxval 	== 	Inf) 	||	(numval == [])
# 			numval	= 	Inf
# 		else
# 			maxval 	= 	minval+(numval-1)*alpha
# 		end
# 	end
#
# 	p 						= sortperm(vec(A))
#
# 	# Compute the ocf
# 	val						= minval
# 	ocg2rad 				= Array{Float64}(size(A,1)*size(A,2)) #the plus 1 covers values taken from the diagonal
# 	ocg2rad[1]				= val
# 	post 					= 1
# 	exceededmax 			= false
# 	ocf						= fill(-1,m,m) #will reverse the order on this after it's been filled
# 	stepcounter = 1
# 	for i = 1:length(p)
# 		if A[p[i]] <= val
# 			ocf[p[i]] = post
# 		else
# 			if numval == Inf
# 				val = A[p[i]]
# 			else
# 				if A[p[i]] == Inf
# 					val = Inf
# 				else
# 					while A[p[i]] > val
# 						stepcounter+=1
# 						if stepcounter == numval
# 							val = maxval # must take this rather cumbersome final step b/c of numerical error, since o/w it can and does happen that the (numval)th grain value fails to equal maxval
# 						else
# 							val+=alpha
# 						end
# 					end
# 				end
# 			end
# 			post+=1
# 			ocf[p[i]] 		= post
# 			ocg2rad[post]	= val
# 		end
# 		if val > privatemax
# 			ocf[p[i:end]] 	= post
# 			exceededmax		= true
# 			break
# 		end
# 	end
# 	if exceededmax
# 		cutoff = post
# 	else
# 		cutoff = post+1
# 	end
# 	deleteat!(ocg2rad,cutoff:length(ocg2rad))
# 	ocg2rad = flipdim(ocg2rad,1)
# 	ocf = cutoff - ocf
# end

function ordercanonicalform{Tv}(
	S::Array{Tv};
	minrad=-Inf,
	maxrad=Inf,
	numrad=Inf,
	fastop::Bool=true,
	verbose::Bool=false)

    symmat_float = convert(Array{Float64},copy(S))
    symmat = copy(S)
	m = size(symmat,1);
	convert(Tv,minrad)
	convert(Tv,maxrad)

	effectivemin = -maxrad
	effectivemax = -minrad
	symmat = -symmat

	for i = 1:m
		symmat[i,i]=Inf
	end
	if fastop
		maxmin = -Inf
		for j=1:m
			holdmin = minimum(symmat[:,j])
			if holdmin > maxmin
				maxmin = holdmin
			end
		end
		if maxmin > effectivemin
			effectivemin = maxmin
		end
	end
	if numrad == 1
		for i = 1:m
			for j = (i+1):m
				sij = symmat[i,j]
				if sij>effectivemin
					symmat[i,j]=1
					symmat[j,i]=1
				else
					symmat[i,j] = 0
					symmat[j,i] = 0
				end
			end
		end
		for i = 1:m
			symmat[i,i] = 0
		end
		ocg2rad = [1]
		return round.(Int65,symmat),ocg2rad
	end
	numfilt = binom(m,2)
	for i = 1:m
		for j = (i+1):m
			sij = symmat[i,j]
			if sij >= effectivemax
				symmat[i,j] = Inf
				symmat[j,i] = Inf
			elseif sij < effectivemin
				# note that loose inequality
				# here could cause some bars
				# that disappear at the last
				# grain to appear to
				# live forever
				symmat[i,j] = -Inf
				symmat[j,i] = -Inf
			end
		end
	end
	for i = 1:m
		symmat[i,i] = -Inf
	end
	ocg2rad = zeros(Float64,binom(m,2))
	p = sortperm(symmat[:],alg=MergeSort)
	if verbose
		print("done sorting")
	end
	ordervalue = 0
	floatingvalue = symmat[p[1]]
	for i = 1:m^2
		ii = p[i]
		if symmat[ii] == floatingvalue
			symmat[ii] = ordervalue
		else
			ordervalue+=1
			floatingvalue = symmat[ii]
			ocg2rad[ordervalue]=symmat_float[ii]
			symmat[ii]=ordervalue
		end
	end
	deleteat!(ocg2rad,(ordervalue+1):binom(m,2))
	return round.(Int64,symmat),ocg2rad
end

function trueordercanonicalform(	M;
									version		=	2,
									rev			=	false,
									firstval	=	1,
									factor		= 	false )
	if version 				== 	1
		# NB: this version only defined for direction == "up"
		m 					= 	length(M)
		n 					= 	length(unique(M))
		ocf 				= 	zeros(Int64,size(M)...)
		val 				= 	zeros(Float64,n)
		p 					= 	sortperm(M[:],alg=MergeSort)
		k 					= 	1
		val[k] 				= 	M[p[1]]
		ocf[p[1]] 			= 	k
		for i 				= 	2:m
			if 	M[p[i]]		> 	M[p[i-1]]
				k 				+=1
				val[k] 		= 	M[p[i]]
			end
			ocf[p[i]] 	= 	k
		end
	elseif version 			== 	2

		m 					= 	length(M)
		perm 				= 	sortperm(M[:],rev = rev,alg=MergeSort)
		oca 				= 	zeros(Int64,size(M)...)

		if m 				==	0
			return 				zeros(Int64,0),zeros(Int64,0)
		end

		if factor
			numvals 			= 	1
			post 				= 	1
			for p 				= 	1:m
				if M[perm[p]] 	!= 	M[perm[post]]
					post 		= 	p
					numvals 	+= 	1
				end
			end
			oca2rad 			= 	Array{Float64}(numvals)
			oca2rad[1] 			= 	M[perm[1]]
		end

		post 					= 	[1]
		k 						= 	[1]
		trueordercanonicalform_subr!(M,perm,oca,oca2rad,post,k,m,factor)
		if factor
			return oca,oca2rad
		else
			return oca
		end
	end
end

function trueordercanonicalform_subr!(M,perm,oca,oca2rad,post,k,m,factor)
	for p 	 			= 	1:m
		if M[perm[p]] 	!= 	M[perm[post[1]]]
			post[1] 	= 	p
			k[1] 		=  	k[1]+1
			if factor
				oca2rad[k[1]] 	= 	M[perm[post[1]]]
			end
		end
		oca[perm[p]] 	= 	k[1]
	end
end

function checktrueordercanonicalform(numits)
	for p 	= 	1:numits
		m 	= 	rand(50:300,1)
		m 	= 	m[1]
		x 	= 	rand(m,m)
		if 		isodd(p)
			x 	= 	x+x';
		end
		ocf,val 	= 	trueordercanonicalform(x,factor=true)
		check1 		= 	x == val[ocf]
		check2 		= 	val == sort(unique(x))

		k 			= 	rand(1:100,1)
		k 			= 	k[1]
		ocau,valu	=   trueordercanonicalform(x,firstval=1,factor=true,rev=false)
		ocad,vald	=   trueordercanonicalform(x,firstval=1,factor=true,rev=true)
		check3 		= 	ocau == maximum(ocad)+1-ocad
		check4 		= 	x == valu[ocau]
		check5 		= 	x == vald[ocad]
		check6 		= 	length(valu) == maximum(ocau)
		check7 		= 	length(vald) == maximum(ocad)

		if !all([check1, check2, check3, check4, check5, check6, check7])
			println("error: please check trueordercanonicalform")
			println([check1, check2, check3, check4, check5, check6, check7])
			return x
		end
	end
	return []
end


################################################################################
	# 	The following two functions appear to be unused as of 2018-04-15.

function fv2ocff_1(	fv=fv,
					dp=dp,
					minrad=minrad,
					maxrad=maxrad,
					numrad=numrad)

	complexdim 		= 	p -> length(fv[p])
	maxsd 		= 	maxdim+2

	### Format the grain data
	# Concatenate the grain vectors
	numcells 		= 	0
	for p 			= 	1:maxsd
	numcells   += 	complexdim(p)
	end
	filt1 = Array{Float64}(numcells)
	numcells = 0
	for p 								= 	1:maxsd
	l 								= 	complexdim(p)
	filt1[numcells+1:numcells+l] 	= 	fv[p]
	numcells   					   += 	l
	end

	# Convert to order canonical form
	filt2 = integersinoppositeorder_nonunique(filt1)
	ocff = fill(Array{Int64}(0),maxsd+3)
	numcells = 0
	for i = 1:maxsd
	l = length(fv[i])
	ocff[i] = filt2[numcells+1:numcells+l]
	numcells += l
	end

	# Compute the grain translator
	ocg2rad = sort(unique(filt1))
	ocg2rad = flipdim(ocg2rad,1)
end

function fv2ocff_2(	fv=fv,
					dp=dp,
					minrad=minrad,
					maxrad=maxrad,
					numrad=numrad)

	if 	typeof(fv) <= Array{Float64}
		fvo 	= 	copy(fv)
	else
		fvo 	= 	cat(1,fv...)
	end

	fvo 	= 	minmaxceil!( 	fvo,
								minrad=minrad,
								maxrad=maxrad,
								numrad=numrad)

	ocg,ocg2rad 	= 	trueordercanonicalform(fvo,factor=true)
	ocg 			= 	maximum(ocg)+1-ocg
	ocg2rad			= 	flipdim(ocg2rad,1)

end
################################################################################

function checkoffdiagmean(numits)
	for p = 1:numits
		numpts 	= 	rand(50:100,1)
		numpts 	= 	numpts[1]
		S 		= 	rand(numpts,numpts)
		T 		= 	copy(S)
		mu 		= 	offdiagmean(T,defaultvalue = 0)
		u 		= 	zeros(numpts)
		for i 	= 	1:numpts
			ran =   setdiff(1:numpts,i)
			u[i]= 	mean(S[ran,i])
		end
		if 			u[:] != mu[:]
			println()
			println("error: please check <offdiagmean>")
			return S,T,mu,u
		end
	end
	return []
end

function getstarweights(symmat)
	m = size(symmat,1)
	w = zeros(Int64,m)
	getstartweights_subr2(symmat::Array{Int64,2},w::Array{Int64,1},m::Int64)
	return w
end

function getstartweights_subr2(symmat::Array{Int64,2},w::Array{Int64,1},m::Int64)
	s = copy(symmat)
	for i = 1:m
		s[i,i]=0
	end
	l = Array{Int64}(m)
	lDown = Array{Int64}(m)
	val 		= Array{Array{Int64,1}}(m)
	supp 		= Array{Array{Int64,1}}(m)
	suppDown 	= Array{Array{Int64,1}}(m)
	for i = 1:m
		supp[i] = find(s[:,i])
		suppDown[i] = i+find(s[(i+1):end,i])
		val[i] = s[supp[i],i]
		l[i] = length(supp[i])
		lDown[i] = length(suppDown[i])
	end

	for i = 1:m
		Si = supp[i]
		Vi = val[i]
		for jp = 1:l[i]
			j = Si[jp]
			dij = Vi[jp]
			Sj = suppDown[j]
			Vj = val[j]
			for kp = 1:lDown[j]
				k = Sj[kp]
				if k == i
					continue
				end
				dkj = Vj[kp]
				dki = s[i,k]
				if dki >= dkj && dij >= dkj
					w[i]+=1
				end
			end
		end
	end
	return w
end

##########################################################################################

####	COMBINATIONS, PERMUTATIONS, AND SET OPERATIONS

##########################################################################################

function intervalcomplementuniquesortedinput(sortedVecOfUniqueIntegers,intervalEndpoint)
	const v = sortedVecOfUniqueIntegers
	const n = intervalEndpoint
	const L = length(v)
	if L==0
		return 1:n
	elseif L==n
		return Array{Int64}(0)
	else
		boundMarker = 1
		upperBound = v[boundMarker]
		complement = Array{Int64}(n-L)
		marker = 0
		for i = 1:(v[end]-1)
			if i<upperBound
				marker+=1
				complement[marker] = i
			else
				boundMarker+=1
				upperBound = v[boundMarker]
			end
		end
		complement[(marker+1):end]=(v[end]+1):n
	end
	return complement
end

function intervalcomplementuniqueunsortedinput(uniquepositiveintegers,intervalEndpoint)
	const v = uniquepositiveintegers
	const n = intervalEndpoint
	const L = length(v)
	if L==0
		return 1:n
	elseif L==n
		return Array{Int64}(0)
	else
		complementsupport = trues(n)
		complementsupport[v]=false
		complement = Array{Int64}(n-L)
		marker = 0
		for i = 1:n
			if complementsupport[i]
				marker+=1
				complement[marker]=i
			end
		end
	end
	return complement
end

function integersinsameorder!(v::Array{Int64,1},maxradue::Int64)
	m = length(v)
	x = zeros(Int64,maxradue)
	for i = 1:m
		x[v[i]]+=1
	end
	y = Array{Int64}(maxradue+1)
	y[1] = 1
	for i = 1:maxradue
		y[i+1]=y[i]+x[i]
	end
	for i = 1:length(v)
		u = v[i]
		v[i] = y[u]
		y[u]+=1
	end
	return v
end

function integersinsameorder(v::Array{Int64,1})
	# Returns the permutation z on {1,...,length(v)} such z[i]<z[j] iff either
	# (a) v[i] < v[j], or
	# (b) v[i] = v[j] and i < j
	if isempty(v)
		z = Array{Int64}(0)
		return z
	else
		m = length(v)
		maxv = maximum(v)
		minv = minimum(v)
		minv = minv-1;
		x = zeros(Int64,maxv-minv)
		z = Array{Int64}(length(v))
		for i = 1:m
			x[v[i]-minv]+=1
		end
		prevsum = 1
		for i = 1:length(x)
			sum = prevsum + x[i]
			x[i] = prevsum
			prevsum = sum
		end
		for i = 1:m
			u = v[i]
			z[i] = x[u-minv]
			x[u-minv]+=1
		end
		return z
	end
end

function integersinsameorder!(v::Array{Int64,1})
	# Replaces v with the permutation z on {1,...,length(v)} such that z[i]<z[j] iff either
	# (a) v[i] < v[j], or
	# (b) v[i] = v[j] and i < j
	if isempty(v)
		z = Array{Int64}(0)
		return z
	else
		m = length(v)
		maxv = maximum(v)
		minv = minimum(v)
		minv = minv-1;
		x = zeros(Int64,maxv-minv)
		z = Array{Int64}(length(v))
		for i = 1:m
			x[v[i]-minv]+=1
		end
		prevsum = 1
		for i = 1:length(x)
			sum = prevsum + x[i]
			x[i] = prevsum
			prevsum = sum
		end
		for i = 1:m
			u = v[i]
			v[i] = x[u-minv]
			x[u-minv]+=1
		end
	end
end

function integersinoppositeorder_nonunique(v)
	if isempty(v)
		return Array{Int64,1}(0)
	end
	p 			= sortperm(v,alg=MergeSort)
	u 			= Array{Int64}(length(v))
	epsilon 	= v[p[end]]
	c			= 1
	for i = length(v):-1:1
		if v[p[i]] == epsilon
			u[p[i]] = c
		else
			c+=1
			u[p[i]] = c
			epsilon = v[p[i]]
		end
	end
	return u
end

function integersinoppositeorder_nonunique_test()
	## diagnostic test to ensure integersinoppositeorder_nonunique worksp properly
	for i = 1:100
		v = rand(10000)
		u = unique(v)
		u = sort(u)
		u = flipdim(u,1)
		l = length(v)
		w = Array{Int64}(l)
		for i = 1:l
			w[i] = findfirst(u,v[i])
		end
		if w!=integersinoppositeorder_nonunique(v)
			print("hmmm")
			return
		end
	end
end

function integersinsameorderbycolumn(v::Array{Int64,1},maxradue::Int64,colptr)
	# Returns a permutation z on {1,...,length(v)} so that
	# (a) cran(colptr,j) maps to cran(colptr,j) for all j, and
	# (b) crows(colptr,v[z],j) is an array in sorted order
	numcols = length(colptr)-1
	m = length(v)
	x = Array{Int64}(maxradue)
	y = Array{Int64}(maxradue+1)
	z = Array{Int64}(length(v))
	for j = 1:numcols
		x[:] = 0
		for i = colptr[j]:(colptr[j+1]-1)
			x[v[i]]+=1
		end
		y[1] = colptr[j]
		for i = 1:maxradue
			y[i+1]=y[i]+x[i]
		end
		for i = colptr[j]:(colptr[j+1]-1)
			u = v[i]
			z[i] = y[u]
			y[u]+=1
		end
	end
	return z
end

#=
- 	In beta; should be compared with integersinsameorderbycolumn3.  See
	/Users/greghenselman/Google Drive/GregDirectory/julia_gd/Julia/workshop/workshop_Oct2017.jl
-   Functionally equivalent to integersinsameorderbycolumn; returns a
	permutation z on {1,...,length(v)} so that for all j
	- cran(colptr,j) maps to cran(colptr,j), and
	- crows(colptr,v[z],j) is an array in sorted order
=#
function integersinsameorderbycolumn2(v::Array{Int64,1},colptr)
	numcols = length(colptr)-1
	m = length(v)
	v = v-minimum(v)+1
	x = zeros(Int64,maximum(v))
	z = Array{Int64}(length(v))
	for j = 1:numcols
		if colptr[j] == colptr[j+1]
			continue
		end
		for i = colptr[j]:(colptr[j+1]-1)
			x[v[i]]+=1
		end
		maxv = v[colptr[j]];   minv = maxv
		for i = (colptr[j]+1):(colptr[j+1]-1)
			if v[i] > maxv
			   maxv = v[i]
			elseif v[i] < minv
			   minv = v[i]
			end
		end
		prevsum = colptr[j]
		for i = minv:maxv
			sum = prevsum + x[i]
			x[i] = prevsum
			prevsum = sum
		end
		for i = colptr[j]:(colptr[j+1]-1)
			u = v[i]
			z[i] = x[u]
			x[u]+=1
		end
		for i = minv:maxv
			x[i] = 0
		end
	end
	return z
end

#=
- 	In beta; should compare with integersinsameorderbycolumn2. See
	/Users/greghenselman/Google Drive/GregDirectory/julia_gd/Julia/workshop/workshop_Oct2017.jl
- 	Functionally equivalent to integersinsameorderbycolumn; returns a
	permutation z on {1,...,length(v)} so that
	- cran(colptr,j) maps to cran(colptr,j) for all j, and
	- crows(colptr,v[z],j) is an array in sorted order
=#
function integersinsameorderbycolumn3(v::Array{Int64,1},colptr)
	numcols = length(colptr)-1
	z 	    = Array{Int64}(length(v))
	for i = 1:numcols
		z[cran(colptr,i)] = colptr[i]-1+integersinsameorder(crows(colptr,v,i))
	end
	return z
end

function integersortperm(v::Array{Int64,1},maxradue::Int64)
	l = length(v)
	u = integersinsameorder(v)
	w = Array{Int64}(l)
	for i = 1:l
		w[u[i]] = i
	end
	return w
end

##########################################################################################

####	SEARCH SUBROUTINES

##########################################################################################

function getPairsLightWrite2!{Tv<:Integer}(
	rowval::Array{Tv,1},
	colptr::Array{Tv,1},
	rowfilt::Array{Tv,1},
	colfilt::Array{Tv,1},
	m::Integer,
	n::Integer,
	prows::Array{Tv,1},
	pcols::Array{Tv,1},
	numpairs::Array{Tv,1};
	verbose = false)

	col2firstplace = zeros(Tv,n)
	rowwisesum = zeros(Tv,m)

	if verbose
		println("starting search for pairs")
	end
	for j = 1:n
		firstplace = colptr[j]
		if firstplace < colptr[j+1]
			firstrow = rowval[firstplace]
			rowwisesum[firstrow]+=1
			if firstplace < colptr[j+1]-1
				filt = rowfilt[firstrow]
				for newplace = (firstplace+1):(colptr[j+1]-1)
					newrow = rowval[newplace]
					newfilt = rowfilt[newrow]
					rowwisesum[newrow]+=1
					if newfilt > filt
						filt = newfilt
						firstplace = newplace
					end
				end
			end
			col2firstplace[j]=firstplace
		end
	end
	colwisesum = colsupportsum(colptr,n) # note allow extra on end
	colfiltptr = getcolptr2(colfilt,n)	# note allow extra on end
	colwisesumlinearized = integersinsameorderbycolumn2(colwisesum,colfiltptr)
	colnamesinorder = Array{Tv}(n)
	colnamesinorder[colwisesumlinearized]=1:n
	ncoveredsupp = trues(m)
	pairmarker = 0
	for jp = 1:n
		j = colnamesinorder[jp]
		if col2firstplace[j]>0
			firstplace = col2firstplace[j]
			firstrow = rowval[firstplace]
			if firstplace==colptr[j+1]-1
				if ncoveredsupp[firstrow]	#######&&  (pairmarker==0 || rowwisesum[firstrow]<10000)
					pairmarker+=1
					prows[pairmarker]=firstrow
					pcols[pairmarker]=j
				end
			else
				firstweight = rowwisesum[firstrow]
				filt = rowfilt[firstrow]
				for newplace = (firstplace+1):(colptr[j+1]-1)
					newrow = rowval[newplace]
					newweight = rowwisesum[newrow]
					if rowfilt[newrow]==filt && newweight <= firstweight && !ncoveredsupp[firstrow] && ncoveredsupp[newrow]
						firstweight = newweight
						firstrow = newrow
					end
				end
				if ncoveredsupp[firstrow] ######&&  (pairmarker==0 || rowwisesum[firstrow]<10000)
					pairmarker+=1
					prows[pairmarker]=firstrow
					pcols[pairmarker]=j
				end
			end
			for ip = cran(colptr,j)
				ncoveredsupp[rowval[ip]]=false
			end
		end
	end

	numpairs[1]=pairmarker
	if verbose
		print("done finding pairs")
	end
end

function finddownstreamelements_embeddedupperunitriangularmatrix(
	Mrv,
	Mcp,
	Mm,
	initialelements::Array{Int64,1},
	prows::Array{Int64,1},
	pcols::Array{Int64,1};
	verbose=false)

	if verbose
		print("starting to get downstream elements")
	end
	if length(prows) != length(pcols)
		print("length of p doesn't match length of q")
		return
	elseif length(prows)==0
		return Array{Int64}(0)
	end
	n = length(prows)
	rowtranslator = Array{Int64}(Mm)
	for i = 1:n
		rowtranslator[prows[i]]=i
	end
	prowsupp = falses(Mm)
	prowsupp[prows]=true
	downstreamsupport = falses(n)
	for i = 1:length(initialelements)
		row = initialelements[i]
		if prowsupp[row]
			downstreamsupport[rowtranslator[row]] = true
		end
	end
	for jp = n:-1:1
		j = pcols[jp]
		ran = cran(Mcp,j)
		for ip in ran
			rawrow = Mrv[ip]
			if prowsupp[rawrow] && downstreamsupport[rowtranslator[rawrow]]
				for kp in ran
					rawrow = Mrv[kp]
					if prowsupp[rawrow]
						downstreamsupport[rowtranslator[rawrow]] = true
					end
				end
				break
			end
		end
	end
	counter = 0
	for i = 1:n
		if downstreamsupport[i]
			counter+=1
		end
	end
	downstreamelements = Array{Int64}(counter)
	counter = 0
	for i = 1:n
		if downstreamsupport[i]
			counter+=1
			downstreamelements[counter] = i
		end
	end
	if verbose
		print("done with downstream elements")
	end
	return downstreamelements
end

##########################################################################################

####	BARCODE UTILITIES

##########################################################################################


function nnzbars(D::Dict;dim = 1)
	sd = dim+2
	if !(0<sd<=length(D["plo"]))
		print("error: requested dimension is outside the range of calculated bars")
		return
	end
	plo = D["plo"][sd]
	phi = D["phi"][sd]
	lowfilt = D["grain"][sd-1]
	higfilt = D["grain"][sd]
	counter::Int64 = 0
	for i = 1:length(plo)
		if higfilt[phi[i]]!=lowfilt[plo[i]]
			counter+=1
		end
	end
	counter += length(D["tid"][sd])-length(plo)
	return counter
end

function nnzbars_test()
	for p = 1:20
		C 	= 	eirene(rand(20,50),model="pc",maxdim=2)
		cr 	= 	C["cyclerep"]
		for q 	=	1:3

		end
	end
end

function barname2cyclename(D::Dict,barnumber = [1];dim = 1)
	if typeof(barnumber) <: Array
		sd = dim+2
		tid = D["tid"][sd]
		plo = D["plo"][sd]
		phi = D["phi"][sd]
		nummortals = length(plo)
		nzcycles = find(D["grain"][sd][phi].!= D["grain"][sd-1][plo])
		append!(nzcycles,Array((nummortals+1):length(tid)))
		return nzcycles[barnumber]
	elseif typeof(barnumber) <: UnitRange
		sd = dim+2
		tid = D["tid"][sd]
		plo = D["plo"][sd]
		phi = D["phi"][sd]
		nummortals = length(plo)
		nzcycles = find(D["grain"][sd][phi].!= D["grain"][sd-1][plo])
		append!(nzcycles,nummortals+1:length(tid))
		return nzcycles[barnumber]
	elseif typeof(barnumber)<:Number
		sd = dim+2
		tid = D["tid"][sd]
		plo = D["plo"][sd]
		phi = D["phi"][sd]
		numclasses = length(tid)
		nummortals = length(plo)
		counter = 0
		cyclename = 0
		for i = 1:nummortals
			if D["grain"][sd][phi[i]] != D["grain"][sd-1][plo[i]]
				counter+=1
			end
			if counter == barnumber
				cyclename = i
				break
			end
		end
		if cyclename == 0
			for i = (nummortals+1):length(tid)
				counter+=1
				if counter == barnumber
					cyclename = i
					break
				end
			end
		end
		return cyclename
	end
end

function getbetticurve(D::Dict,sd;ocf = false)
	# NB: D["ocg2rad"]: {grains > 0} --> rad
	# NB: the ocf barcode takes values in [0, #{grains}-1]

	if length(D["farfaces"][sd])==0
		return Array{Float64}(0,2)
	end

	# the plus 1 is for the zero grain, the "s" in ngrains is for "shifted",
	# since we shift up to avoid zero indices
	ngrains = length(D["ocg2rad"])
	v = zeros(Int64,ngrains)

	bco = barcode(D;dim = sd-2,ocf=true)
	# bco[bco.==Inf] = ngrains
	bco = convert(Array{Int64},bco)

	for i = 1:size(bco,1)
		ran = 1+ (bco[i,1]:(bco[i,2]-1))
		v[ran]+=1
	end

	if ocf == false
		u = sort(D["ocg2rad"])
	else
		u = Array(1:ngrains)
	end
	return hcat(u,v)
end

function getcycle(farfaces,firstv,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plo,phi,tid,sd,cyclenumber::Number)
	if sd == 2
		numlowlows = 0
	else
		numlowlows = length(farfaces[sd-2])
	end
	numnlpl = length(farfaces[sd-1])-length(plo[sd-1])

	summands = tid[sd][crows(Licp[sd],Lirv[sd],cyclenumber)]
	append!(summands,[tid[sd][cyclenumber]])
	brv = ff2aflight(farfaces,firstv,sd-1,summands)
	supp = falses(length(farfaces[sd-2]))
	for k in brv
		supp[k] = !supp[k]
	end

	brv = find(supp[tid[sd-1]])
	bcp = [1,length(brv)+1]
	brv,bcp = spmmF2silentLeft(Lrv[sd-1],Lcp[sd-1],brv,bcp,numnlpl)
	brv,bcp = spmmF2silentLeft(Rrv[sd-1],Rcp[sd-1],brv,bcp,numnlpl)

	plow2phigtranslator = zeros(Int64,numlowlows)
	plow2phigtranslator[plo[sd-1]]=phi[sd-1]
	brv		= plow2phigtranslator[tid[sd-1][brv]] # some of the nonzero entries might lie in non-basis rows
	brv 	= brv[find(brv)]

	brv = append!(brv,summands)

	return brv
end

function getcycle_cell(rv,cp,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plo,phi,tid,sd,cyclenumber::Number)
	if sd == 2
		# this case must be treated separately b/c julia indexing starts at 1
		numlowlows = 0
		numnlpll = 0
	else
		numlowlows = length(cp[sd-2])-1
		numnlpll = numlowlows-length(plo[sd-2])
	end
	numnlpl = length(cp[sd-1])-1-length(plo[sd-1])	# the -1 accounts for the fact that length(cp[sd-1]) = (# cells of dimension secard-1) - 1

	summands = tid[sd][crows(Licp[sd],Lirv[sd],cyclenumber)]
	append!(summands,[tid[sd][cyclenumber]])

	if sd == 2
		return summands
	end

	supp = falses(numlowlows)
	sc = sd-1
	for j in summands	# this is a bit ridiculus; it's here bc i don't have a verified form of spmmF2 on hand
		for k in cran(cp[sc],j)
			i = rv[sc][k]
			supp[i] = !supp[i]
		end
	end
	brv 				= find(supp[tid[sd-1]])
	bcp 				= [1,length(brv)+1]
	brv,bcp 			= spmmF2silentLeft(Lrv[sd-1],Lcp[sd-1],brv,bcp,numnlpll)
	brv,bcp 			= spmmF2silentLeft(Rrv[sd-1],Rcp[sd-1],brv,bcp,numnlpll)
	plow2phigtranslator = zeros(Int64,numlowlows)
	plow2phigtranslator[plo[sd-1]]=phi[sd-1]
	brv					= plow2phigtranslator[tid[sd-1][brv]] # some of the nonzero entries might lie in non-basis rows
	brv 				= append!(brv,summands)
	return brv
end

function getcycle_cell(rv,cp,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plo,phi,tid,sd,cyclenumber::Array{Int64,1})
	numclasses 	= length(cyclenumber)
	rep 	 	= Array{Array{Int64,1},1}(numclasses)
	for p 		= 1:numclasses
		rep[p] 	= getcycle_cell(rv,cp,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plo,phi,tid,sd,cyclenumber[p])
	end
	return rep
end

function getcycle(farfaces,firstv,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plo,phi,tid,sd,cyclenumber)
	if sd == 2
		numlowlows = 0
	else
		numlowlows = length(farfaces[sd-2])
	end
	numlows    = length(farfaces[sd-1])
	numnlpl = length(farfaces[sd-1])-length(plo[sd-1])

	numclasses = length(cyclenumber)
	summands = Array{Array{Int64,1},1}(numclasses)
	rep 	 = Array{Array{Int64,1},1}(numclasses)
	summandsupp = falses(numlows)
	for i = 1:numclasses
		summands[i] = tid[sd][crows(Licp[sd],Lirv[sd],cyclenumber[i])]
		append!(summands[i],[tid[sd][cyclenumber[i]]])
		summandsupp[summands[i]]=true
	end

	lowgenerators = find(summandsupp)
	numlowgenerators = length(lowgenerators)
	translator = zeros(Int64,numlows)
	translator[lowgenerators] = 1:length(lowgenerators)

	lowfacemat = ff2aflight(farfaces,firstv,sd-1,lowgenerators)

	supp = falses(numlowlows)
	m = size(lowfacemat,1)
	plow2phigtranslator = Array{Int64}(numlowlows)
	plow2phigtranslator[plo[sd-1]]=phi[sd-1]
	for i = 1:numclasses

		supp[:] = false
		for j = 1:length(summands[i])
			for k = 1:m
				kk = lowfacemat[k,translator[summands[i][j]]]
				supp[kk] = !supp[kk]
			end
		end

		brv = find(supp[tid[sd-1]])
		bcp = [1,length(brv)+1]
		brv,bcp = spmmF2silentLeft(Lrv[sd-1],Lcp[sd-1],brv,bcp,numlowlows)
		brv,bcp = spmmF2silentLeft(Rrv[sd-1],Rcp[sd-1],brv,bcp,numlowlows)

		rep[i] = append!(plow2phigtranslator[tid[sd-1][brv]],summands[i])
	end
	return rep
end

function getcycle(D::Dict,sd,cyclenumber)
	if !haskey(D,"Lirv")
		if !haskey(D,"cyclerep")
			println("This object does not store a complete cycle basis.")
		else
			println("This object does not store a complete cycle basis, only those cycles that represent persistent homology classes.")
		end
		return
	end
	farfaces = D["farfaces"];firstv = D["firstv"];Lirv = D["Lirv"];Licp=D["Licp"];Lrv=D["Lrv"]
	Lcp=D["Lcp"];Rrv=D["Rrv"];Rcp=D["Rcp"];plo=D["plo"];phi=D["phi"];tid=D["tid"]
	rrv = getcycle(farfaces,firstv,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plo,phi,tid,sd,cyclenumber)
	return rrv
end

function getcycle(D::Dict,cyclenumber;dim = 1)
	if !haskey(D,"Lirv")
		println("Error: the function <getcycle(D::Dict,cyclenumber;dim = 1)> assumes a key value for \"Lirv\" in the input object D.  This key value is absent.")
		return
	end
	sd = dim+2
	Lirv = D["Lirv"];Licp=D["Licp"];Lrv=D["Lrv"]
	Lcp=D["Lcp"];Rrv=D["Rrv"];Rcp=D["Rcp"];
	plo=D["plo"];phi=D["phi"];tid=D["tid"]
	if haskey(D,"farfaces")
		farfaces = D["farfaces"];firstv = D["firstv"];
		rrv = getcycle(farfaces,firstv,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plo,phi,tid,sd,cyclenumber)
	else
		rrv = getcycle_cell(D["rv"],D["cp"],Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plo,phi,tid,sd,cyclenumber)
	end
	return rrv
end

function getcyclesize(farfaces,firstv,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plo,phi,tid,sd,cyclenumber)
	if sd == 2
		numlowlows = 0
	else
		numlowlows = length(farfaces[sd-2])
	end
	numlows    = length(farfaces[sd-1])
	numnlpl = length(farfaces[sd-1])-length(plo[sd-1])

	numclasses = length(cyclenumber)
	summands = Array{Array{Int64,1},1}(numclasses)
	rep 	 = Array{Int64}(numclasses)
	summandsupp = falses(numlows)
	for i = 1:numclasses
		summands[i] = tid[sd][crows(Licp[sd],Lirv[sd],cyclenumber[i])]
		append!(summands[i],[tid[sd][cyclenumber[i]]])
		summandsupp[summands[i]]=true
	end

	lowgenerators = find(summandsupp)
	numlowgenerators = length(lowgenerators)
	translator = zeros(Int64,numlows)
	translator[lowgenerators] = 1:length(lowgenerators)

	lowfacemat = ff2aflight(farfaces,firstv,sd-1,lowgenerators)

	supp = falses(numlowlows)
	m = size(lowfacemat,1)
	plow2phigtranslator = Array{Int64}(numlowlows)
	plow2phigtranslator[plo[sd-1]]=phi[sd-1]
	for i = 1:numclasses

		supp[:] = false
		for j = 1:length(summands[i])
			for k = 1:m
				kk = lowfacemat[k,translator[summands[i][j]]]
				supp[kk] = !supp[kk]
			end
		end

		brv = find(supp[tid[sd-1]])
		bcp = [1,length(brv)+1]
		brv,bcp = spmmF2silentLeft(Lrv[sd-1],Lcp[sd-1],brv,bcp,numnlpl)
		brv,bcp = spmmF2silentLeft(Rrv[sd-1],Rcp[sd-1],brv,bcp,numnlpl)

		rep[i] = length(brv)+length(summands[i])
	end
	return rep
end

function getcyclesize(D::Dict,cyclenumber;dim = 1)
	sd = dim+2
	if !haskey(D,"Lirv")
		if !haskey(D,"cyclerep")
			println("This object does not store a complete cycle basis.")
		else
			println("This object does not store a complete cycle basis, only those cycles that represent persistent homology classes.")
		end
		return
	end
	farfaces = D["farfaces"];firstv = D["firstv"];Lirv = D["Lirv"];Licp=D["Licp"];Lrv=D["Lrv"]
	Lcp=D["Lcp"];Rrv=D["Rrv"];Rcp=D["Rcp"];plo=D["plo"];phi=D["phi"];tid=D["tid"]

	rrv = getcyclesize(farfaces,firstv,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plo,phi,tid,sd,cyclenumber)
	return rrv
end

function getrepsize(D::Dict,classnumber;dim=1)
	sd = dim+2
	if !haskey(D,"cyclerep")
		println("This object does not contain data about cycle representatives.")
		return
	elseif typeof(classnumber)<: Number
		return length(D["cyclerep"][dim+2][classnumber])
	else
		l = length(classnumber)
		rsize = Array{Int64}(l)
		for i = 1:l
			rsize[i] = length(D["cyclerep"][dim+2][classnumber[i]])
		end
		return rsize
	end
end

function birthtime(C;chain=zeros(Int64,0),dim=1)
	if isempty(chain)
		return -Inf
	else
		translator 	= 	C["ocg2rad"]
		sd 	= 	dim+1
		ocg 		= 	C["grain"][sd] # ocg stands for order-canonical grain
		return 		empteval(maximum,translator[ocg[chain]],-Inf)
	end
end

function deathtime(C;chain=zeros(Int64,0),dim=1)
	if isempty(chain)
		return -Inf
	elseif !isempty(chainboundary(C,chain=chain,dim=dim))
		return Inf
	else
		sd 	= 	dim+1;
		Lrv 		= 	C["Lrv"][sd+1]
		Lcp 		= 	C["Lcp"][sd+1]
		Rrv 		= 	C["Rrv"][sd+1]
		Rcp 		= 	C["Rcp"][sd+1]
		numrows 	= 	boundarycorank(C,dim=dim)
		translator 	= 	zeros(Int64,complexrank(C,dim=dim))
		translator[C["tid"][sd+1]] = 1:numrows
		rv 			= 	find(translator[chain])
		rv 			=  	translator[chain[rv]]
		cp 			= 	[1,length(rv)+1]
		rv,cp 	=  	spmmF2silentLeft(Lrv,Lcp,rv,cp,numrows)
		rv,cp 	=  	spmmF2silentLeft(Rrv,Rcp,rv,cp,numrows)
		#
		# recall that plo = the (ordered) subvector consisting of
		# the first (rank(boundary operator)) elements of tid
		#
		if maximum(rv) > length(C["plo"][sd+1])
			return 		Inf
		else
			ocg2rad = 	C["ocg2rad"]
			grain 	= 	C["grain"][sd+1]
			names 	= 	C["phi"][sd+1]
			return 		empteval(maximum,ocg2rad[grain[names[rv]]],-Inf)
		end
	end
end

function chainboundary(C;chain=zeros(Int64,0),dim=1)
	m 			= 	length(chain)
	if m == 0
		return zeros(Int64,0)
	end
	crv 		= 	convert(Array{Int64,1},1:m);
	ccp 		= 	[1,m+1]
	brv,bcp 	= 	boundarymatrix(C,dim=dim,cols=chain)
	brv,bcp 	= 	spmmF2(brv,bcp,crv,ccp,empteval(maximum,brv,0))
	return 			brv
end

function numcols(cp)
	return length(cp)-1
end

function complexrank(C;dim=1)
	sd 		= 	dim+1
	if 	dim > C["input"]["maxdim"]+1 || dim < 0
		return 0
	elseif C["input"]["model"] == "complex"
		return length(C["cp"][sd])-1
	else
		return length(C["farfaces"][sd])
	end
end

function boundaryrank(C;dim=1)
	sd 	= 	dim+1;
	if complexrank(C,dim=dim) == 0
		return 0
	else
		return length(C["plo"][sd])
	end
end

function boundarycorank(C;dim=1)
	sd 	= 	dim+1;
	if complexrank(C,dim=dim) == 0
		return 0
	else
		return complexrank(C,dim=dim)-boundaryrank(C,dim=dim)
	end
end

function empteval(f,a,c)
	if isempty(a)
		return c
	else
		return f(a)
	end
end

##########################################################################################

####	USER-FRIENDLY UTILITIES

##########################################################################################

function boundarymatrix(C;dim=1,rows="a",cols="a")
	crr 					= 	complexrank(C,dim=dim-1)
	crc 					= 	complexrank(C,dim=dim)
	if rows == "a"
		rows 				= 	1:crr#complexrank(C,dim=dim-1)
	end
	if cols == "a"
		cols 				= 	1:crc#complexrank(C,dim=dim)
	end
	if empteval(maximum,cols,0) > crc
		println()
		println("error: keyword argument <cols> contains an integer greater than the rank of the complex in dimension <dim>")
		println()
		return
	elseif empteval(maximum,rows,0) > crc
		println()
		print("error: keyword argument <rows> contains an integer greater than the rank of the complex in dimension <dim-1>")
		println()
		return
	end
	if isempty(rows) || isempty(cols)
		rv 					= 	zeros(Int64,0)
		cp 					= 	ones(Int64,length(cols)+1)
		return 					rv,cp
	end
	ncols 					= 	length(cols)
	nrows 					= 	length(rows)
	sd 						= 	dim+1;
	if haskey(C,"farfaces")
		rv 					= 	ff2aflight(C,dim+1,cols)
		rv  				= 	reshape(rv,length(rv))
		cp  				= 	convert(Array{Int64,1},1:sd:(ncols*sd+1))
		cols 				= 	1:ncols
	else
		rv 					= 	C["rv"][sd]
		cp 					= 	C["cp"][sd]
	end
	rv,dummrv,cp,dummycp 	= 	stackedsubmatrices(
								rv,
								cp,
								rows,
								Array{Int64}(0),
								cols,
								max(empteval(maximum,rows,0),empteval(maximum,rv,0))
								)
	return 						rv,cp
end

function boundarymatrices(C)
	if haskey(C,"farfaces")
		rv,cp = ff2complex(C["farfaces"],C["firstv"])
	else
		rv = C["rv"]
		cp = C["cp"]
	end
	return rv,cp
end

function classrep(
	D::Dict;
	dim = 1,
	class = 1,
	format = "vertex x simplex")

	if any(class.>nnzbars(D,dim=dim))
		print("error: the value for keyword argument <class> has an integer greater than the number of nonzero bars in the specified dimension")
		return
	elseif !(0<=dim<=D["input"]["maxdim"])
		print("error: barcodes were not computed in the specified dimension")
		return
	end
	if !haskey(D,"farfaces")
		format = "index"
	end

	if format == "vertex x simplex"
		return classrep_faces(D,dim = dim,class = class)
	elseif format == "vertex"
		return classrep_vertices(D,dim = dim, class = class)
	elseif format == "index"
		return classrep_cells(D,dim=dim,class=class)
	end
end

function classrep_cells(
	D::Dict;
	dim = 1,
	class = 1)

	sd = dim+2

	if haskey(D,"cyclerep")
		rep = D["cyclerep"][sd][class]
	else
		cyclename = barname2cyclename(D,class;dim=dim)
		rep = getcycle(D,sd,class)
	end

	return rep
end

function classrep_faces(
	D::Dict;
	dim = 1,
	class = 1)

	sd 				= 	dim+2
	rep 			= 	classrep_cells(D,dim=dim,class=class)

	vrealization 	= 	vertexrealization(D::Dict,sd-1,rep)
	vrealization 	= 	D["nvl2ovl"][vrealization]
	return vrealization
end

function classrep_vertices(
	D::Dict;
	dim = 1,
	class = 1)

	sd 			= 	dim+2
	rep 		= 	classrep_cells(D,dim=dim,class=class)

	vertices 	= 	incidentverts(D::Dict,sd-1,rep)
	vertices 	= 	D["nvl2ovl"][vertices]
	return vertices
end

function cyclevertices(
	D::Dict;
	dim = 1,
	cycle = 1)

	sd 			= dim+2
	rep 		= getcycle(D,sd,cycle)
	vertices 	= incidentverts(D::Dict,sd-1,rep)
	vertices 	= D["nvl2ovl"][vertices]
	return vertices
end

function printsize(var,varname)
	println(string("size(",varname,") = ",size(var)))
end

function printval(var,varname)
	println(string(varname," = $(var)"))
end

function barcode(D::Dict;dim = 1,ocf = false)
	if 	haskey(D,:perseusjlversion)
		return	barcode_perseus(D,dim=dim)
	elseif haskey(D,:barcodes)
		return  D[:barcodes][dim+1]
	elseif !haskey(D,"cp") & !haskey(D,"farfaces")
		print("unrecognized object:")
		display(D)
		return
	elseif dim > D["input"]["maxdim"]
		maxdim 	= 	D["input"]["maxdim"]
		println("error: barcodes were not computed in dimensions greater than $(maxdim).")
		return
	end
	sd = dim+2
	plo = D["plo"][sd]
	phi = D["phi"][sd]
	tid = D["tid"][sd]
	lg 	= D["grain"][sd-1]
	hg  = D["grain"][sd]

	mortalprimagrain 	= 	lg[plo]
	mortalultragrain 	= 	hg[phi]

	finind 				= 	find(mortalprimagrain .!= mortalultragrain)
	numfin 				= 	length(finind)
	numinf 				= 	length(tid)-length(plo)

	mortalprimagrain 	= 	mortalprimagrain[finind]
	mortalultragrain 	= 	mortalultragrain[finind]

	mortalran 			= 	1:numfin
	evergrran 			= 	numfin+1:numfin+numinf
	finran 				= 	1:(2*numfin+numinf)
							# stands for finite range; these are the linear
							# indices of the barcode array that take finite
							# values
	evrgrbran 			= 	length(tid)-numinf+1:length(tid)
							# stands for evergreen birth range; this satisfies
							# tid[ebergrbran] = [array of evergreen classes]

	bc 					= 	zeros(Int64,numfin+numinf,2)
	bc[mortalran,1] 	= 	mortalprimagrain
	bc[mortalran,2] 	= 	mortalultragrain
	bc[evergrran,1] 	= 	lg[tid[evrgrbran]]

	if !ocf
		bcc 				= 	copy(bc)
		bc 					= 	Array{Float64}(bc)
		bc[finran] 			= 	D["ocg2rad"][bcc[finran]]
		bc[evergrran,2] 	= 	Inf
	else
		bc 					= 	length(D["ocg2rad"])-bc
	end

	return bc
end

function getpersistencediagramprimitives(
	C;
	dim = 1,
	ocf = false,
	descriptivetext = true,
	showsize = true)

	if haskey(C,"cyclerep")
		showsize = true
	else
		showsize = false
	end

	bco = barcode(C,dim = dim,ocf = ocf)
	rows = unique(bco,1)
	numrows = size(rows,1)
	numbrs = size(bco,1)
	sd = dim+2

	if numbrs == 0
		x0=[];y0=[];l0=[];x1=[];y1=[];l1=[];x2=[];y2=[]
		return x0,y0,l0,x1,y1,l1,x2,y2
	end

	if showsize
		barsizes = Array{Int64}(numbrs)
		for i = 1:numbrs
			barsizes[i] = length(C["cyclerep"][sd][i])
		end
	end

	if descriptivetext
		if showsize
			labels = fill("class/size  ",numrows)
		else
			labels = fill("class  ",numrows)
		end
	else
		labels = fill("",numrows)
	end

	D = Dict()
	for i = 1:numrows
		D[Symbol(rows[i,:])] = i
	end
	for j = 1:numbrs
		uniquerowind = D[Symbol(bco[j,:])]
		if labels[uniquerowind] in ["", "class  ","class/size  "]
			if showsize
				labels[uniquerowind] = string(labels[uniquerowind],"$(j)/$(barsizes[j])")
			else
				labels[uniquerowind] = string(labels[uniquerowind],"$(j)")
			end
		else
			if showsize
				labels[uniquerowind] = string(labels[uniquerowind],", $(j)/$(barsizes[j])")
			else
				labels[uniquerowind] = string(labels[uniquerowind],", $(j)")
			end
		end
	end
	if any(rows[:,2].!=Inf)
		topheight = 1.1*maximum(rows[:,2][rows[:,2].!=Inf])
	else
		topheight = 0
	end
	infrows = find(rows[:,2].==Inf)
	finrows = find(rows[:,2].!=Inf)

	x0 = rows[infrows,1]
	y0 = x0
	l0 = labels[infrows]
	x1 = rows[finrows,1]
	y1 = rows[finrows,2]
	l1 = labels[finrows]
	x2 = [minimum(rows[:,1]),maximum(rows[:,1])]
	y2 = [topheight,topheight]
	return x0,y0,l0,x1,y1,l1,x2,y2
end

function plotpersistencediagram_pjs(C;dim = 1,showlabels = false,ocf = false)
	if showlabels == true
		mmode = "markers+text"
	else
		mmode = "markers"
	end

	x0,y0,l0,x1,y1,l1,x2,y2 = getpersistencediagramprimitives(
		C;
		dim = dim,
		ocf = ocf)

	if length(x0)==0 && length(x1) == 0
		print("There are no persistent homology classes in dimension $(dim).")
		return
	end

    trace0 = PlotlyJS.scatter(
    	x=x0,
    	y=y0,
        text=l0,
        mode=mmode,
        marker_color = "red",
        textposition="top center",
        hoverinfo = "x+text",
        marker_size=6,
        textfont_family="Raleway, sans-serif")
    trace1 = PlotlyJS.scatter(
    	x=x1,
    	y=y1,
        text=l1,
        mode=mmode,
        marker_color = "black",
        textposition="top center",
        hoverinfo = "x+y+text",
        marker_size=5,
        textfont_family="Raleway, sans-serif")

    L = makelayout_pjs(trace1,showlegend = false)
    PlotlyJS.plot([trace1,trace0],L)
end

function classrep_pjs(
	D::Dict;
	dim = 1,
	class=1,
	showcloud= [],
	coords = [],
	model= D["input"]["model"],
	embeddingdim = 3,
	embeddingobj = [],
	specrange = [-Inf,Inf],
	classcolor = "spectral",
	cloudcolor = [],
	textlabels = [],
	showlabels = false,
	alwaysshowcyclelabels = false)

	###
	if !(embeddingdim in [2,3])
		println("Invalid key value: <embeddingdim> must be either 2 or 3.")
		return
	end
	if embeddingobj == "hop" && showcloud == true
		println("Error: keyword <embeddingobj> may only take value <\"hop\"> when keyword <showcloud> takes value <false>.")
		println()
		return
	end

	###
	sd 	= dim+2
	facecard 	= dim+1

	##
	if showcloud == []
		showcloud = true
	end
	if embeddingobj == []
		embeddingobj = "dmat"
	end

	###
	cyclename = barname2cyclename(D,class;dim = dim)
	if haskey(D,"cyclerep")
		rep = D["cyclerep"][sd][class]
	else
		rep = getcycle(D,sd,cyclename)
	end
	classvinnewspace = incidentverts(D,sd-1,rep)
	classvinoldspace = D["nvl2ovl"][classvinnewspace]
	if showcloud == true
		vsupp = trues(length(D["farfaces"][1]))
		vsupp[classvinnewspace] = false
		compvinnewspace = find(vsupp)
		compvinoldspace = D["nvl2ovl"][compvinnewspace]
	else
		compvinoldspace = []
	end

	###
	if haskey(D["input"],"pointlabels")
		textlabels = D["input"]["pointlabels"]
	else
		m = length(D["farfaces"][1])
		textlabels = Array{String,1}(m)
		for i = 1:m
			textlabels[i] = "$(i)"
		end
	end
	if showlabels == true || showlabels == "cycle" || showlabels == "all"
		showlabtemp = true
	else
		showlabtemp = false
	end

	###
	if (classcolor == "spectral") | (embeddingobj == "hop")
		vrealization = vertexrealization(D,dim=dim,class=class)
 		vrealization = D["nvl2ovl"][vrealization]
		vertexinverter = Array{Int64}(maximum(classvinoldspace))
		vertexinverter[classvinoldspace]=1:length(classvinoldspace)
		classedges = d1faces(vrealization)
		edges_orderverts = vertexinverter[classedges]
		L = graphlaplacian_normal(edges_orderverts)
		efact_class = eigfact(L,1:4)
	end

	if showcloud && embeddingobj == "hop"
		hopedges = find(hoprange[1]<=D["grain"][2] & D["grain"][2].<=hoprange[2])
		cloudedges = vertexrealization(D,dim=1,hopedges)
		cloudedges_orderverts = vetexinverter[cloudedges]
	end

	###
	if coords == []
		if D["input"]["pc"] == "n/a"
			print("No point cloud is available.  Please consider using the mds keyword argument to generate a Euclidean embedding from the distance matrix (see documentation).")
			return "nopointcloud","nopointcloud"
		elseif D["input"]["model"]=="pc"
			coords = D["input"]["genera"]
		elseif typeof(D["input"]["pc"]) <: Array{Float64} || typeof(D["input"]["pc"]) <: Array{Int64}
			coords = D["input"]["genera"]
		end
		if size(coords,1) > 3
			print("The input coordinates have dimension greater than 3.  The generated plot will use the first three to represent each point. For more options re: high-dimensional representations, please see the documentation for multidimensional scaling.")
			coords = coords[1:3,:]
		end
		if !showcloud
			coords = coords[:,classvinoldspace]
		end
	elseif coords == "mds"
		if showcloud
			if embeddingobj == "dmat"
				metricmatrix = D["input"]["dmat"]
				metricmatrix = metricmatrix - minimum(metricmatrix)
				for i = 1:size(metricmatrix,1)
					metricmatrix[i,i]=0
				end
			elseif embeddingobj == "hop"
				metricmatrix = hopdistance(cloudedges_orderverts,inputis = "edges")
			end
		else
			if embeddingobj == "dmat"
				metricmatrix = D["input"]["dmat"][classvinoldspace,classvinoldspace]
				metricmatrix = metricmatrix - minimum(metricmatrix)
				for i = 1:size(metricmatrix,1)
					metricmatrix[i,i]=0
				end
			elseif embeddingobj == "hop"
				metricmatrix = hopdistance(edges_orderverts,inputis = "edges")
			end
		end
		coords = classical_mds(metricmatrix,embeddingdim)
		coords = round.(coords,10)
		model = "pc"
	end

	###
	if !showcloud
		textlabels = textlabels[classvinoldspace]
		classvinoldspace = 1:length(classvinoldspace)
	end

	###
	if classcolor == "spectral"
		classcolor1 = efact_class[:vectors][:,2]
		classcolor1 = (classcolor1-minimum(classcolor1))/(maximum(classcolor1)-minimum(classcolor1))
		classcolor2 = efact_class[:vectors][:,3]
		classcolor2 = (classcolor2-minimum(classcolor2))/(maximum(classcolor2)-minimum(classcolor2))
	end

	###
	T1 = maketrace_pjs(
		coords;
		model = model,
		subset = classvinoldspace,
		textlabels = textlabels,
		showlabels = showlabtemp)
	T1["marker_cmin"] = 0
	T1["marker_cmax"] = 1
	T1["marker_color"] = classcolor1
	T1["marker_line_color"] = classcolor2
	T1["marker_line_width"] = 2
	T1["marker_colorscale"] = "Jet"
	T1["marker_line_colorscale"] = "Jet"
	data = [T1]

	if !isempty(compvinoldspace)
		if showlabels == "all"
			showlabtemp = true
		else
			showlabtemp = false
		end
		T2 = maketrace_pjs(
			coords;
			model = model,
			subset = compvinoldspace,
			color = "rgb(31,119,180)",
			opacity = 0.5,
			textlabels = textlabels,
			showlabels = showlabtemp)
		append!(data,[T2])
	end

	if model == "pc"
		dim = size(coords,1)
	else
		dim = size(coords,2)
	end
	xspan = maximum(T1["x"])-minimum(T1["x"])
	yspan = maximum(T1["y"])-minimum(T1["y"])
	xa = xspan/xspan
	ya = yspan/xspan
	if dim == 3
		zspan = maximum(T1["z"])-minimum(T1["z"])
		za = zspan/xspan
	end

	layout = makelayout_pjs(data)

	return data, layout

end

function plotclassrep_pjs(
	D::Dict;
	dim = 1,
	class=1,
	showcloud = [],
	coords = [],
	model= "pc",
	embeddingdim = 3,
	embeddingobj = [],
	specrange = [-Inf,Inf],
	classcolor = "spectral",
	cloudcolor = [],
	textlabels = [],
	showlabels = "cycle")

	if D["input"]["model"]!="pc" && D["input"]["pc"] == "n/a" && coords == []
		print("No point cloud is available.  Coordinates may be supplied by the user with the <coords> keyword argument, or generated automatically via the mds keyword argument.  Please see documentation.")
		return	"nopointcloud","nopointcloud"
	end

	data,layout = classrep_pjs(
		D;
		dim = dim,
		class=class,
		showcloud = showcloud,
		coords = coords,
		model="pc",
		embeddingdim = embeddingdim,
		embeddingobj = embeddingobj,
		specrange = specrange,
		classcolor = classcolor,
		cloudcolor = cloudcolor,
		textlabels = textlabels,
		showlabels = showlabels)

	if data != "nopointcloud"
		return PlotlyJS.plot(data,layout)
	else
		return
	end
end

function betticurve(D::Dict;dim = 1,ocf = false)
	sd = dim+2
	return getbetticurve(D,sd,ocf = ocf)
end

function plotbetticurve_pjs(D::Dict;dim=1,ocf = false)
	bcu = betticurve(D;dim = dim, ocf = ocf)
	T = PlotlyJS.scatter(x = bcu[:,1],y=bcu[:,2],mode = "line",line_width=1)
	L = makelayout_pjs(T)
	PlotlyJS.plot(T,L)
end

##########################################################################################

####	BARCODE PLOTTING
####	Contributed by Paul Breiding and Sara Kališnik Verovšek, 25 January 2018

##########################################################################################

"""
	plotbarcode_pjs(C; <keyword arguments>)

Generate barcode diagram in PlotlyJS for a dictionary object `C` returned by
the function `eirene`.

# Keyword Arguments
- `dim = 0:C["input"]["maxdim"]`: (homological) homological dimensions to plot
- `sortby = "birth"`: determines whether bars appear in order of `"birth"` (that is, birth time) or `"age"` (that is, death time - birth time)
- `minage = zeros(Float64,length(dim))`: array or range of real numbers (tolerances), one for each dimension to be plotted; bars that do not meet this minimum lifetime requirement will not be plotted
- `lw = 2`: line width of bars to be plotted
"""
function plotbarcode_pjs(C::Dict; dim = 0:C["input"]["maxdim"], sortby = "birth", minage = zeros(Float64,length(dim)), lw=2)

        @assert length(dim) == length(minage) "Number of dimensions (was $(length(dim))) must be the same as the number of minageerance values (was $(length(minage)))."

        # colors = Colors.distinguishable_colors(length(dim)+1,[RGB(1,1,1)])[2:end]
        cols = Colors.colormap("Blues", mid = 0.5)
        range = Int.(round.(collect(linspace(50,100,length(dim)+1))))
        colors = map(r -> cols[r], range)

        B = [barcode(C, dim = d) for d in dim]
        #upper_limit = maximum([maximum(b[b.< Inf]) for b in B])
		upper_limit = [empteval(maximum,b[b.< Inf],0) for b in B]
		upper_limit = empteval(maximum,upper_limit,0)

        traces = map(1:length(B)) do j
             b = B[j]
             lengths = b[:,2]-b[:,1]
             b = b[lengths .> minage[j],:]
             if size(b,1) == 0
                 return [PlotlyJS.scatter(;x=[0,0], y=[0,0], mode="lines",  line_width = lw, line_color = colors[j], name = "dimension $(dim[j])")]
             end
             s = sortperm(b[:,2]-b[:,1],alg=MergeSort)
             b = b[s,:]

             i = find(x->x==Inf, b[:,2])
             b[i,2] .= 2 * upper_limit

             if sortby == "age"
             elseif sortby == "birth"
                s = sortperm(b[:,1],alg=MergeSort)
                b = b[s,:]
             else
                println("The second argument must be either \"length\" or \"lowerlimit\".")
                return 0
             end

             return [PlotlyJS.scatter(;x=b[1,:], y=[1,1], mode="lines",  line_width = lw, line_color = colors[j], name = "Dimension $(dim[j])");
             [PlotlyJS.scatter(;x=b[i,:], y=[i,i], mode="lines",  line_width = lw, line_color = colors[j], showlegend=false) for i in 2:size(b,1)]]
         end

         for i = 2:length(traces)
             for j = 1:length(traces[i])
                 traces[i][j][:y] = traces[i][j][:y] + traces[i-1][end][:y] + 10
             end
         end
         traces = vcat(traces...)

         x = maximum(vcat([t[:x] for t in traces]...))
         y = maximum(vcat([t[:y] for t in traces]...))

         layout = PlotlyJS.Layout(;
         xaxis = attr(range = [-.001, x+0.001], showgrid=false, zeroline =false, title = "ϵ"),
         yaxis = attr(range = [0,y+0.1], showgrid=false, ticks = false))
         return PlotlyJS.plot(traces[end:-1:1], layout)
end

##########################################################################################

####	SPECTRAL FUNCTIONS

##########################################################################################

function submatrixsublevellaplacianeivenstats(A;indices=1:size(A,1),threshold = Inf,statrange=1:2)
	m = length(indices)
	L = A[indices,indices]
	for alpha in L
		if alpha < 0
			print("Please note, it appears the Laplacian is being generated from a matrix with one or more negative entries.")
		end
	end
	for j = 1:m
		for i = 1:m
			if L[i,j] > threshold
				L[i,j] = 0
			end
		end
	end
	for i = 1:m
		L[i,i]=0
	end
	c = sum(L,1)
	for i = 1:m
		if c[i] != 0
			c[i] = c[i]^(-1/2)
		end
	end
	for i = 1:m
		L[i,i] = -c[i]
	end
	L = -L
	L = broadcast(*,L,c)
	L = broadcast(*,L,c')
	L = Symmetric(L)

	F = eigfact(L,statrange)
	return F
end

function graphlaplacian_normal(edges)
	numverts = maximum(edges)
	numedges = size(edges,2)

	d = zeros(Int64,numverts)
	for k in edges
		d[k]+=1
	end
	d = d.^(-1/2)
	L = zeros(numverts,numverts)
	for k = 1:numedges
		i = edges[1,k]
		j = edges[2,k]
		coeff = -d[i]*d[j]
		L[i,j]=coeff
		L[j,i]=coeff
	end
	for k = 1:numverts
		L[k,k]=1
	end
	return Symmetric(L)
end

function pcloudevec(pcloud;indices=1:size(pcloud,2),threshold = Inf,eval = 2)
	pcloud = pcloud[:,indices]

	l = length(indices)
	if isodd(l)
		A = zeros(l+1,l+1)
		A[1:l,1:l] = Distances.pairwise(Euclidean(),pcloud)
		eval = 3
	else
		A = Distances.pairwise(Euclidean(),pcloud)
	end

	F = submatrixsublevellaplacianeivenstats(A,threshold = threshold,statrange = eval:eval)
	v = F[:vectors]
	if isodd(l)
		v = v[1:end-1]
	end
	v = v[:]
	v = v-minimum(v)
	v = v/maximum(v)
	return v
end

##########################################################################################

####	GRAPHING

##########################################################################################

function maketrace_pjs(
	pcloud;
	model = "pc",
	subset = 1:size(pcloud,2),
	color = "rgb[250,250,250]",
	colorscale = "Jet",
	threshold = Inf,
	opacity = 1,
	outlineonly = false,
	textlabels = [],
	markersize = 5,
	showlabels = false)

	if model == "pc"
		dim = size(pcloud,1)
		x = pcloud[1,subset]
		y = pcloud[2,subset]
		if dim >= 3
			z = pcloud[3,subset]
			if dim > 3
				print("It appears the dimension of the input point cloud exceeds 3. Using the first three coordinates only.\n")
			end
		end
	elseif model == "points"
		dim = size(pcloud,2)
		x = pcloud[subset,1]
		y = pcloud[subset,2]
		if dim >= 3
			z = pcloud[subset,3]
			if dim > 3
				print("It appears the dimension of the input point cloud exceeds 3. Using the first three coordinates only.\n")
			end
		end
	end
	x = x[:]
	y = y[:]
	if dim == 3
		z = z[:]
	end

	if outlineonly
		symb = "circle-open"
	else
		symb = "circle"
	end

	if dim == 2
		T = PlotlyJS.scatter(x = x,y = y)
	else
		T = PlotlyJS.scatter3d(;x = x, y= y, z=z,marker_line_width=2)
	end

	T["marker_size"] = markersize
	T["marker_opacity"] = opacity
	T["marker_symbol"] = symb
	T["autocolorscale"] = false

	if !isempty(textlabels)
		T["text"] = textlabels[subset]
	elseif showlabels
		T["text"] = Array{String,1}(length(x))
		for i = 1:length(x)
			T["text"][i] = "$(i)"
		end
	end
	if showlabels
		T["mode"] = "markers+text"
	else
		T["mode"] = "markers"
	end

	if isa(color,String) && color != "spectral"
		T["marker_color"] = color
		T["marker_line_color"] = color
	else
		if color == "spectral"
			if model == "points"
				pcloud = pcloud'
			end
			v = pcloudevec(pcloud;indices=subset,threshold = threshold,eval = 2)
		elseif isa(color,Array)
			v = color
		end
		v = v - minimum(v)
		v = v/maximum(v)
		T["marker_cmin"] = 0
		T["marker_cmax"] = 1
		T["marker_color"] = v
		T["marker_colorscale"] = colorscale
	 	T["marker_line_color"] = "green"
	end
	return T
end

function makelayout_pjs(data;showlegend = false,equalaxis = true,scenecolor = "rgb(2,44,82)")
	if !(typeof(data)<:Array)
		T = data
		data = Array{Any}(1)
		data[1] = T
	end

	if data[1]["model"] == "scatter3d"
		dim=3
	else
		dim=2
	end

	if equalaxis
		xmax = -Inf
		xmin = Inf
		ymax = -Inf
		ymin = Inf
		if dim == 3
			zmax = -Inf
			zmin = Inf
		end

		for i = 1:length(data)
			T = data[i]
			if !isempty(T["x"])
				xmax = max(xmax,maximum(T["x"]))
				xmin = min(xmin,minimum(T["x"]))
				ymax = max(ymax,maximum(T["y"]))
				ymin = min(ymin,minimum(T["y"]))
				if dim == 3
					zmax = max(zmax,maximum(T["z"]))
					zmin = min(zmin,minimum(T["z"]))
				end
			end
		end
		if xmin == Inf
			L = PlotlyJS.Layout()
			return L
		end
		xspan = xmax - xmin
		yspan = ymax - ymin
		xa = xspan/xspan
		ya = yspan/xspan
		if dim == 3
			zspan = zmax-zmin
			za = zspan/xspan
		end
	else
		xa=1
		ya=1
		if dim==3
			za=1
		end
	end

	if dim==3
		L = PlotlyJS.Layout(
			showlegend = showlegend,
			width=1000,
			height=800,
			margin=attr(l=50, r=50, b=50, t=50),
			scene = attr(
				aspectratio=attr(x=xa,y=ya,z=za),
				aspectmode = "manual"
			),
		)
	else
		L = PlotlyJS.Layout(
			showlegend = showlegend,
			hovermode = "closest",
			width=1000,
			height=800,
			margin=attr(l=50, r=50, b=50, t=50),
			scene = attr(
				aspectratio=attr(
					x=xa,y=ya
				),
				aspectmode = "manual"
			),
		)
	end

	return L
end

function ezplot_pjs(
	pcloud;
	model = "pc",
	subset = 1:size(pcloud,2),
	color = "rgb[250,250,250]",
	colorscale = "Jet",
	threshold = Inf,
	opacity = 1,
	outlineonly = true,
	textlabels = [],
	markersize = 5,
	showlabels = false)

	T = maketrace_pjs(
		pcloud;
		model = model,
		subset = subset,
		color = color,
		colorscale = colorscale,
		threshold = threshold,
		opacity = opacity,
		outlineonly = outlineonly,
		textlabels = textlabels,
		markersize = markersize,
		showlabels = showlabels)

	L = makelayout_pjs(T)

	PlotlyJS.plot(T,L)
end

function edgetrace_pjs(coordinates,edges;model="pc")
	if model != "pc"
		coordinates = coordinates'
	end

	edgetraces = []
	if size(coordinates,1) == 2
		for i = 1:size(edges,2)
			verts = edges[:,i]
			if edgetraces == []
				trace = PlotlyJS.scatter(
					x = coordinates[1,verts],
					y = coordinates[2,verts],
					mode = "lines",
					name = "Dim 1 Faces")
				edgetraces = [trace]
			else
				trace = PlotlyJS.scatter(
					x = coordinates[1,verts],
					y = coordinates[2,verts],
					mode = "lines",
					name = "edge ($(verts[1]),$(verts[2]))",
					showlegend = false)
				append!(edgetraces,[trace])
			end
		end
	elseif size(coordinates,1) == 3
		for i = 1:size(edges,2)
			verts = edges[:,i]
			if edgetraces == []
				trace = PlotlyJS.scatter3d(
					x = coordinates[1,verts],
					y = coordinates[2,verts],
					z = coordinates[3,verts],
					mode = "lines",
					opacity = 0.5,
					name = "Dim 1 Faces")
				edgetraces = [trace]
			else
				trace = PlotlyJS.scatter3d(
					x = coordinates[1,verts],
					y = coordinates[2,verts],
					z = coordinates[3,verts],
					name = "edge ($(verts[1]),$(verts[2]))",
					showlegend = false,
					mode = "lines",
					opacity = 0.5)
				append!(edgetraces,[trace])
			end
		end
	else
		print("It appears the coordinates provided have dimension other than 2 or 3; these are the only two currently supported.")
	end
	return edgetraces
end

function d1faces(facesbycol)
	facecard = size(facesbycol,1)
	numfaces = size(facesbycol,2)
	M = maximum(facesbycol)
	supp = falses(M)
	for m in facesbycol
		supp[m] = true
	end
	vertices = find(supp)
	numverts = length(vertices)
	translator = Array{Int64}(M)
	translator[vertices]=1:numverts
	supp = falses(numverts,numverts)
	for i = 1:facecard-1
		for j = (i+1):facecard
			for k = 1:numfaces
				row = translator[facesbycol[i,k]]
				col = translator[facesbycol[j,k]]
				supp[max(row,col),min(row,col)]=true
			end
		end
	end
	numedges = countnz(supp)
	edges = Array{Int64}(2,numedges)
	counter = 1
	for j = 1:numverts
		for i = (j+1):numverts
			if supp[i,j]
				edges[1,counter]=j
				edges[2,counter]=i
				counter+=1
			end
		end
	end
	return vertices[edges]
end

##########################################################################################

####	MISC

##########################################################################################

function pairwiseisequal(X;under=identity)
	l 					= 	length(X)
	for p 				= 	1:l
		if !isequal(under(copy(X[1])),under(copy(X[p])))
			return 			false
		end
	end
	return 					true
end

function csvimport2linends(M)
	m,n = size(M)
	endpoints = zeros(Int64,m)
	for p = 1:m
		for q = 1:n
			if M[p,q] == ""
				endpoints[p] = q-1
				break
			elseif q == n
				endpoints[p] = n
			end
		end
	end
	return endpoints
end

# stands for extension-by-constant
function ec(v,p,k)
	if 0 < p <= length(v)
		return v[p]
	else
		return k
	end
end

# stands for number of simplicies of cardinality less than or equal to k
# A is an array of arrays
# k is an integer
# this function is not used at the time of this writing (jan 14, 2018)
function numsimcardlek(A,k)
	c = 0;
	for p = 1:k
		c 	+= length(A[p])
	end
	return c
end

function undercat(X)
	l = length(X);
	if l <= 1
		return X
	else
		m = l+l-1;
		Y = Array{Any}(m)
		Y[2:2:m] = "_"
		Y[1:2:m] = X
		return string(Y...)
	end
end

function ezread(s)
	if s[end-2:end] == "csv"
		return readdlm(s,',','\r')
	elseif s[end-2:end] == "txt"
		return readdlm(s)
	else
		println("Please ensure the input file is either comma

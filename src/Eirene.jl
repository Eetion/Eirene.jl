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

#### 	SIMPLICIAL CONSTRUCTIONS

##########################################################################################

# NOTA BENE: eirene permutes vertex labels prior to calculation;
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

# NOTA BENE: eirene permutes vertex labels prior to calculation;
# all versions of <vertexrealization> must be interpreted with
# respect to the PERMUTED labeling scheme
function vertexrealization(D::Dict,facecardinality,facenames)
	return vertexrealization(D["farfaces"],D["firstv"],facecardinality,facenames)
end

# NOTA BENE: eirene permutes vertex labels prior to calculation;
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

# NOTA BENE: eirene permutes vertex labels prior to calculation;
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

# NOTA BENE: eirene permutes vertex labels prior to calculation;
# the output of <incidentverts> must be interpreted with respect
# to the PERMUTED labeling scheme
function incidentverts(D::Dict,facecardinality,facenames)
	return incidentverts(D["farfaces"],D["firstv"],facecardinality,facenames)
end

# NOTA BENE: eirene permutes vertex labels prior to calculation;
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

function ff2complex(farfaces,firstv;maxcard = length(farfaces))
	Nrv 	= fill(Array{Int64}(0),maxcard)
	Ncp 	= fill(Array{Int64}(0),maxcard)
	Nrv		= convert(Array{Array{Int64,1}},Nrv)
	Ncp		= convert(Array{Array{Int64,1}},Ncp)
	Nrv[1] 	= Array{Int64}(0)
	Ncp[1]	= fill(1,length(farfaces[1])+1)
	for sd = 2:maxcard
		Nrv[sd],Ncp[sd] = ff2boundary(farfaces,firstv,sd=sd)
	end
	return Nrv,Ncp
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

NOTA BENE
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

	if !isempty(nplows) && sd > 2
		npfilt = grain[sd-1][nplows]
		nporder = integersinsameorder(npfilt)
		addinteger!(nporder,numppair)
	else
		nporder = (numppair+1):(numppair+length(nplows))
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

function grain2maxcard(grain)
	c = 0
	for i = 1:length(grain)
		if !isempty(grain[i])
			c = i
		end
	end
	return c
end

function skelcount(numvertices,maxcardinality)
	c = 0
	for i = 1:maxcardinality
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
-	the columns of M must be ordered by grain
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

function persistF2_core_cell(
	Nrv,
	Ncp,
	grain;
	maxcard = length(Nrv),
	record="cyclerep",
	verbose=false,
	prepairs = fill(Array{Int64}(0),maxcard+1)
	)
	if record == "all" || record == "cyclerep"
		storetransform = true
	else
		storetransform = false
	end
	tcp::Array{Array{Int64,1},1}			=Array{Array{Int64,1},1}(maxcard+1)
	trv::Array{Array{Int64,1},1}			=Array{Array{Int64,1},1}(maxcard+1)
	phi::Array{Array{Int64,1},1}			=Array{Array{Int64,1},1}(maxcard+1)
	plo::Array{Array{Int64,1},1}			=Array{Array{Int64,1},1}(maxcard+1)
	tid::Array{Array{Int64,1},1}			=Array{Array{Int64,1},1}(maxcard+1)
	for i in [1,maxcard+1]
		tcp[i] 				=[1]
		trv[i]				=Array{Int64}(0)
		phi[i]				=Array{Int64}(0)
		plo[i]				=Array{Int64}(0)
		tid[i]				=Array{Int64}(0)
	end
	maxnzs 									=zeros(Int64,maxcard+1)
	m 										=length(Ncp[1])-1
	for sd = 2:maxcard
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
			trv[sd] 	= Array{Int64}(0)
			tcp[sd] 	= [1]
			tid[sd] 	= Array{Int64}(0)
			plo[sd] 	= Array{Int64}(0)
			phi[sd] 	= Array{Int64}(0)
			continue
		end
		Mm0				= maximum(Mrv) # only used for the definition of the low-translator, which I believe is only there for debugging purposes
		Mm				= [maximum(Mrv)]  # temporary
		Mn 				= [length(Mcp)-1] # temporary
		higlab			= convert(Array{Int64},1:Mn[1])	# we'll assume prepairs to be empty, for now
		lowlab	 		= intervalcomplementuniqueunsortedinput(lowbasisnames,Mm[1])	# temporary, and we'll assume prepairs is empty for now
		nporder			= sortperm(grain[sd-1][lowlab])
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

function persistF2_core_vr(
	farfaces::Array{Array{Int64,1},1},
	firstv::Array{Array{Int64,1},1},
	prepairs::Array{Array{Int64,1},1},
	grain::Array{Array{Int64,1},1},
	maxcard::Integer;
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

	trv::Array{Array{Int64,1},1}			=Array{Array{Int64,1},1}(maxcard+1);
	tcp::Array{Array{Int64,1},1}			=Array{Array{Int64,1},1}(maxcard+1);
	phi::Array{Array{Int64,1},1}			=Array{Array{Int64,1},1}(maxcard+1);
	plo::Array{Array{Int64,1},1}			=Array{Array{Int64,1},1}(maxcard+1);
	tid::Array{Array{Int64,1},1}			=Array{Array{Int64,1},1}(maxcard+1);
	for i in [1,maxcard+1]
		tcp[i]  			=[1];
		trv[i]				=Array{Int64}(0)
		tid[i]				=Array{Int64}(0)
		phi[i]				=Array{Int64}(0)
		plo[i]				=Array{Int64}(0)
	end

	maxnzs 									=zeros(Int64,maxcard+1);

	for sd = 2:maxcard
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

function persistF2!(
	D::Dict;maxcard=0,
	dictionaryoutput::Bool = true,
	verbose::Bool = false,
	record = "cyclerep")

	farfaces = D["farfaces"]
	firstv = D["firstv"]
	prepairs = D["prepairs"]
	grain = D["grain"]
	if maxcard == 0
		maxcard = length(farfaces)-1
	end

	trv,tcp,plo,phi,tid,maxnzs =
	persistF2_core_vr(farfaces,firstv,prepairs,grain,maxcard::Integer;record=record,verbose = verbose)
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

function persistF2_vr(
	s,
	maxcard;
	model 			= "vr",
	filetype 		= "textfile",
	minrad			= -Inf,
	maxrad			= Inf,
	numrad			= Inf,
	fastop			= true,
	vscale			= "default",
	pointlabels 	= [],
	verbose 		= false,
	record 			= "cyclerep")

	#### Start timer
# 	tic()

	#### Extract data as necessary
	inputisfile = false
	if typeof(s) == String
		inputisfile = true
		filename = s #modified 12/29/2017; note that strings are immutable
		if filetype == "textfile"
			if typeof(readdlm(filename,','))<:Array{Float64}
				s = readdlm(s,',')
			elseif typeof(readdlm(filename,' '))<:Array{Float64}
				s = readdlm(s,' ')
			else
				print("Error reading text file.  Input files of this format must be either comma or space delimited.")
				return
			end
		elseif filetype == "perseus_distmat"
			s,minrad,maxrad,maxcard = parse_perseusdistmat(filename)
			model = "clique_perseus"
			numrad = Inf
		elseif filetype == "perseus_brips"
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
		"maxdim" 		=> maxcard-2,
		"maxrad"		=> maxrad,
		"minrad"		=> minrad,
		"numrad"		=> numrad,
		"fastop"		=> fastop,
		"record"	 	=> record,
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
	elseif model == "vr"
		d = convert(Array{Float64,2},s)
	end

	(t,ocg2rad) = ordercanonicalform_3(
		d;
		minrad=minrad,
		maxrad=maxrad,
		numrad=numrad,
		fastop=fastop,
		vscale=vscale,
		verbose = verbose)

	#### Build the complex
	D = buildcomplex3(t,maxcard;verbose = verbose)
	D["ocg2rad"]=ocg2rad

	#### Compute persistence
	persistF2!(D;verbose = verbose,record = record)

	#### Store input data
	D["input"] = input

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
	- line 4: coundary matric column pattern
=#
function persistF2_cell(filepath::String;
						maxdim=Inf,
						record = "cyclerep",
						verbose=false)
	rv,cp,fv,dp 	= 	filepath2unsegmentedfilteredcomplex(filepath)
	C  				= 	persistF2_cell(
						rv,cp,fv,dp,
						maxdim=maxdim,
						record=record,
						verbose=verbose)
	C["input"]["source"]	= 	filepath;
	return C
end

function filepath2unsegmentedfilteredcomplex(filepath)
	M = readcsv(filepath)
	if size(M,1) != 4
		if M == ones(Float64,2,1)
			rv 	= zeros(Int64,0)
			cp 	= ones(Int64,1)
			fv	= zeros(Float64,0)
			dp  = ones(Int64,1)
			return rv,cp,fv,dp
		else
			print("Error: input file should have .csv format with four lines.")
			return
		end
	end
	n = size(M,2)
	endpoints = zeros(Int64,4)
	for p = 1:4
		for q = 1:n
			if M[p,q] == ""
				endpoints[p] = q-1
				break
			elseif q == n
				endpoints[p] = n
			end
		end
	end
	dp 	    = M[1,1:endpoints[1]]  # stands for dimension pattern
	fv      = M[2,1:endpoints[2]]  # stands for filtration values
	rv 	    = M[3,1:endpoints[3]]  # stands for row values
	cp      = M[4,1:endpoints[4]]  # stands for column pattern
	dp 		= convert(Array{Int64,1},dp)
	fv      = convert(Array{Float64,1},fv)
	rv 		= convert(Array{Int64,1},rv)
	cp 		= convert(Array{Int64,1},cp)
	return rv,cp,fv,dp
end

function unsegmentedfilteredcomplex2segmentedfilteredcomplex(rv,cp,fv,dp;ndim=Inf)
	# ndim stands for number of dimensions
	# nsd stands for number of stored dimensions
	nsd 	= 	length(dp)-1
	if ndim == Inf
		ndim  	= nsd
	end
	m 		= 	min(nsd,ndim)

	fvc     = Array{Array{Float64,1}}(ndim+1)
	rvc     = Array{Array{Int64,1}}(ndim+1)
	cpc     = Array{Array{Int64,1}}(ndim+1)

	# remarks:
	# (a) #{cells of dimension ≤ (p   = k+1)} 				= dp[k+2]-1		= 	dp[p+1]-1
	# (b) #{cells of dimension ≤ (p-2 = (k-2)+1 = k-1)} 	= dp[k]-1 		=	dp[p-1]-1
	for p = 1:m
		rvc[p],cpc[p]   = 	copycolumnsubmatrix(rv,cp,cran(dp,p))
		rvc[p]			= 	rvc[p] - ec(dp,p-1,0) + 1  # we subtract off the number of cells of dimension 2 less than those of interest, since starts the _faces_ of the cells of interest off at index 1; we get away with this in dimension 0 bc rvc[1] is the empty set
		fvc[p] 			=   convert(Array{Float64,1},fv[cran(dp,p)])
	end

	for p = (m+1):(ndim+1)
		rvc[p]		=	zeros(Int64,0)
		cpc[p]		=	ones(Int64,1)
		fvc[p] 	=	zeros(Int64,0)
	end
	return rvc,cpc,fvc
end

# version with 4 (integer array) non-keyword arguments
# rv:row values, cp: column pattern, fv: filtration values, dp: dimension pattern
function persistF2_cell(rv::Array{Int64,1},
						cp::Array{Int64,1},
						fv::Array{Float64,1},
						dp::Array{Int64,1};
						maxdim=length(dp),
						record = "cyclerep",
						verbose=false)
	rvc,cpc,fvc = unsegmentedfilteredcomplex2segmentedfilteredcomplex(rv,cp,fv,dp;ndim=maxdim+2)
	return persistF2_cell(rvc,cpc,fvc;maxdim=maxdim,record = record,verbose=false)
end

# version with 3 (array of arrays) non-keyword arguments
function persistF2_cell(rv,cp,filt;maxdim=length(rv),record = "cyclerep",verbose=false)

	### Record the input parameters
	input = Dict(
		"model"			=> "complex",
		"version" 		=> Pkg.installed("Eirene"),
		"date"			=> string(Dates.Date(now())),
		"time"			=> string(Dates.Time(now())),
		"maxdim"		=> maxdim,
		"record" 		=> record
		)

	n = length(rv)

	### Format the grain data
	# Concatenate the grain vectors
	numcells = 0
	for i = 1:n
		numcells+=length(filt[i])
	end
	filt1 = Array{Float64}(numcells)
	numcells = 0
	for i = 1:n
		l = length(filt[i])
		filt1[numcells+1:numcells+l] = filt[i]
		numcells += l
	end

	# Convert to order canonical form
	filt2 = integersinoppositeorder_nonunique(filt1)
	ocff = fill(Array{Int64}(0),n+3)
	numcells = 0
	for i = 1:n
		l = length(filt[i])
		ocff[i] = filt2[numcells+1:numcells+l]
		numcells += l
	end

	# Compute the grain translator
 	ocg2rad = sort(unique(filt1))
 	ocg2rad = flipdim(ocg2rad,1)

	### Perform the persistence computation
	trv,tcp,plo,phi,tid =
	persistF2_core_cell(
		rv,
		cp,
		ocff;
		maxcard = maxdim+2,
		record=record,
		verbose=false,
		prepairs = convert(Array{Array{Int64,1},1},fill(Array{Int64}(0),maxdim+2))
		)

	### Create the dictionary that stores all relevant data
	D = Dict(
		"rv" 		=> rv,
		"cp" 		=> cp,
		"grain"		=> ocff,
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
	maxcard = convert(Int64,A[2,4])+1
	return s,minrad,maxrad,maxcard
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
	numnlpl 	= numl-length(plo[sd-1])
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
		Lrv,Lcp,Lirv,Licp,Rrv,Rcp = boundarylimit_cell(rv,cp,trv,tcp,plo,phi,tid,sd;verbose=false)
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

	lowtranslator = zeros(Int64,numl)
	lowtranslator[tid[sd]] = 1:numnlpl
	dummy0 = zeros(Int64,0)

	brv,dummy1,bcp,dummy2 = stackedsubmatrices(rv[sd],cp[sd],tid[sd],dummy0,phi[sd],numl)
	yafterx!(lowtranslator,brv)
	append!(bcp,fill(bcp[end],numnlpl-nump)) # <-- note there's no +1 b/c bcp is already 1 elt. longer than the # of columns
	return brv,bcp
end

function unpack!(D::Dict)
	l = length(D["grain"])
	maxcard = D["input"]["maxdim"]+2 # grain2maxcard(D["grain"])

	Lirv = Array{Array{Int64,1},1}(l);  Licp = Array{Array{Int64,1},1}(l)
	Lrv  = Array{Array{Int64,1},1}(l);  Lcp  = Array{Array{Int64,1},1}(l)
	Rrv  = Array{Array{Int64,1},1}(l);  Rcp  = Array{Array{Int64,1},1}(l)

	if 	D["input"]["record"] == "all"
		N 	= 	maxcard
	else
		N 	= 	maxcard-1
		if maxcard >= 2
			Lirv[maxcard],Licp[maxcard] = boundarylimit_Lionly(D["trv"],D["tcp"],D["tid"],maxcard)
		elseif maxcard == 1
			Lirv[maxcard] = Array{Int64,1}(0)
			Licp[maxcard] = ones(Int64,1)
		end
		Lrv[maxcard]=Array{Int64,1}(0)
		Lcp[maxcard]=Array{Int64,1}(0)
		Rrv[maxcard]=Array{Int64,1}(0)
		Rcp[maxcard]=Array{Int64,1}(0)
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

	D["cyclerep"] = fill(Array{Array{Int64,1},1}(0),maxcard)  # for each dimension one stores an array of arrays

	for i = 2:maxcard
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

function copycolumnsubmatrix{Tv<:Integer}(Arv::Array{Tv,1},Acp::Array{Tv,1},columnindices::Array{Tv,1})
	allocationspace = 0
	for j in columnindices
		allocationspace+= Acp[j+1]-Acp[j]
	end
	Brv = Array{Tv}(allocationspace)
	Bcp = Array{Tv}(length(columnindices)+1)
	Bcp[1]=1
	for jp = 1:length(columnindices)
		j = columnindices[jp]
		Bcp[jp+1] = Bcp[jp]+Acp[j+1]-Acp[j]
		Brv[Bcp[jp]:(Bcp[jp+1]-1)]=Arv[Acp[j]:(Acp[j+1]-1)]
	end
	return Brv,Bcp
end

function copycolumnsubmatrix{Tv<:Integer}(Arv::Array{Tv,1},Acp::Array{Tv,1},columnindices::UnitRange{Int64})
	allocationspace = 0
	for j in columnindices
		allocationspace+= Acp[j+1]-Acp[j]
	end
	Brv = Array{Tv}(allocationspace)
	Bcp = Array{Tv}(length(columnindices)+1)
	Bcp[1]=1
	for jp = 1:length(columnindices)
		j = columnindices[jp]
		Bcp[jp+1] = Bcp[jp]+Acp[j+1]-Acp[j]
		Brv[Bcp[jp]:(Bcp[jp+1]-1)]=Arv[Acp[j]:(Acp[j+1]-1)]
	end
	return Brv,Bcp
end

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
	x = cos(theta)
	y = sin(theta)
	z = zeros(100)
	pcloud = hcat(x,y,z)'
	pcloud = hcat(pcloud,pcloud,pcloud)
	pcloud = pcloud + 0.3*rand(3,300)
	return pcloud
end

function torus(;m = 100,n = 50,mrad=1,nrad = 0.2)
	theta = (1:m)*2*pi/m;
	torus1 = mrad*repmat(cos(theta),1,n)
	torus2 = mrad*repmat(sin(theta),1,n)
	torus3 = nrad*repmat(sin(2*pi*(1:n)./n),1,m)'
	torus4 = repmat(cos(2*pi*(1:n)./n),1,m)'
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
		numpts = round(Int64,40*cos(t))
		theta = 2*pi*(1:numpts)/numpts
		append!(x,cos(t)*cos(theta))
		append!(y,cos(t)*sin(theta))
		append!(z,fill(sin(t),numpts))
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

function chessboardcomplex_symmat(m,n)
	#### 1-A-I, where A is the adjacency matrix of the chessboard complex C(m,n)
	#### and I is identity
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
	theta1 = A[1,:]
	theta2 = A[2,:]
	x = cos(theta1)
	y = sin(theta1)
	z = 0.25*sin(theta2)
	alpha = 1-0.25*cos(theta2)
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
	x = cosd(theta2)
	y = sind(theta2)
	z = sind(theta1)
	w = cosd(theta1)
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

function simplecityfilepath()
	fp = joinpath(@__DIR__,"examples/simplemapscitydata.csv")
	return fp
end

##########################################################################################

####	MATRIX WEIGHTS AND FORMATTING

##########################################################################################

function customceil(N,a,b,numrad)
	S = copy(N)
	ran = linspace(a,b,numrad)
	ran = convert(Array{Float64},ran)
	ran[1] 		= a
	ran[end] 	= b
	append!(ran,[Inf])
	for j = 1:(m-1)
		for i = (j+1):m
			post = 1
			while ran[post] < N[i,j]
				post+=1
			end
			S[i,j] = ran[post]
			S[j,i] = ran[post]
		end
	end
	return S
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

			Q = customceil(copy(M),minentry,upperlim,numrad)
			N = customceil(copy(M),minentry,maxentry,numrad)
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
end

function offdiagmin(S,i)
	if i == 1
		return(minimum(S[2:end,i]))
	elseif i == size(S,1)
		return minimum(S[1:end-1,i])
	else
		return min(minimum(S[1:i-1,i]),minimum(S[i+1:end,i]))
	end
end

function ordercanonicalform_3{Tv}(
	S::Array{Tv};
	maxrad=Inf,
	minrad=-Inf,
	numrad=Inf,
	vscale="default",
	fastop::Bool=true,
	verbose::Bool=false)

	# Assumes that the input array S has only finite entries.
	# The value for keyword <numrad> must be either a positive integer or Inf

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
	p 						= sortperm(vec(symmat))

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
		return round(Int65,symmat),ocg2rad
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
	p = sortperm(symmat[:])
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
	p 			= sortperm(v)
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
	if length(D["farfaces"][sd])==0
		return Array{Float64}(0,2)
	end

	maxrad = Int(maximum(D["grain"][2]))
	v = zeros(Int64,maxrad)

	bco = barcode(D;dim = sd-2,ocf=true)
	bco[bco.==Inf] = maxrad
	bco = convert(Array{Int64},bco)

	for i = 1:size(bco,1)
		v[bco[i,1]:bco[i,2]]+=1
	end

	if ocf == false
		u = sort(D["ocg2rad"])
	else
		u = Array(1:maxrad)
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
	for p in cyclenumber
		rep[p] 	= getcycle_cell(rv,cp,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plo,phi,tid,sd,p)
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

function barcode(D::Dict;dim = 1,ocf = false,verbose = false, givenztidindices = false)
	sd = dim+2
	plo = D["plo"][sd]
	phi = D["phi"][sd]
	higfilt = D["grain"][sd][phi]
	lowfilt = D["grain"][sd-1][plo]
	nump = length(plo)
	numnzmortalbars = countnz(higfilt.!=lowfilt)

	tid = D["tid"][sd]
	numbars = length(tid)
	numnzevergreenbars = numbars-nump
	numnzbars = numnzmortalbars + numnzevergreenbars
	numlows = length(D["grain"][sd-1])

	translate2plowindex			= zeros(Int64,numlows)
	translate2plowindex[plo]  = 1:nump
	nzbarcounter = 0
	deathtimes = fill(Inf64,numbars)
	for i = 1:numbars
		k = translate2plowindex[tid[i]]
		if k>0
			deathtimes[i] = D["grain"][sd][phi[k]]
		end
		if deathtimes[i] != D["grain"][sd-1][tid[i]]
			nzbarcounter+=1
		end
	end

	if verbose
		println(["typeof(deathtimes)" typeof(deathtimes)])
		println(["length(tid)" length(tid)
		"length(D[plo])" length(plo)
		"numinf(deathtimes)" countnz(deathtimes.==Inf)])
	end
	summary = Array{Any}(nzbarcounter,2)
	nzbarcounter = 0
	for i = 1:numbars
		if deathtimes[i] != D["grain"][sd-1][tid[i]]
			nzbarcounter+=1
			summary[nzbarcounter,1] = D["grain"][sd-1][tid[i]]
			summary[nzbarcounter,2] = deathtimes[i]
		end
	end
	if ocf == false
		summary[:,1]=D["ocg2rad"][convert(Array{Int64},summary[:,1])]
		summary[1:numnzmortalbars,2] = D["ocg2rad"][convert(Array{Int64},summary[1:numnzmortalbars,2])]
	else
		summary = 1+maximum(D["grain"][2])-summary
		summary = convert(Array{Float64},summary)
		summary[summary.==-Inf] = Inf
	end
	if givenztidindices
		tidindices = vcat(find(higfilt.!=lowfilt),Array(length(plo)+1:length(tid)))
		return (summary,tidindices)
	else
		return summary
	end
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
		coords = round(coords,10)
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
             s = sortperm(b[:,2]-b[:,1])
             b = b[s,:]

             i = find(x->x==Inf, b[:,2])
             b[i,2] .= 2 * upper_limit

             if sortby == "age"
             elseif sortby == "birth"
                s = sortperm(b[:,1])
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

# stands for extension-by-constant
function ec(v,p,k)
	if 0 < p <= length(v)
		return v[p]
	else
		return k
	end
end

function ezread(s)
	if s[end-2:end] == "csv"
		return readdlm(s,',','\r')
	elseif s[end-2:end] == "txt"
		return readdlm(s)
	else
		println("Please ensure the input file is either comma separated (.csv)
		or space delimited (.prn)")
	end
end

function ezlabel(y)
	l = length(y)
	x = Array{String,1}(l)
	for i = 1:l
		x[i] = try
			convert(String,y[i])
		catch
			"$(i)"
		end
	end
	return x
end

function separatelabels(s,side)
	if side == "left" || side == "right"
		numlabels = size(s,1)
	else
		numlabels = size(s,2)
	end
	labels = Array{String,1}(numlabels)

	if side == "none"
		for i = 1:numlabels
			labels[i] = "$i"
		end
	end

	if side == "left"
		for i = 1:numlabels
			labels[i] = "$(s[i,1])"
		end
		s = s[:,2:end]
	elseif side == "right"
		for i = 1:numlabels
			labels[i] = "$(s[i,end])"
		end
		s = s[:,1:end-1]
	elseif side == "top"
		for i = 1:numlabels
			labels[i] = "$(s[1,i])"
		end
		s = s[2:end,:]
	elseif side == "bottom"
		for i = 1:numlabels
			labels[i] = "$(s[end,i])"
		end
		s = s[1:end-1,:]
	end
	return s,labels
end

function binom(x,y)
	k = 1;
	for i = x:-1:(x-y+1)
		k = i*k
	end
	for i = 2:y
		k = k/i
	end
	k = convert(Int64,k)
	return k
end

function binom_float(x,y)
	k = 1;
	a = x
	b = 1
	for i = 1:y
		k = k*a/b
		a-=1
		b+=1
	end
	return k
end

function yafterx!{Tv<:Integer}(y::Array{Tv,1},x::Array{Tv,1})
	for i = 1:length(x)
		x[i] = y[x[i]]
	end
end

function yafterx{Tv}(y::AbstractVector{Tv},x)
	z = Array{Tv}(length(x))
	for i = 1:length(x)
		z[i] = y[x[i]]
	end
	return z
end

function ss2full(rowval,colptr,m)
# 	renamed from 'showfull' on 12/27/2017; ss stands for 'sparse support'
	n = length(colptr)-1
	M = zeros(Int8,m,n)
	for j = 1:n
		M[rowval[cran(colptr,j)],j]=1
	end
	return M
end

function ss2full(rowval,colptr,m,n)
# 	renamed from 'showfull' on 12/27/2017; ss stands for 'sparse support'
	M = zeros(Int8,m,n)
	for j = 1:n
		M[rowval[cran(colptr,j)],j]=1
	end
	return M
end

function full2ss(A)
# 	added 12/27/2017; ss stands for 'sparse support'
# 	input: a full array A
# 	output: the support of A, encoded in sparse column format
	m,n = size(A)
	rv  = find(A)
	rv  = mod.(rv-1,m)+1
	cp  = zeros(Int64,n+1)
	cp[1] = 1
	for p = 1:n
		cp[p+1] = cp[p]+countnz(A[:,p])
	end
	return rv,cp
end

function sparsifydesparsifytest(m,n)
# 	added 12/27/2017
	for p = 1:m
		A = rand(n,n).<0.1
		A = convert(Array{Int64},A)
		rv,cp = full2ss(A)
		B = ss2full(rv,cp,n)
		if A != B
			print("error on iteration $(p)")
			return A,B
		end
	end
	print("test successful")
end

function colsupportsum(colptr,n::Integer)
	x = Array{Int64}(n)
	@inbounds begin
	for i = 1:n
		x[i] = colptr[i+1]-colptr[i]
	end
	end
	return x
end

function rowsupportsum(Arv,Acp,Am::Int64,cols)
	x = zeros(Int64,Am)
	for j in cols
		for ip in cran(Acp,j)
			x[Arv[ip]]+=1
		end
	end
	return x
end

function getcolptr2{Tv<:Integer}(orderedpositiveintegerlist::Array{Tv,1},howfartolookbeforestopping::Tv)
	#### please note: order can be ascending or descending
	v = orderedpositiveintegerlist
	if isempty(v)
		return []
	end
	colptr = Array{Int64}(length(v)+1)
	colptr[1] = 1
	transitioncounter = 1
	currentvalue = v[1]
	for i = 2:howfartolookbeforestopping
		if v[i]!=currentvalue
			transitioncounter+=1
			colptr[transitioncounter]=i
			currentvalue = v[i]
		end
	end
	colptr[transitioncounter+1]=howfartolookbeforestopping+1
	deleteat!(colptr,(transitioncounter+2):(length(v)+1))
	return colptr
end

function addinteger!{Tv}(v::Array{Tv,1},k::Int64)
	for i = 1:length(v)
		v[i]+=k
	end
end

function sparseadjacencymatrix(A;inputis = "adjacencymatrix")
	if inputis == "adjacencymatrix"
		if size(A,1) != size(A,2)
			print("Error: unless the <inputis> keywork argument has value 'edges', the input array must be square.")
			return
		end
		m = size(A,1)
		rv = Array{Int64}(0)
		cp = zeros(Int64,m+1)
		cp[1] = 1
		for i = 1:m
			adjverts = find(A[:,i])
			append!(rv,adjverts)
			cp[i+1] = cp[i]+length(adjverts)
		end
		return rv, cp
	elseif inputis == "edges"
		if size(A,1) != 2
			print("Error: when the <inputis> keywork argument has value 'edges', the input array must have exactly two rows.")
		end
		m = maximum(A)
		adjmat = falses(m,m)
		for i = 1:size(A,2)
			adjmat[A[1,i],A[2,i]]=true
		end
		adjmat[find(transpose(adjmat))] = true
		return sparseadjacencymatrix(adjmat)
	end
end

function hopdistance_sparse(rv,cp)
	m = length(cp)-1
	H = zeros(Int64,m,m)
	for i = 1:m
		c = 0
		metnodes = falses(m)
		metnodes[i] = true
		fringenodes = falses(m)
		fringenodes[i] = true
		fringelist = [i]

		while !isempty(fringelist)
			c+=1
			for j in fringelist
				for k in crows(cp,rv,j)
					if !metnodes[k]
						metnodes[k] = true
						fringenodes[k] = true
						H[k,i] = c
					end
				end
			end
			fringelist = find(fringenodes)
			fringenodes[:] = false
		end
		H[!metnodes,i]=m+1
	end
	return H
end

function hopdistance(rv,cp)
	return hopdistance_sparse(rv,cp)
end

function hopdistance(A;inputis = "fulladj")
	if inputis == "fulladj"
		rv,cp = sparseadjacencymatrix(A)
	elseif inputis == "edges"
		rv,cp = sparseadjacencymatrix(A,inputis="edges")
	else
		rv,cp = A
	end
	return hopdistance_sparse(rv,cp)
end


##########################################################################################

####	VIETORIS-RIPS CONSTRUCTION

##########################################################################################

function buildcomplex3{Tv}(symmat::Array{Tv},maxcard; dictionaryoutput = true, verbose = false)

	m = size(symmat,1)
	w = vec(mean(symmat,1))

	vperm = sortperm(-w)
	symmat = symmat[vperm,vperm]

	grain = Array{Array{Int64,1}}(maxcard+1)
	farfaces = Array{Array{Int64,1}}(maxcard+1)
	prepairs = Array{Array{Int64,1}}(maxcard+1)
	firstv = Array{Array{Int64,1}}(maxcard+1)

	farfaces[maxcard+1] = Array{Int64}(0)
	firstv[maxcard+1] = ones(Int64,1)
	grain[maxcard+1] = Array{Int64}(0)
	prepairs[maxcard+1] = Array{Int64}(0)

	farfaces[1] = convert(Array,1:m)
	firstv[1] = convert(Array,1:(m+1))
	grain[1] = diag(symmat)
	prepairs[1] = Array{Int64}(0)

	r,c,z = generate2faces(symmat)
	farfaces[2] = r
	firstv[2] = c
	grain[2] = z
	prepairs[2] = Array{Int64}(0)

	if maxcard == 3
		generate3faces!(farfaces,firstv,grain,prepairs,m,symmat;verbose = verbose)
		if dictionaryoutput == true
			D = Dict(
				"farfaces" => farfaces,
				"firstv" => firstv,
				"grain" => grain,
				"prepairs" => prepairs,
				"symmat" => symmat,
				"nvl2ovl"=>vperm)
			return D
		else
			return farfaces,firstv,grain,prepairs,symmat,vperm
		end
	end

	fpi = Array{Int64}(0)
	ff2pv = Array{Int64}(0)
	pmhist = zeros(Int64,m,m)

	for sd = 3:maxcard
		if verbose
			print(["set cardinality = " sd])
			println(["num sd-1 cells" length(farfaces[sd-1])])
		end

		nl = length(farfaces[sd-1])
		nll = length(farfaces[sd-2])

		startlength = nl
		stepsize = min(10^7,Int(ceil(nl/4)))

		npsupp = trues(nl)
		pflist = Array{Int64}(nl)
		jrv = farfaces[sd-1]
		jcp = firstv[sd-1]
		jz = grain[sd-1]
		zll= grain[sd-2]
		izfull = Array{Int}(nll)
		r = Array{Int64}(startlength)
		z = Array{Int64}(startlength)
		c = Array{Int64}(m+1)
		c[1]=1
		numpairs = [0]
		facecount = [0]
		if sd == maxcard-1
			ff2pv = Array{Int64}(nl)
			ff2pv[:] = m+1
		end
		if sd == maxcard
			#### sort j-matrix by grain
			alterweight = Array{Int64}(length(zll));
			maxweight = maximum(zll);
			for i = 1:length(alterweight)
				alterweight[i] = 1+maxweight-zll[i]
			end
			lowfilt = yafterx(alterweight,jrv)
			invertiblevec = integersinsameorderbycolumn2(lowfilt,jcp)
			inversevec0 = Array{Int64}(nl)
			inversevec0[invertiblevec]=1:nl
			jrv = yafterx(jrv,inversevec0)
			jz = yafterx(jz,inversevec0)

			lowfilt = yafterx(ff2pv,jrv)
			invertiblevec = integersinsameorderbycolumn2(lowfilt,jcp)
			inversevec1 = Array{Int64}(nl)
			inversevec1[invertiblevec]=1:nl
			jrv = yafterx(jrv,inversevec1)
			jz = yafterx(jz,inversevec1)
			translatorvecb = yafterx(inversevec0,inversevec1)
			inversevec0 = [];inversevec1 = [];lowfilt = [];invertiblevec = [];
			#gc()
			(rt,ct,zt) = transposeLighter(jrv,jcp,jz,nll)
			colsum = ct-1

			pmhist = zeros(Int64,m+1,m) #for sloth (apologies) we'll leave some unsed stuff in row m+1
			fpi = zeros(Int64,m+1,m)
			processfpi!(pmhist,fpi,jcp,jrv,ff2pv,m)

			#### reset ff2pv for next round
			ff2pvold = copy(ff2pv)
			ff2pv = Array{Int64}(nl)
			ff2pv[:] = m+1

			oldclaw = Array{Int64}(m)
		end

		for i = 1:m
			izfull[:]=0
			lrange = cran(jcp,i)
			izfull[jrv[lrange]] = jz[lrange]

			for j = (i+1):m
				dij = symmat[j,i]
				if dij == 0
					continue
				end
				if sd < maxcard-1
					process_sd_lt_maxcard!(
						i::Int64,j::Int64,dij::Int64,stepsize::Int64,
						facecount::Array{Int64,1},numpairs::Array{Int64,1},
						jrv::Array{Int64,1},jcp::Array{Int64,1},jz::Array{Int64,1},
						r::Array{Int64,1},z::Array{Int64,1},pflist::Array{Int64,1},
						izfull::Array{Int64,1},ff2pv::Array{Int64,1},pmhist::Array{Int64,2},
						npsupp::BitArray{1})
				elseif sd == maxcard-1
					process_sd_onelt_maxcard_1!(
						i::Int64,j::Int64,dij::Int64,stepsize::Int64,
						facecount::Array{Int64,1},numpairs::Array{Int64,1},
						jrv::Array{Int64,1},jcp::Array{Int64,1},jz::Array{Int64,1},
						r::Array{Int64,1},z::Array{Int64,1},pflist::Array{Int64,1},
						izfull::Array{Int64,1},ff2pv::Array{Int64,1},pmhist::Array{Int64,2},
						npsupp::BitArray{1})
				else
					for l = 1:(i-1)
						oldclaw[l] = minimum(symmat[l,[i,j]])
					end
					process_maxcard_one2i!(
						i::Int64,j::Int64,dij::Int64,stepsize::Int64,
						facecount::Array{Int64,1},numpairs::Array{Int64,1},
						jrv::Array{Int64,1},jcp::Array{Int64,1},jz::Array{Int64,1},
						r::Array{Int64,1},z::Array{Int64,1},pflist::Array{Int64,1},
						oldclaw::Array{Int64,1},zll::Array{Int64,1},colsum::Array{Int64,1},
						rt::Array{Int64,1},ct::Array{Int64,1},zt::Array{Int64,1},
						izfull::Array{Int64,1},ff2pv::Array{Int64,1},
						pmhist::Array{Int64,2},fpi::Array{Int64,2},
						npsupp::BitArray{1})

					process_maxcard_i2i!(
						i::Int64,j::Int64,dij::Int64,stepsize::Int64,
						facecount::Array{Int64,1},numpairs::Array{Int64,1},
						jrv::Array{Int64,1},jcp::Array{Int64,1},jz::Array{Int64,1},
						r::Array{Int64,1},z::Array{Int64,1},pflist::Array{Int64,1},
						oldclaw::Array{Int64,1},zll::Array{Int64,1},colsum::Array{Int64,1},
						rt::Array{Int64,1},ct::Array{Int64,1},zt::Array{Int64,1},
						izfull::Array{Int64,1},ff2pv::Array{Int64,1},
						pmhist::Array{Int64,2},fpi::Array{Int64,2},
						npsupp::BitArray{1})

					process_maxcard_i2j!(
						i::Int64,j::Int64,dij::Int64,stepsize::Int64,
						facecount::Array{Int64,1},numpairs::Array{Int64,1},
						jrv::Array{Int64,1},jcp::Array{Int64,1},jz::Array{Int64,1},
						r::Array{Int64,1},z::Array{Int64,1},pflist::Array{Int64,1},
						oldclaw::Array{Int64,1},zll::Array{Int64,1},colsum::Array{Int64,1},
						rt::Array{Int64,1},ct::Array{Int64,1},zt::Array{Int64,1},
						izfull::Array{Int64,1},ff2pv::Array{Int64,1},
						pmhist::Array{Int64,2},fpi::Array{Int64,2},
						npsupp::BitArray{1})

					process_maxcard_j2j!(
						i::Int64,j::Int64,dij::Int64,stepsize::Int64,
						facecount::Array{Int64,1},numpairs::Array{Int64,1},
						jrv::Array{Int64,1},jcp::Array{Int64,1},jz::Array{Int64,1},
						r::Array{Int64,1},z::Array{Int64,1},pflist::Array{Int64,1},
						oldclaw::Array{Int64,1},zll::Array{Int64,1},colsum::Array{Int64,1},
						rt::Array{Int64,1},ct::Array{Int64,1},zt::Array{Int64,1},
						izfull::Array{Int64,1},ff2pv::Array{Int64,1},
						pmhist::Array{Int64,2},fpi::Array{Int64,2},
						npsupp::BitArray{1})

					process_maxcard_j2end!(
						i::Int64,j::Int64,dij::Int64,stepsize::Int64,
						facecount::Array{Int64,1},numpairs::Array{Int64,1},
						jrv::Array{Int64,1},jcp::Array{Int64,1},jz::Array{Int64,1},
						r::Array{Int64,1},z::Array{Int64,1},pflist::Array{Int64,1},
						oldclaw::Array{Int64,1},zll::Array{Int64,1},colsum::Array{Int64,1},
						rt::Array{Int64,1},ct::Array{Int64,1},zt::Array{Int64,1},
						izfull::Array{Int64,1},ff2pv::Array{Int64,1},
						pmhist::Array{Int64,2},fpi::Array{Int64,2},
						npsupp::BitArray{1})
				end
			end
			# update the column pattern and the total number of nonzeros
			# encountered per codim2 face
			c[i+1] = facecount[1]+1
			if sd == maxcard
				colsum[jrv[cran(jcp,i)]]+=1
			end
		end
		delrange = c[end]:length(r)
		deleteat!(r,delrange)
		deleteat!(z,delrange)
		deleteat!(pflist,(numpairs[1]+1):nl)
		if sd == maxcard
			r = translatorvecb[r]
		end
		firstv[sd] = c
		farfaces[sd] = r
		prepairs[sd] = pflist
		grain[sd] = z
		if isempty(farfaces[sd])
			for nextcard = (sd+1):maxcard
				firstv[nextcard] = [1;1]
				farfaces[nextcard] = Array{Int64}(0)
				prepairs[nextcard] = Array{Int64}(0)
				grain[nextcard] = Array{Int64}(0)
			end
			if verbose
				println("no simplices of cardinality $(sd) or higher")
			end
			break
		end
	end
	if verbose
		println("collecting garbage")
		println(["number of edges" length(farfaces[2])])
	end
	#gc()
	if dictionaryoutput == true
		D = Dict(
			"farfaces" => farfaces,
			"firstv" => firstv,
			"grain" => grain,
			"prepairs" => prepairs,
			"symmat" => symmat,
			"nvl2ovl"=> vperm)
		return D
	else
		return farfaces,firstv,grain,prepairs,symmat,vperm
	end
end

function process_sd_lt_maxcard!(
	i::Int64,j::Int64,dij::Int64,stepsize::Int64,
	facecount::Array{Int64,1},numpairs::Array{Int64,1},
	jrv::Array{Int64,1},jcp::Array{Int64,1},jz::Array{Int64,1},
	r::Array{Int64,1},z::Array{Int64,1},pflist::Array{Int64,1},
	izfull::Array{Int64,1},ff2pv::Array{Int64,1},pmhist::Array{Int64,2},
	npsupp::BitArray{1})
	for k = cran(jcp,j)
		kk = jrv[k]
		farfilt = jz[k]
		if izfull[kk]>0
			claw = min(izfull[kk],dij)
			faceupdate!(facecount,r,z,k,min(farfilt,claw),stepsize)
			if claw >= farfilt && npsupp[k]
				pairupdate!(k,facecount,pflist,numpairs,npsupp,1)
			end
		end
	end
end

function process_sd_onelt_maxcard_1!(
	i::Int64,j::Int64,dij::Int64,stepsize::Int64,
	facecount::Array{Int64,1},numpairs::Array{Int64,1},
	jrv::Array{Int64,1},jcp::Array{Int64,1},jz::Array{Int64,1},
	r::Array{Int64,1},z::Array{Int64,1},pflist::Array{Int64,1},
	izfull::Array{Int64,1},ff2pv::Array{Int64,1},pmhist::Array{Int64,2},
	npsupp::BitArray{1})
	for k = cran(jcp,j)
		kk = jrv[k]
		farfilt = jz[k]
		if izfull[kk]>0
			claw = min(izfull[kk],dij)
			faceupdate!(facecount,r,z,k,min(farfilt,claw),stepsize)
			if npsupp[k] && (claw >= farfilt)
				pairupdatedeluxe!(k,i,j,numpairs,facecount,pflist,ff2pv,npsupp,pmhist)
			end
		end
	end
end

function process_maxcard_one2i!(
	i::Int64,j::Int64,dij::Int64,stepsize::Int64,
	facecount::Array{Int64,1},numpairs::Array{Int64,1},
	jrv::Array{Int64,1},jcp::Array{Int64,1},jz::Array{Int64,1},
	r::Array{Int64,1},z::Array{Int64,1},pflist::Array{Int64,1},
	oldclaw::Array{Int64,1},zll::Array{Int64,1},colsum::Array{Int64,1},
	rt::Array{Int64,1},ct::Array{Int64,1},zt::Array{Int64,1},
	izfull::Array{Int64,1},ff2pv::Array{Int64,1},
	pmhist::Array{Int64,2},fpi::Array{Int64,2},
	npsupp::BitArray{1})
	for l = 1:(i-1)
		if fpi[l,j]<fpi[l+1,j]
			ocl = oldclaw[l]
			if ocl < dij
				process_maxcard_one2i_subroutine!(
					i::Int64,j::Int64,dij::Int64,stepsize::Int64,
					facecount::Array{Int64,1},numpairs::Array{Int64,1},
					jrv::Array{Int64,1},jcp::Array{Int64,1},jz::Array{Int64,1},
					r::Array{Int64,1},z::Array{Int64,1},pflist::Array{Int64,1},
					oldclaw::Array{Int64,1},zll::Array{Int64,1},colsum::Array{Int64,1},
					rt::Array{Int64,1},ct::Array{Int64,1},zt::Array{Int64,1},
					izfull::Array{Int64,1},ff2pv::Array{Int64,1},
					pmhist::Array{Int64,2},fpi::Array{Int64,2},
					npsupp::BitArray{1},
					l::Int64,ocl::Int64)
			end
		end
	end
end

function process_maxcard_one2i_subroutine!(
	i::Int64,j::Int64,dij::Int64,stepsize::Int64,
	facecount::Array{Int64,1},numpairs::Array{Int64,1},
	jrv::Array{Int64,1},jcp::Array{Int64,1},jz::Array{Int64,1},
	r::Array{Int64,1},z::Array{Int64,1},pflist::Array{Int64,1},
	oldclaw::Array{Int64,1},zll::Array{Int64,1},colsum::Array{Int64,1},
	rt::Array{Int64,1},ct::Array{Int64,1},zt::Array{Int64,1},
	izfull::Array{Int64,1},ff2pv::Array{Int64,1},
	pmhist::Array{Int64,2},fpi::Array{Int64,2},
	npsupp::BitArray{1},
	l::Int64,ocl::Int64)
	for k = fpi[l,j]:(fpi[l+1,j]-1)
		kk = jrv[k]	## may have to reindex this
		farfilt = jz[k]
		if zll[kk] <= ocl
			break
		elseif oldclaw[l] < min(farfilt,izfull[kk])
			claw = min(izfull[kk],dij)
			if claw >= farfilt
				if npsupp[k]
					faceupdate!(facecount,r,z,k,farfilt,stepsize)
					pairupdate!(k,facecount,pflist,numpairs,npsupp,3)
					ff2pv[k] = i
				elseif oldclaw[ff2pv[k]]>=farfilt
					continue
				elseif saveface(ct,kk,colsum,farfilt,oldclaw,rt,zt)
					faceupdate!(facecount,r,z,k,farfilt,stepsize)
				end
			elseif (claw>0) && saveface(ct,kk,colsum,claw,oldclaw,rt,zt)
				faceupdate!(facecount,r,z,k,claw,stepsize)
			end
		end
	end
end

function process_maxcard_i2i!(
	i::Int64,j::Int64,dij::Int64,stepsize::Int64,
	facecount::Array{Int64,1},numpairs::Array{Int64,1},
	jrv::Array{Int64,1},jcp::Array{Int64,1},jz::Array{Int64,1},
	r::Array{Int64,1},z::Array{Int64,1},pflist::Array{Int64,1},
	oldclaw::Array{Int64,1},zll::Array{Int64,1},colsum::Array{Int64,1},
	rt::Array{Int64,1},ct::Array{Int64,1},zt::Array{Int64,1},
	izfull::Array{Int64,1},ff2pv::Array{Int64,1},
	pmhist::Array{Int64,2},fpi::Array{Int64,2},
	npsupp::BitArray{1})
	for k = fpi[i,j]:(fpi[i+1,j]-1)
		kk = jrv[k]
		farfilt = jz[k]
		if dij >= farfilt && npsupp[k]
			faceupdate!(facecount,r,z,k,farfilt,stepsize)
			pairupdate!(k,facecount,pflist,numpairs,npsupp,4)
			ff2pv[k] = i
		else
			farfilt = min(dij,farfilt)
			if saveface(ct,kk,colsum,farfilt,oldclaw,rt,zt)
				faceupdate!(facecount,r,z,k,farfilt,stepsize)
			end
		end
	end
end

function process_maxcard_i2j!(
	i::Int64,j::Int64,dij::Int64,stepsize::Int64,
	facecount::Array{Int64,1},numpairs::Array{Int64,1},
	jrv::Array{Int64,1},jcp::Array{Int64,1},jz::Array{Int64,1},
	r::Array{Int64,1},z::Array{Int64,1},pflist::Array{Int64,1},
	oldclaw::Array{Int64,1},zll::Array{Int64,1},colsum::Array{Int64,1},
	rt::Array{Int64,1},ct::Array{Int64,1},zt::Array{Int64,1},
	izfull::Array{Int64,1},ff2pv::Array{Int64,1},
	pmhist::Array{Int64,2},fpi::Array{Int64,2},
	npsupp::BitArray{1})
	for k = fpi[i+1,j]:(fpi[j,j]-1)
		kk = jrv[k]
		if izfull[kk]>0
			farfilt = jz[k]
			claw = min(izfull[kk],dij)
			if claw >= farfilt && npsupp[k]
				faceupdate!(facecount,r,z,k,farfilt,stepsize)
				pairupdate!(k,facecount,pflist,numpairs,npsupp,5)
				ff2pv[k] = i
			else
				farfilt = min(claw,farfilt)
				if saveface(ct,kk,colsum,farfilt,oldclaw,rt,zt)
					faceupdate!(facecount,r,z,k,farfilt,stepsize)
				end
			end
		end
	end
end

function process_maxcard_j2j!(
	i::Int64,j::Int64,dij::Int64,stepsize::Int64,
	facecount::Array{Int64,1},numpairs::Array{Int64,1},
	jrv::Array{Int64,1},jcp::Array{Int64,1},jz::Array{Int64,1},
	r::Array{Int64,1},z::Array{Int64,1},pflist::Array{Int64,1},
	oldclaw::Array{Int64,1},zll::Array{Int64,1},colsum::Array{Int64,1},
	rt::Array{Int64,1},ct::Array{Int64,1},zt::Array{Int64,1},
	izfull::Array{Int64,1},ff2pv::Array{Int64,1},
	pmhist::Array{Int64,2},fpi::Array{Int64,2},
	npsupp::BitArray{1})
	for k = fpi[j,j]:(fpi[j+1,j]-1)
		kk = jrv[k]
		if izfull[kk]>0
			claw = min(izfull[kk],dij)
			if saveface(ct,kk,colsum,claw,oldclaw,rt,zt)
				faceupdate!(facecount,r,z,k,claw,stepsize)
			end
		end
	end
end

function process_maxcard_j2end!(
	i::Int64,j::Int64,dij::Int64,stepsize::Int64,
	facecount::Array{Int64,1},numpairs::Array{Int64,1},
	jrv::Array{Int64,1},jcp::Array{Int64,1},jz::Array{Int64,1},
	r::Array{Int64,1},z::Array{Int64,1},pflist::Array{Int64,1},
	oldclaw::Array{Int64,1},zll::Array{Int64,1},colsum::Array{Int64,1},
	rt::Array{Int64,1},ct::Array{Int64,1},zt::Array{Int64,1},
	izfull::Array{Int64,1},ff2pv::Array{Int64,1},
	pmhist::Array{Int64,2},fpi::Array{Int64,2},
	npsupp::BitArray{1})
	for k = fpi[j+1,j]:(jcp[j+1]-1)
		kk = jrv[k]
		if izfull[kk]>0
			farfilt = jz[k]
			claw = min(izfull[kk],dij)
			if claw >= farfilt && npsupp[k]
				faceupdate!(facecount,r,z,k,farfilt,stepsize)
				pairupdate!(k,facecount,pflist,numpairs,npsupp,6)
				ff2pv[k] = i
			else
				farfilt = min(claw,farfilt)
				if saveface(ct,kk,colsum,farfilt,oldclaw,rt,zt)
					faceupdate!(facecount,r,z,k,farfilt,stepsize)
				end
			end
		end
	end
end

function pairupdate!(k::Int64,facecount::Array{Int64,1},pflist::Array{Int64,1},numpairs::Array{Int64,1},npsupp::BitArray{1},iterateNumber)
	numpairs[1]+=1
	pflist[numpairs[1]] = facecount[1]
	npsupp[k]=false
end

function pairupdatedeluxe!(k::Int64,i::Int64,j::Int64,numpairs::Array{Int64,1},facecount::Array{Int64,1},pflist::Array{Int64,1},ff2pv::Array{Int64,1},npsupp::BitArray{1},pmhist::Array{Int64,2})
	numpairs[1]+=1
	pmhist[i,j]+=1
	npsupp[k]=false
	pflist[numpairs[1]]=facecount[1]
	ff2pv[k] = i
end

function faceupdate!(facecount::Array{Int64,1},r::Array{Int64,1},z::Array{Int64,1},k::Int64,farfilt::Int64,stepsize::Int64)
	facecount[1]+=1
	if facecount[1]>length(r)
		append!(r,Array{Int64}(stepsize))
		append!(z,Array{Int64}(stepsize))
	end
	r[facecount]= k
	z[facecount]= farfilt
end

function faceupdatedeluxe!(facecount::Array{Int64,1},r::Array{Int64,1},z::Array{Int64,1},k::Int64,farfilt::Int64,stepsize::Int64,s::Array{Int64,1},i::Int64)
	facecount[1]+=1
	if facecount[1]>length(r)
		append!(r,Array{Int64}(stepsize))
		append!(z,Array{Int64}(stepsize))
		append!(s,Array{Int64}(stepsize))
	end
	r[facecount]= k
	z[facecount]= farfilt
	s[facecount]= i
end

function saveface(ct::Array{Int64,1},kk::Int64,colsum::Array{Int64,1},farfilt::Int64,oldclaw::Array{Int64,1},rt::Array{Int64,1},zt::Array{Int64,1})
	keep = true
	for l = ct[kk]:colsum[kk]
		if  zt[l]>= farfilt && oldclaw[rt[l]]>=farfilt
			keep = false
			break
		end
	end
	return keep
end

function savefacespecial(ct::Array{Int64,1},kk::Int64,colsum::Array{Int64,1},farfilt::Int64,oldclaw::Array{Int64,1},rt::Array{Int64,1},zt::Array{Int64,1})
	keep = true
	for l = ct[kk]:colsum[kk]
		if  zt[l]>= farfilt && oldclaw[rt[l]]>=farfilt
			keep = false
			println(["rt[l]" rt[l]])
			break
		end
	end
	if keep
		println("kept")
	end
	return keep
end

function processfpi!(pmhist::Array{Int64,2},fpi::Array{Int64,2},jcp::Array{Int64,1},jrv::Array{Int64,1},ff2pv::Array{Int64,1},m::Integer)
	for p = 1:m
		for q = jcp[p]:(jcp[p+1]-1)
			pmhist[ff2pv[jrv[q]],p]+=1
		end
	end
	for p = 1:m
		fpi[1,p] = jcp[p]
		for q = 1:m
			fpi[q+1,p] = fpi[q,p]+pmhist[q,p]
		end
	end
end

function generate2faces(symmat)
	m = size(symmat,1)
	if issparse(symmat)
		return symmat
	else
		L = 0
		for i = 1:m
			for j = (i+1):m
				if symmat[j,i]>0
					L+=1
				end
			end
		end
		rowval = Array{Int64}(L)
		nzval = Array{Int64}(L)
		colptr = Array{Int64}(m+1)
		marker = 0
		colptr[1] = 1
		for i = 1:m
			colptr[i+1]=colptr[i]
			for j = (i+1):m
				if symmat[j,i]>0
					colptr[i+1]+=1
					rowval[colptr[i+1]-1] = j
					nzval[colptr[i+1]-1] = symmat[j,i]
				end
			end
		end
	end
	return rowval,colptr,nzval
end

function generate3faces!(
	farfaces_cell,
	firstv_cell,
	grain_cell,
	prepairs_cell,
	m,
	symmat;
	verbose = false)

	grain::Array{Int64,1} = grain_cell[2]
	farfaces::Array{Int64,1} = farfaces_cell[2]
	firstv::Array{Int64,1} = firstv_cell[2]

	numverts = length(firstv)-1
	numedges = length(farfaces)
	stepsize = 10^7
	facecount= [0]
	numpairs = 0

	closefaces = Array{Int64}(numedges)
	for i = 1:m
		closefaces[cran(firstv,i)]=i
	end
	iso = integersinsameorder(farfaces)
	closefaces_higsorted = 	Array{Int64}(numedges)
	grain_higsorted = 	Array{Int64}(numedges)
	closefaces_higsorted[iso] = closefaces
	grain_higsorted[iso] = grain

	firstv_hs = zeros(Int64,m+1)
	for face in farfaces
		firstv_hs[face+1]+=1
	end
	firstv_hs[1] = 1
	for i = 2:(m+1)
		firstv_hs[i] = firstv_hs[i-1]+firstv_hs[i]
	end

	adist = Array{Int64}(m)
	idist = Array{Int64}(m)
	r = Array{Int64}(numedges)
	z = Array{Int64}(numedges)
	s = Array{Int64}(numedges)

	clawvec = Array{Int64}(m)
	ncheckedges = trues(numedges)

	for a = 1:m
		adist[:]=0
		adist[crows(firstv,farfaces,a)] = crows(firstv,grain,a)
		for ip = cran(firstv,a)
			i = farfaces[ip]
			dai = grain[ip]
			idist[:]=0
			idist[crows(firstv,farfaces,i)]	= crows(firstv,grain,i)
			idist[crows(firstv_hs,closefaces_higsorted,i)] = crows(firstv_hs,grain_higsorted,i)
			for jp = cran(firstv,i)
				if ncheckedges[jp]
					j = farfaces[jp]
					dij = grain[jp]
					if dij <= dai && dij <= adist[j] # note this condition bakes in the req. that j be adjacent to a
						numpairs+=1
						ncheckedges[jp] = false
						clawvec[1:i] = 0
						for lp = cran(firstv_hs,j)
							l = closefaces_higsorted[lp]
							if l >= i
								break
							elseif idist[l]!=0
								clawvec[l] = min(idist[l],grain_higsorted[lp])
							end
						end
						for kp = cran(firstv,j)
							k = farfaces[kp]
							djk = grain[kp]
							dak = adist[k]
							dik = idist[k]
							if dak < dij && dak<djk && dak<dik	# this bakes in req. that dik>0
								dijk = min(dij,dik,djk)
								keepface = true
								for bp = cran(firstv_hs,k)
									b = closefaces_higsorted[bp]
									if b >= i
										break
									elseif min(clawvec[b],grain_higsorted[bp]) >= dijk
										keepface = false
										break
									end
								end
								if keepface
									faceupdatedeluxe!(facecount,r,z,kp,dijk,stepsize,s,i)
								end
							end
						end
					end
				end
			end
		end
	end
	holdi = 0
	for edge = find(ncheckedges)
		i = closefaces[edge]
		j = farfaces[edge]
		dij = grain[edge]
		if i != holdi
			idist[:]=0
			idist[crows(firstv,farfaces,i)]	= crows(firstv,grain,i)
			idist[crows(firstv_hs,closefaces_higsorted,i)] = crows(firstv_hs,grain_higsorted,i)
			holdi = i
		end
		clawvec[1:i] = 0
		for lp = cran(firstv_hs,j)
			l = closefaces_higsorted[lp]
			if l >= i
				break
			elseif idist[l]!=0
				clawvec[l] = min(idist[l],grain_higsorted[lp])
			end
		end
		#### a facsimile of above
		for kp = cran(firstv,j)
			k = farfaces[kp]
			dik = idist[k]
			if dik==0
				continue
			end
			djk = grain[kp]
			dijk = min(dij,dik,djk)
			keepface = true
			for bp = cran(firstv_hs,k)
				b = closefaces_higsorted[bp]
				if b >= i
					break
				elseif min(clawvec[b],grain_higsorted[bp]) >= dijk
					keepface = false
					break
				end
			end
			if keepface
				faceupdatedeluxe!(facecount,r,z,kp,dijk,stepsize,s,i)
			end
		end
		####
	end

	num3faces = facecount[1]
	holderlengths = length(r)
	deletionrange = (num3faces+1):holderlengths
	deleteat!(r,deletionrange)
	deleteat!(z,deletionrange)
	deleteat!(s,deletionrange)

	iso = integersinsameorder(s)
	r[iso]=r
	z[iso]=z
	fv3 = zeros(Int64,numverts+1)
	fv3[1] = 1
	for face in s
		fv3[face+1]+=1
	end
	for i = 2:(numverts+1)
		fv3[i] = fv3[i-1]+fv3[i]
	end

	pairmarker = 0
	npes = trues(numedges)
	prepairs = Array{Int64}(numedges)
	for i = 1:num3faces
		edge = r[i]
		if npes[edge] && z[i] == grain[edge]
			npes[edge]=false
			pairmarker+=1
			prepairs[pairmarker]=i
		end
	end
	deleteat!(prepairs,(pairmarker+1):numedges)

	farfaces_cell[3]=r
	firstv_cell[3]=fv3
	grain_cell[3]=z
	prepairs_cell[3]=prepairs

	buffer1 = Array{Array{Int64,1},1}(1)
	buffer2 = Array{Array{Int64,1},1}(1)
	buffer1[1] = Array{Int64}(0)
	buffer2[1] = ones(Int64,numverts+1)
	append!(farfaces_cell,buffer1)
	append!(grain_cell,buffer1)
	append!(prepairs_cell,buffer1)
	append!(firstv_cell,buffer1)

	return r,fv3,z,prepairs,numpairs
end

##########################################################################################

####	BATCH OPERATIONS

##########################################################################################

function eirene_batchcsv(
	inputdirectory,
	outputdirectory;
	maxdim = 1,
	model="dmat",
	filetype="textfile",
	lowerlim=-Inf,
	upperlim=Inf,
	numrad=Inf,
	fastop=true,
	record="cyclerep",
	pointlabels=[],
	verbose=false)

	filenames = readdir(inputdirectory)

	for i = 2:length(filenames)
		filename = filenames[i]
		filepath = "$(inputdirectory)/$(filename)"
		C = eirene(readcsv(filepath),
				maxdim 	=maxdim,
				model		=model,
				filetype	=filetype,
				lowerlim	=lowerlim,
				upperlim	=upperlim,
				numrad	=numrad,
				fastop	=fastop,
				record		=record,
				pointlabels	=pointlabels,
				verbose		=verbose)
		savepath = "$(outputdirectory)/$(filename).jld"
		JLD.save(savepath,"C",C)
	end
end

##########################################################################################

####	MAIN

##########################################################################################

"""

    eirene(X[, keyword arugemts])

Computes the persistent homology of a filtered complex.

"""
function eirene(
	s;
	filetype	= "textfile",
	model		= "vr",
	maxdim 		= 1,
	minrad		= -Inf,
	maxrad		= Inf,
	numrad		= Inf,
	fastop		= true,
	vscale		= "default",
	record		= "cyclerep",
	pointlabels	= [],
	verbose		= false)

	if in(model,["vr","pc"])
		maxcard = 		maxdim+2
		D = persistF2_vr(
			s,
			maxcard;
			model 		= model,
			minrad 		= minrad,
			maxrad 		= maxrad,
			numrad 		= numrad,
			fastop 		= fastop,
			record 		= record,
			filetype 	= filetype,
			pointlabels = pointlabels,
			verbose 	= verbose)

		return 	D
	elseif model == "complex"
		D   = 	persistF2_cell(
				s;
				maxdim=maxdim,
				record = record,
				verbose=false)
		return 	D
	else
		println()
		println("Error: the only valid values for keyword <model> are \"vr\", \"pc\", and \"cell\".")
		println("user input:")
		println(model)
	end
end

function eirene(
	rv,
	cp,
	filt;
	maxdim = length(rowvalues),
	record="cyclerep",
	pointlabels=[],
	verbose=false)

	D = persistF2_cell(
		rv,
		cp,
		filt;
		maxdim = maxdim,
		record = record,
		verbose = verbose)

	return D
end

function eirene(
	rv,
	cp,
	filt,
	dim;
	maxdim = length(rowvalues),
	record="cyclerep",
	pointlabels=[],
	verbose=false)

	D = persistF2_cell(
		rv,
		cp,
		filt,
		dim,
		maxdim = maxdim,
		record = record,
		verbose = verbose)

	return D
end

##########################################################################################

####	TESTING AND DIAGNOSTICS // FUNCTIONWISE

##########################################################################################

function persistencestats(x)
	L = length(x)
	A = Array{Float64}(8,L)
	for ip = 1:L
		i = x[ip]
		println(i)
		println(i)
		pcloud = rand(20,i)
		tic()
		D = persistF2_vr(pcloud,5,model = "pc",fastop=false)
 		t = toc()
		A[1,ip] = length(D["prepairs"][5])
		A[2,ip] = length(D["farfaces"][5])
		A[3,ip] = binom(i-1,4)
		A[4,ip] = length(D["trv"][5])
		A[5,ip] = t
		A[6,ip] = i
		A[7,ip] = length(D["farfaces"][5]) - length(D["prepairs"][5]) / binom(i-1,4)
		A[8,ip] = length(D["trv"][5]) / binom(i-1,4)
	end
	return A
end

function persistencemultistats(x;numtrials = 10)
	trials = Array{Any,1}(numtrials)
	for i = 1:numtrials
		trials[i] = persistencestats(x)
	end
	return trials
end

function construction_sanitycheck(;numtrials = 10,samplesize = 50,sd=4)
	for i = 1:numtrials
		pcloud = rand(20,samplesize)
		d = Distances.pairwise(Euclidean(),pcloud)
		(t,ocg2rad) = ordercanonicalform(d;fastop=false)
		construction_sanitycheck_subroutine(t,sd,samplesize)
		#gc()
	end
end

function construction_sanitycheck_subroutine(t,n,samplesize)
	N = n+2
	D2 = buildcomplex3(t,N,dictionaryoutput = true)
	D3 = buildcomplex3(t,n,dictionaryoutput = true)
	D2supp = trues(length(D2["farfaces"][n]))
	D2supp[D2["farfaces"][n+1][D2["prepairs"][n+1]]]=false
	checkthesearch = false
	for i = 1:samplesize
		D2ran = cran(D2["firstv"][n],i)
		D2ran = D2ran[D2supp[D2ran]]
		D2ff = D2["farfaces"][n][D2ran]
		D3ff = crows(D3["firstv"][n],D3["farfaces"][n],i)
		if sort(D2ff) != sort(D3ff)
			println("	sort(D2ff) != sort(D3ff) ")
			checkthesearch = true
			break
		end
		if checkthesearch == true
			break
		end
	end
	if sort(D2["farfaces"][n][D2["prepairs"][n]])!=sort(D3["farfaces"][n][D3["prepairs"][n]])
		println("please check prepairs")
		sleep(10)
	elseif checkthesearch
		println("checkthesearch evaluated to true")
		sleep(10)
	else
		println("ok so far as these checks are concerned :)")
	end
end

function cellcheck()
	c = 0
	for i = 1:20
		print([i])
		D 	= eirene(rand(20,50),model="pc",maxdim=4)
		N 	= ff2complex(D["farfaces"],D["firstv"])
		Nf 	= ocff2of(D["grain"],D["ocg2rad"])
		N1 = copy(N[1]);
		N2 = copy(N[2]);
		Nf_copy = copy(Nf)
		F = persistF2_cell(N[1],N[2],Nf)
		if N1 != N[1] || N2 != N[2]
			print("changed N")
			break
		elseif Nf_copy != Nf
			print("changed Nf")
			break
		end
		for k = 1:4
			if sortrows(barcode(D,dim=k))!= sortrows(barcode(F,dim=k))
				c+=1
			end
		end
	end
	return c
end

##########################################################################################

####	TESTING AND DIAGNOSTICS // SYSTEM LEVEL

##########################################################################################

# 12/30/2017
# File nameing systems.
function 	inputsuffix(model,iteration)
	if model == "complex"
		pref = "gdb_da_a_a_cell"
	elseif model == "vr"
		pref = "gdb_da_a_b_vr"
	elseif model == "pc"
		pref = "gdb_da_a_c_pc"
	end
	return pref
end

function 	pathsuffix(model,iteration)
	pref	= 	inputsuffix(model,iteration)
	w		= 	string(	pref,"/",
						pref,"$(iteration)","/",
						pref,"$(iteration)","_input.csv"
						)
	return w
end

function 	modit2filepath(model,iteration)
	# fileroot = 	"/Users/greghenselman/Google Drive/GregDirectory/julia_gd/gdb_eirene/gdb_data/gdb_da_a_calib/"
	fileroot = 	"/Users/gh10/Google Drive/gregtaylorgoogledrive_nongooglefiles/GregDirectory/julia_gd/gdb_eirene/gdb_data/gdb_da_a_calib/"
	return(string(fileroot,pathsuffix(model,iteration)))
end

# 12/30/2017
# This function runs through a sampling of different keyword parameters, and compares the
# output of the current Eirene version against that of a previous one, which has been stored
# to files.  The specific outputs checked are barcodes and (for pre-cooked examples) cycle
# reps.
function unittest()

	filepath 							= 	joinpath(@__DIR__,"solutions/testsolutions.jld")
	K 									= 	JLD.load(filepath)			# K for "calibrate" (letter C was taken)
	K									=  	K["K"]
	errorindices						= 	Array{Any,1}(0)

	for filetype 	= ["textfile"]
		for model 		= ["complex"]#["vr" "pc" "complex"]
			for maxdim 		= [0 1 2]
				for minrad 		= [-Inf 0]
					for maxrad 		= [Inf, 100]
						for numrad 		= [1 10 Inf]
							for fastop 		= [true,false]
								for vscale 		= [[]]
									for record 		=  ["all" "cyclerep" "none"]
										for pointlabels	= [[]]
											for iteration	= [1 2]

												filepath 	= 	modit2filepath(model,iteration-1)

												C = eirene(
													filepath;
													filetype	= filetype,
													model		= model,
													maxdim 		= maxdim,
													minrad		= minrad,
													maxrad		= maxrad,
													numrad		= numrad,
													fastop		= fastop,
													vscale		= vscale,
													record		= record,
													pointlabels	= pointlabels,
													verbose		= false)

												D	=	Dict(
													:filetype	=> filetype,
													:model		=> model,
													:maxdim 	=> maxdim,
													:minrad		=> minrad,
													:maxrad		=> maxrad,
													:numrad		=> numrad,
													:fastop		=> fastop,
													:vscale		=> vscale,
													:record		=> record,
													:pointlabels	=> pointlabels
													)

												# nli stands for number-level-index
												if [minrad maxrad numrad] == [-Inf Inf Inf]
													nli	=	1;
												elseif [minrad maxrad numrad] == [0 100 10]
													nli	=	2;
												else
													nli = 	0
												end

												xbc = 	(nli > 0) && in(model,["vr" "complex"]) && in(record,["cyclereps","all"]) # stands for cross-check barcode
												xcr = 	xbc && model == "complex" # stands for cross-check cycle rep

												if xbc
													xkey = undercat([model,"$(nli)","$(iteration)"])
													for r 	= 	0:maxdim
														if r+1 > length(K[xkey][:barcodes])
															println("")
															println("error 2")
															println(xkey)
															println(K[xkey][:barcodes])
															println(r)
															return C,D
														end
														Ba	= 	K[xkey][:barcodes][r+1]
														Bw  = 	barcode(C,dim=r)
														if record == "all"
															for p 	= 	1:size(Bw,1)
																rep 	= 	classrep(C,dim=r,class=p,format="index")
																case1 	= 	birthtime(C,dim=r,chain=rep) == Bw[p,1]
																case2 	= 	deathtime(C,dim=r,chain=rep) == Bw[p,2]
																case3 	= 	isempty(chainboundary(C,dim=r,chain=rep))
																if !(case1 & case2 & case3)
																	println()
																	println("error: please check birth and death times")
																end
															end
														end
														if sortrows(Ba) != sortrows(Bw)
															println("[]") ##############
															println("error 1")
															println(Ba)##############
															println(Bw)##############
															append!(errorindices,[xkey])
															return C,D,K,xkey
														end
														if xcr 	& 	size(Ba,1) 	==	1
															Ra 	= 	K[xkey][:cyclerep][r+1]
															Rb	= 	classrep(C,dim=r,class=1)
															if 	sort(Ra) 	   != 	sort(Rb)
																println("")
																println("error 3")
																println(sort(Ra))
																println(sort(Rb))
																println(K[xkey]) ##############
																println([r])
																append!(errorindices,[xkey])
																return C,r,xkey,D
															end
														end
													end
												end
											end
										end
									end
								end
							end
						end
					end
				end
			end
		end
	end
	return errorindices
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

# 12/28/2017
# This function is meant to generate crosscheck data for the version of <unittest>
# defined 12/30/2017.
function generatecrosscheckdata_perseus()

	K 	= 	Dict()			# K for "calibrate" (letter C was taken)

	for filetype 	= ["textfile"]
		for model 		= ["vr" "pc"] # cellular examples must be handled separately
			for maxdim 		= [0 1 2]
				for minrad 		= [-Inf 0]
					for maxrad 		= [Inf, 100]
						for numrad 		= [1 10 Inf]
							for fastop 		= [true,false]
								for vscale 		= [[]]
									for record 		= ["all" "cyclerep" "none"]
										for pointlabels	= [[]]
											for iteration	= [1 2]

												filepath 	= 	modit2filepath(model,iteration-1)

												# 	nli stands for number-level-index
												if [minrad maxrad numrad] == [-Inf Inf Inf]
													nli	=	1;
												elseif [minrad maxrad numrad] == [0 100 10]
													nli	=	2;
												else
													nli = 	0
												end

												#		xbc stands for cross-check barcode
												#		xcr	stands for cross-check cycle rep
												xbc = 	(nli > 0) &&
														in(model,["vr" "complex"]) &&
														in(record,["cyclereps","all"]) &&
														(maxdim == 2)
												xcr = 	xbc && model == "complex"

												if xbc
													if 		numrad 	==	10
														stepsz 	= 	10
													elseif 	numrad 	== 	Inf
														stepsz	= 	1
													end

													E	= 	perseusvrjl(
															filepath;				# 	filepaths should end with .txt
															model					= 	"vr",
															rowsare 				= 	"distances",
															# datapath				= 	"/Users/greghenselman/Julia/testtext5.txt",
															datapath 				= 	"/Users/gh10/Google Drive/gregtaylorgoogledrive_nongooglefiles/GregDirectory/julia_gd/gdc_agora/gdc_a_peresuswrappers/perseustemporaydatasilo.txt.rtf",
															outpath					= 	"perseusoutput_caliber",
															maxdim 					= 	maxdim,
															minrad					= 	0,
															stepsz					= 	stepsz,
															nsteps					= 	Inf,
															# perseusfilepath 		= 	"/Users/greghenselman/Downloads/perseusMac")
															perseusfilepath 		= 	"/Users/gh10/Google Drive/gregtaylorgoogledrive_nongooglefiles/GregDirectory/julia_gd/gdc_agora/gdc_a_peresuswrappers/perseusMac")

													xkey 					= 	undercat([model,"$(nli)","$(iteration)"])
													K[xkey]					=	Dict()
													K[xkey][:barcodes]		= 	Array{Any,1}(maxdim+1)
													K[xkey][:cyclerep]		= 	Array{Any,1}(maxdim+1)  # this will be filled in manually
													for r 	= 	1:(maxdim+1)
														K[xkey][:barcodes][r] 		= perseus_barcode(E,dim=r-1)
													end
												end
											end
										end
									end
								end
							end
						end
					end
				end
			end
		end
	end
	return K
end

function manualcycleadditions(K)

    # SPHERE / MACHINE PRECISION

	maxdim		= 	2
	model 		= 	"complex"
	nli			=	1
	iteration	=   2
	xkey 		= 	undercat([model,"$(nli)","$(iteration)"])

	K[xkey]					=	Dict()
	K[xkey][:barcodes]		= 	Array{Any,1}(maxdim+1)
	K[xkey][:cyclerep]		= 	Array{Any,1}(maxdim+1)  # this will be filled in manually

	K[xkey][:barcodes][1] 	=	[20.0 30.0;10.0 Inf]
	K[xkey][:barcodes][2] 	=	[40.0 50.0]
	K[xkey][:barcodes][3] 	=	[60.0 70.0]

	K[xkey][:cyclerep][1]	= 	[1; 2]
	K[xkey][:cyclerep][2]	= 	[1; 2]
	K[xkey][:cyclerep][3]	= 	[1; 2]

    # SPHERE / 10:100 PRECISION
    # (identical to sphere with machine precision, just a different nli & key value)

	maxdim		= 	2
	model 		= 	"complex"
	nli			=	2
	iteration	=   2
	xkey 		= 	undercat([model,"$(nli)","$(iteration)"])

	K[xkey]					=	Dict()
	K[xkey][:barcodes]		= 	Array{Any,1}(maxdim+1)
	K[xkey][:cyclerep]		= 	Array{Any,1}(maxdim+1)  # this will be filled in manually

	K[xkey][:barcodes][1] 	=	[20.0 30.0;10.0 Inf]
	K[xkey][:barcodes][2] 	=	[40.0 50.0]
	K[xkey][:barcodes][3] 	=	[60.0 70.0]

	K[xkey][:cyclerep][1]	= 	[1; 2]
	K[xkey][:cyclerep][2]	= 	[1; 2]
	K[xkey][:cyclerep][3]	= 	[1; 2]

    # EMPTY SPACE / MACHINE PRECISION
	maxdim		= 	2

	model 		= 	"complex"
	nli			=	1
	iteration	=   1
	xkey 		= 	undercat([model,"$(nli)","$(iteration)"])

	K[xkey]					=	Dict()
	K[xkey][:barcodes]		= 	Array{Any,1}(maxdim+1)
	K[xkey][:cyclerep]		= 	Array{Any,1}(maxdim+1)  # this will be filled in manually

	for p = 1:maxdim+1
		K[xkey][:barcodes][p] 	=	Array{Int64,2}(0,2)
	end

    # EMPTY SPACE / 10:100 PRECISION
	maxdim		= 	2

	model 		= 	"complex"
	nli			=	2
	iteration	=   1
	xkey 		= 	undercat([model,"$(nli)","$(iteration)"])

	K[xkey]					=	Dict()
	K[xkey][:barcodes]		= 	Array{Any,1}(maxdim+1)
	K[xkey][:cyclerep]		= 	Array{Any,1}(maxdim+1)  # this will be filled in manually

	for p = 1:maxdim+1
		K[xkey][:barcodes][p] 	=	Array{Int64,2}(0,2)
	end
end

function savesolutions()
	print(
	"""

	====================================================

	Please note: the file

	/Users/gh10/Google Drive/gregtaylorgoogledrive_nongooglefiles/GregDirectory/julia_gd/gdc_agora/gdc_a_peresuswrappers/gdc_a_a_perseuswrapper.jl

	must be loaded for this function to work properly.

	====================================================
	""")
	K 	= 	generatecrosscheckdata_perseus()
	manualcycleadditions(K)
	filepath = joinpath(@__DIR__,"solutions/testsolutions.jld")
	JLD.save(filepath,"K",K)
	# JLD.save("/Users/gh10/Google Drive/gregtaylorgoogledrive_nongooglefiles/GregDirectory/julia_gd/gdb_eirene/gdb_data/gdb_da_a_calib/gdb_da_a_d_solns/gdb_da_a_dsoln_a.jld","K",K)
	# JLD.save("/Users/greghenselman/Google Drive/GregDirectory/julia_gd/gdb_eirene/gdb_data/gdb_da_a_calib/gdb_da_a_d_solns/gdb_da_a_dsoln_a.jld","K",K)
	return K
end

##########################################################################################

####	PERSEUS WRAPPER

##########################################################################################

#=
This file is a simple Julia wrapper for the Perseus library for persistent
homology, which is produced and maintained by Vidit Nanda:

    http://people.maths.ox.ac.uk/nanda/perseus/

Note that keyword argument <filepath> should be a filepath to an executable
form of Perseus.
=#

function perseusvrjl(
					datum;					# 	filepaths should end with .txt
					model					= 	"vr",
					rowsare 				= 	"distances",
					datapath				= 	"/Users/greghenselman/Julia/testtext4.txt",
					outpath					= 	"perseusoutput_caliber",
					maxdim 					= 	1,
					minrad					= 	0,
					stepsz					= 	0.1,
					nsteps					= 	10,
					fr						= 	false, # stands for full resolution
					perseusfilepath 		= 	"/Users/gh10/Google Drive/gregtaylorgoogledrive_nongooglefiles/GregDirectory/julia_gd/gdc_agora/gdc_a_peresuswrappers/perseusMac"
					)

	if model == "vr"
		if typeof(datum)==String
			filename = datum
			if typeof(readdlm(filename,','))<:Array{Float64}
				s = readdlm(filename,',')
			elseif typeof(readdlm(filename,' '))<:Array{Float64}
				s = readdlm(filename,' ')
			else
				print("Error reading text file.  datum files of this format must be either comma or space delimited.")
				return
			end
		elseif typeof(datum)<:Array{Int64} || typeof(datum)<:Array{Float64}
			s = datum
		else
			print("Error: datum must be a string or an array of type Array{Float64} or Array{Int64}")
		end

		if rowsare == "dimensions"
			d = Distances.pairwise(Euclidean(),s)
		elseif rowsare == "points"
			d = Distances.pairwise(Euclidean(),s')
		elseif rowsare == "distances"
			d = convert(Array{Float64},s)
		end

		if fr
			d,ocg2rad 	= 	ordercanonicalform_3{Tv}(
							d;
							maxrad=Inf,
							minrad=-Inf,
							numrad=Inf,
							vscale="default",
							)
			d 			= 	maximum(d)-d
			nsteps		= 	Inf
			stepsz		= 	1
			minrad		= 	0
			maxrad		= 	Inf
		end

		if nsteps == Inf
			nsteps = 1 + ceil(Int64,maximum(d)/stepsz)
		end

		samplesize = size(d,1)

		datastream = open(datapath,"w+")
		close(datastream)  				# this clears the current data file
		datastream = open(datapath,"a+")
		str = "$(samplesize)\n$(minrad) $(stepsz) $(nsteps) $(maxdim+1)"
		write(datastream,str)
		for p = 1:samplesize
			str = string(d[p,:])
			str = replace(str,"[","")
			str = replace(str,"]","")
			str = replace(str,",","")
			str = "\n"*str
			write(datastream,str)
		end
		close(datastream)
	elseif model == "perseusdistmat"
		datapath = datum
	end

	command  = `$perseusfilepath distmat $datapath $outpath`
	run(`$(command)`)

	D = Dict(
			:barcodes => 	Array{Array{Float64,2},1}(maxdim+1),
			:betti   => 	Array{Int64}(0,maxdim+1)
			)

	# si stands for shifted index; since peresus ouptput starts indexing at 0 and arrays are 1-indexed, the true barcode of dimension r is D[:filtvalssi][D[:barcodes][r]+1]
	if fr
		D[:filtvalssi] 	=	flipdim(ocg2rad,1)
	else
		D[:filtvalssi] 	=	minrad:stepsz:(minrad+(1+nsteps)*stepsz)
	end

	g = open(outpath*"_betti.txt")  	# i check that the file is not empty in order to avoid throwing an error with readdlm
	seekend(g)
	if position(g) != 0
		D[:betti]	= 	readdlm(outpath*"_betti.txt")
	else
		D[:betti] 	=	Array{Int64,2}(0,maxdim+1)
	end

	for p = 1:maxdim+1
		g = open(outpath*"_$(p-1).txt")  	# i check that the file is not empty in order to avoid throwing an error with readdlm
		seekend(g)
		if position(g) != 0
			D[:barcodes][p] = readdlm(outpath*"_$(p-1).txt")
		else
			D[:barcodes][p] = Array{Float64,2}(0,2)
		end
	end

	return D
end

function perseus_barcode(D;dim=1)
	p = dim+1
	if !haskey(D,:barcodes)
		print("Error: input dictionary does not appear to contain barcode data.")
		return
	elseif p > length(D[:barcodes])
		return zeros(Int64,0,2)
	else
		B 					=	D[:barcodes][p]
		B					=	2+round.(Int64,B)
		translator 			=	Array{Float64,1}(1+length(D[:filtvalssi]))
		translator[2:end]	=	D[:filtvalssi]
		translator[1]		=	Inf
		return					translator[B]
	end
end

end # module

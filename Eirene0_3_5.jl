# 	  Eirene Library for Homological Algebra
# 	  Copyright (C) 2016, 2017  Gregory Henselman
#	  www.gregoryhenselman.org
# 	
#     This file is part of the Eirene Library for Homological Algebra (Eirene).
# 
#     Eirene is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     Eirene is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with Eirene.  If not, see <http://www.gnu.org/licenses/>.

print("\n
Eirene Library for Homological Algebra
Copyright (C) 2016, 2017  Gregory Henselman
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
v0.3.5                 	
				
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

# if typeof(Pkg.installed("ParallelAccelerator")) == Void
# usepa = false
# warn("ParallelAccelerator.jl may not be installed. This package can
# significantly enhance performance. To install, enter the foll-
# owing at the Julia prompt
# 
# Pkg.add(\"ParallelAccelerator\")
# using ParallelAccelerator \n\n
# 
# and reload Eirene.")
# else
# 	using ParallelAccelerator
# end
#### @acc begin

##########################################################################################

#### 	SIMPLICIAL CONSTRUCTIONS

##########################################################################################

function vertexrealization(farfaces,firstv,facecardinality,facenames)
	numfaces::Int64 = length(facenames)
	fc::Int64 = facecardinality	
	m::Int64 = length(firstv[2])-1
	preallocationspace = 0	
	loci::Array{Int64,1} = copy(facenames)
	vrealization = Array(Int64,facecardinality,numfaces)	
	post0::Int64 = 1
	post1::Int64 = 1	
	
	for setcard = facecardinality:-1:1
		cp::Array{Int64,1} = firstv[setcard]
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
			loci[i] = farfaces[setcard][locus]
			vrealization[setcard,i]=post0
		end
	end		
	return vrealization
end

function vertexrealization(D::Dict,facecardinality,facenames)
	return vertexrealization(D["farfaces"],D["firstv"],facecardinality,facenames)
end

function vertexrealization(D::Dict;dim = 1, class = 1)
	setcard = dim+2
	facecard = dim+1
	
	if haskey(D,"cyclerep")
		rep = D["cyclerep"][setcard][class]
	else
		cyclename = barname2cyclename(D,class;dim = dim)	
		rep = getcycle(D,setcard,cyclename)		
	end
	
	vrealization = vertexrealization(D::Dict,facecard,rep)
end

function incidentverts(farfaces,firstv,facecardinality,facenames)
	numfaces::Int64 = length(facenames)
	fc::Int64 = facecardinality	
	m::Int64 = length(firstv[2])-1
	preallocationspace = 0	
	vsupp = falses(m)
	loci::Array{Int64,1} = copy(facenames)
	post0::Int64 = 1
	post1::Int64 = 1	
	for setcard = facecardinality:-1:1
		cp::Array{Int64,1} = firstv[setcard]
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
			loci[i] = farfaces[setcard][locus]
			vsupp[post0]=true
		end
	end	
	return find(vsupp)
end

function incidentverts(D::Dict,facecardinality,facenames)
	return incidentverts(D["farfaces"],D["firstv"],facecardinality,facenames)
end

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
	hclosefaces = Array(Int64,facecard,n)
	if n == 0
		return hclosefaces
	else
		rowdepth = facecard-1
		rosettacol = Array(Int64,maximum(lrowval))
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

function buildclosefromclose(
	lrowval,
	lcolptr,
	lclosefaces,
	hrowval,
	hcolptr,
	selectedcolumnindices;
	verbose=false,
	facecard = size(lclosefaces,1)+1)
	if verbose
		print("PLEASE NOTE: COLUMNS MUST BE IN SORTED ORDER FOR THIS TO WORK PROPERLY")
	end
	m = length(hcolptr)-1
	numhigs = length(hrowval)
	numselected = length(selectedcolumnindices)
	hclosefaces = Array(Int64,facecard,numselected)
	if numselected == 0
		return hclosefaces
	else
		rowdepth = facecard-1
		rosettacol = Array(Int64,maximum(lrowval))
		columnsupp = falses(numhigs)
		columnsupp[selectedcolumnindices]=true
		columnmarker = 0
		buildclosefromclose_subr(rosettacol::Array{Int64,1},lrowval::Array{Int64,1},lcolptr::Array{Int64,1},hrowval::Array{Int64,1},hcolptr::Array{Int64,1},hclosefaces::Array{Int64,1},columnmarker::Int64,rowdepth::Integer)
		return hclosefaces	
	end
end

function  buildclosefromclose_subr(rosettacol::Array{Int64,1},lrowval::Array{Int64,1},lcolptr::Array{Int64,1},hrowval::Array{Int64,1},hcolptr::Array{Int64,1},hclosefaces::Array{Int64,1},columnmarker::Int64,rowdepth::Integer)
	for i = 1:m	
		rosettacol[lrowval[cran(lcolptr,i)]]=cran(lcolptr,i)
		for j = cran(hcolptr,i)
			if columnsupp[j]		
				columnmarker+=1
				farface = hrowval[j]
				buildclosefromclose_subr_subr(rowdepth::Integer,hclosefaces::Array{Int64,1},columnmarker::Int64,rowettacol::Array{Int64,1},lclosefaces::Array{Int64,1},farface::Int64)
				hclosefaces[setcard,columnmarker] = rosettacol[lrowval[farface]]
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
	setcard = rowdepth+1
	hclosefaces = Array(Int64,setcard+1,numselected)
	if numselected == 0
		return hclosefaces
	end
	rosettacol = Array(Int64,maximum(lrowval))
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
				hclosefaces[setcard,columnmarker] = rosettacol[lrowval[farface]]
				hclosefaces[setcard+1,columnmarker] = farface
			end
		end	
	end
	return hclosefaces
end

function buildclosefaces!(lrowval,lcolptr,lclosefaces,lfarfaces,hrowval,hcolptr,destinationmatrix;verbose = false)
	m = length(hcolptr)-1
	n = length(hrowval)
	rowdepth = size(lclosefaces,1)
	setcard = rowdepth+1
	rosettacol = Array(Int64,maximum(lrowval))
	for i = 1:m
		rosettacol[lrowval[cran(lcolptr,i)]]=cran(lcolptr,i)
		for j = cran(hcolptr,i)
			farface = hrowval[j]
			for k = 1:rowdepth
				destinationmatrix[k,j]=rosettacol[lclosefaces[farface]]
			end
			destinationmatrix[setcard,j] = rosettacol[lrowval[farface]]
		end	
	end
	for j = 1:n
		for i = 1:setcard
			lclosefaces[i,j]=destinationmatrix[i,j]
		end
	end
end

function buildclosefromfar(farfaces,firstv,setcard)
	m = length(firstv[1])-1
	n = length(farfaces[setcard])
	# destinationmatrix = Array(Int64,setcard,n)
	if setcard == 1
		return Array(Int64,0,m)
	end
	lclosefaces = Array(Int64,1,firstv[2][end]-1)
	for i = 1:m
		lclosefaces[cran(firstv[2],i)]=i
	end
	if setcard == 2
		return lclosefaces'
	end
	for i = 3:setcard
		lclosefaces = buildclosefromclose(farfaces[i-1],firstv[i-1],lclosefaces,farfaces[i],firstv[i];facecard=i-1)
	end
	return lclosefaces
end

function buildclosefromfar(farfaces,firstv,setcard,columnsinorder)
	m = length(firstv[1])-1
	n = length(farfaces[setcard])
	if setcard == 1
		return Array(Int64,0,m)
	end
	lclosefaces = Array(Int64,1,firstv[2][end]-1)
	for i = 1:m
		lclosefaces[cran(firstv[2],i)]=i
	end
	if setcard == 2
		return lclosefaces[columnsinorder]'
	end
	for i = 3:(setcard-1)
		lclosefaces = buildclosefromclose(farfaces[i-1],firstv[i-1],lclosefaces,farfaces[i],firstv[i];facecard=i-1)
	end
	lclosefaces = buildclosefromclose(farfaces[setcard-1],firstv[setcard-1],lclosefaces,farfaces[setcard],firstv[setcard],columnsinorder;facecard = setcard-1)
	return lclosefaces
end

function buildallfromfar(farfaces,firstv,setcard,columnsinorder;verbose = false)
	m = length(firstv[1])-1
	n = length(farfaces[setcard])
	# destinationmatrix = Array(Int64,setcard,n)
	if setcard == 1
		return Array(Int64,0,m)
	end
	lclosefaces = Array(Int64,1,firstv[2][end]-1)
	for i = 1:m
		lclosefaces[cran(firstv[2],i)]=i
	end
	if setcard == 2
		return vcat(lclosefaces[columnsinorder]',farfaces[setcard][columnsinorder]')
	end
	for i = 3:(setcard-1)
		lclosefaces = buildclosefromclose(farfaces[i-1],firstv[i-1],lclosefaces,farfaces[i],firstv[i];facecard=i-1)
		gc()
	end
	lclosefaces = buildallfromclose(farfaces[setcard-1],firstv[setcard-1],lclosefaces,farfaces[setcard],firstv[setcard],columnsinorder;verbose = verbose)
	gc()
	return lclosefaces
end

function buildfacematfromfar(farfaces,firstv,setcard,rowsinorder,columnsinorder;verbose=false)
	if verbose
		print("PLEASE NOTE, COLUMN INDICES MUST BE IN SORTED ORDER FOR THIS TO WORK PROPERLY.
		NO MATTER HOW ROW INDICES ARE DELIVERED, THE OUTCOME WILL BE AS THOUGH THE UNIQUE LIST
		OF SORTED ROWS WITHOUT REPETITION WERE INPUT.")
	end
	
	ff = farfaces[setcard]
	nh = length(columnsinorder)	
	nl = length(farfaces[setcard-1])
	rowsupp = falses(nl)
	rowsupp[rowsinorder]=true
	lclosefaces = buildclosefromfar(farfaces,firstv,setcard,columnsinorder)
	
	preallocationspace = 0
	for k in lclosefaces
		if rowsupp[k]
			preallocationspace+=1
		end
	end
	for i = 1:nh
		if rowsupp[ff[columnsinorder[i]]]
			preallocationspace+=1
		end
	end
	
	facerv = Array(Int64,preallocationspace)
	facecp = Array(Int64,nh+1)
	facecp[1] = 1
	marker = 0
	for jp = 1:nh
		j = columnsinorder[jp]
		for k = 1:(setcard-1)
			row = lclosefaces[k,jp]
			if rowsupp[row]
				marker+=1
				facerv[marker]=row
			end
		end
		if rowsupp[ff[j]]
			marker+=1
			facerv[marker]=ff[j]
		end
		facecp[jp+1]=marker+1		
	end
	return facerv,facecp
end

function ff2aflight_sc2(farfaces,firstv,columns)	
	setcard = 2
	f0faces::Array{Int64,1} = farfaces[setcard] 
	colptr::Array{Int64,1} = firstv[2]
	columnpost::Int64   = 1
	columnpostp1::Int64 = 2
	faces::Array{Int64,2} = Array(Int64,2,length(columns))

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
	setcard = 3
	fcfaces::Array{Int64,2} = buildclosefromfar(farfaces,firstv,setcard-1,1:length(farfaces[2]))
	
	f0faces::Array{Int64,1} = farfaces[setcard]
	f1faces::Array{Int64,1} = farfaces[setcard-1]		

	fvscm0::Array{Int64,1}  = firstv[setcard]
	fvscm1::Array{Int64,1}  = firstv[setcard-1]
	fvscm2::Array{Int64,1}  = firstv[setcard-2]	
	
	holdi=[1];holdip1=[2]
	t1::Array{Int64,1} = Array(Int64,fvscm2[end]-1);t1[crows(fvscm1,f1faces,1)]=cran(fvscm1,1)		
	
	faces::Array{Int64,2} = Array(Int64,3,length(columns))
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

function ff2aflight_scgt3(farfaces,firstv,setcard,columns)							
	f0faces::Array{Int64,1} = farfaces[setcard]
	f1faces::Array{Int64,1} = farfaces[setcard-1]
	f2faces::Array{Int64,1} = farfaces[setcard-2]	
	fcfaces::Array{Int64,2} = buildallfromfar(farfaces,firstv,setcard-2,1:(firstv[setcard-2][end]-1))

	fvscm0::Array{Int64,1}  = firstv[setcard]
	fvscm1::Array{Int64,1}  = firstv[setcard-1]
	fvscm2::Array{Int64,1}  = firstv[setcard-2]
	fvscm3::Array{Int64,1}  = firstv[setcard-3]

	holdi=[1];holdip1=[2];holdj=[1];holdjp1=[2]
	t1::Array{Int64,1} = Array(Int64,fvscm2[end]-1);t1[crows(fvscm1,f1faces,1)]=cran(fvscm1,1)
	t2::Array{Int64,1} = Array(Int64,fvscm3[end]-1);t2[crows(fvscm2,f2faces,1)]=cran(fvscm2,1)

	scm0::Int64 = setcard; scm1::Int64 = setcard-1; scm2::Int64 = setcard-2
	faces::Array{Int64,2} = Array(Int64,setcard,length(columns))

	for fp = 1:length(columns)
		f0 = columns[fp]
		f1 = f0faces[f0]
		f2 = f1faces[f1]
		updatetranslator!(f0::Int64,fvscm0::Array{Int64,1},holdi::Array{Int64,1},holdip1::Array{Int64,1},t1::Array{Int64,1},fvscm1::Array{Int64,1},f1faces::Array{Int64,1})
		updatetranslator!(f1::Int64,fvscm1::Array{Int64,1},holdj::Array{Int64,1},holdjp1::Array{Int64,1},t2::Array{Int64,1},fvscm2::Array{Int64,1},f2faces::Array{Int64,1})			
		for i = 1:scm2
			fcfaces[i,f2]
			t2[fcfaces[i,f2]]			
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
			fcfaces[i,f2]
			t2[fcfaces[i,f2]]			
			faces[i,fp] = t1[t2[fcfaces[i,f2]]]
		end
		faces[scm1,fp] = t1[f2]
		faces[scm0,fp] = f1
	end			
end

function ff2aflight(farfaces,firstv,setcard,columns)	
	if setcard == 1
		return Array(Int64,0,length(columns))	
	elseif setcard == 2
		return ff2aflight_sc2(farfaces,firstv,columns)	
	elseif setcard == 3
		return ff2aflight_sc3(farfaces,firstv,columns)	
	else
		return ff2aflight_scgt3(farfaces,firstv,setcard,columns)	
	end		
end

function ff2aflight(D::Dict,setcard,columns)
	farfaces = D["farfaces"]; firstv = D["firstv"]
	faces = ff2aflight(farfaces,firstv,setcard,columns)
	return faces
end

#### end

function filteredmatrixfromfarfaces{Tv}(
	farfaces,
	firstv,
	prepairs,
	filtration,
	setcard::Integer,
	lowbasisnames::Array{Tv,1};
	verbose = false)	
	#### low filtration is to be arranged least to greatest
		
	numhigs = length(farfaces[setcard])
	numlows = length(farfaces[setcard-1])
	numppair= length(prepairs[setcard])	
			
	pphigs = prepairs[setcard]
	pplows = farfaces[setcard][pphigs]		
  	lpls = lowbasisnames
	hphs = farfaces[setcard+1][prepairs[setcard+1]]
	nplows = intervalcomplementuniqueunsortedinput(vcat(lpls,pplows),numlows)	
	nphigs = intervalcomplementuniqueunsortedinput(vcat(hphs,pphigs),numhigs)		
	
	numnhph = numhigs-length(hphs)
	Ml = numlows - length(lpls)
	Mh = numhigs - length(hphs)
	
	higtranslator = zeros(Tv,numnhph)	
	lowtranslator = zeros(Tv,numlows)		
	lowtranslator[pplows] = 1:numppair		
	
	if !isempty(nplows) && setcard > 2		
		npfilt = filtration[setcard-1][nplows]	
		nporder = integersinsameorder(npfilt,maximum(npfilt))
		addinteger!(nporder,numppair)
	else
		nporder = (numppair+1):(numppair+length(nplows))	
	end

	lowtranslator[nplows] = nporder
	higsinpointorder = intervalcomplementuniqueunsortedinput(hphs,numhigs)
	lowlab = Array(Int64,Ml)
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
	allfaces = buildallfromfar(farfaces,firstv,setcard,higsinpointorder;verbose = verbose)	
	if verbose
		print("done building allfromfar")		
	end	
	Mrv,Mcp,Mm = presparsefull2unsortedsparsetranspose(allfaces,lowtranslator,higtranslator;verbose=verbose)	
	higtranslator = [];npfilt = [];ppsupp = [];allfaces = [] 			
	gc()												
	if verbose && length(Mrv)>(Mcp[end]-1)
		print("There was the thought that Mrv should have no extra elements")
		sleep(3)
	end
	return Mrv,Mcp,lowlab,higlab,Mm,lowtranslator	
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
	prows = Array(Int64,keptmarker)
	pcols = Array(Int64,keptmarker)
	for i = 1:keptmarker
		keptindex = keptlist[i]
		prows[i] = Jprows[keptindex]
		pcols[i] = Jpcols[keptindex]
	end

	Arv,Crv,Acp,Ccp = stackedsubmatrices(Mrv,Mcp,prows,comprows,pcols,Mm0) 
	Brv,Drv,Bcp,Dcp = stackedsubmatrices(Mrv,Mcp,prows,comprows,compcols,Mm0) 
	Lrv,Lcp = copycolumnsubmatrix(Trv,Tcp,pcols)
	Rrv,Rcp = copycolumnsubmatrix(Trv,Tcp,compcols)

	translator = Array(Int64,Mm0)
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

function morsereducechaincomplex{Tv<:Integer}(
	Mrv::Array{Tv,1},
	Mrowfiltration::Array{Tv,1},
	Mcp::Array{Tv,1},
	Mcolfiltration::Array{Tv,1},
	lowlab::Array{Tv,1},
	higlab::Array{Tv,1},
	pplow,pphig,
	Mm::Integer,
	setcard,
	lowlab0;
	storetransform = true,
	verbose = false,
	diagnostic = false)
	#### Please note: the columns of M must be ordered by filtration

 	rowlab = higlab;collab = lowlab	
	if isempty(lowlab0)
		maxlowlab0 = 1
	else
		maxlowlab0 = maximum(lowlab0)
	end	
	lowtranslator0 = zeros(Int64,maxlowlab0)
	lowtranslator0[lowlab0] = 1:length(lowlab0)
	plowsupp = falses(maxlowlab0)
	if diagnostic && maxlowlab0<length(lowlab0)
		println("please check lowlab0")
		sleep(3)
	end	
	Mm = [length(higlab)]
	Mn = [length(lowlab)]
	Mn0 = Mn[1]				
	
	maxnumpairs = min(Mm[1],Mn[1]); numjunpairs = [length(pplow)]; numsenpairs = [0]
	Sprows=Array(Tv,maxnumpairs);Spcols=Array(Tv,maxnumpairs);	
	Jprows=Array(Tv,maxnumpairs);Jpcols=Array(Tv,maxnumpairs);
	Jprows[1:numjunpairs[1]]=pphig;Jpcols[1:numjunpairs[1]]=pplow
	comprows = convert(Array{Tv,1},(numjunpairs[1]+1):Mm[1])
	compcols = convert(Array{Tv,1},(numjunpairs[1]+1):Mn[1])	
	
	Trv=Array(Tv,0);Srv = Array(Tv,0)
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
				
	if diagnostic
		if numsenpairs[1]==0
			println(["numjunpairs" numjunpairs[1] "numsenpairs" numsenpairs[1]])
			println("please check numsenpairs")
		else
			println(["num zero collab" countnz(collab.==0)
				"num zero spcols" countnz(Spcols[1:numsenpairs[1]].==0)])
			println(["max Spcols[1:numsenpairs]" maximum(Spcols[1:numsenpairs[1]]) 
				"min Spcols[1:numsenpairs]" minimum(Spcols[1:numsenpairs[1]]) 
				"length(lowtranslator0)" length(lowtranslator0)])
		end
		plowsupp[lowlab0[Spcols[1:numsenpairsOLD]]]=true
		for k in Srv[1:(Scp[numsenpairs[1]+1]-1)]
			if !plowsupp[lowlab0[k]]
				println("it may be time to assess trv")
				sleep(3)
				break
			end
		end
	end
	if verbose	
		println("first shurr finished")
		if Mn[1]>0
			println([Mcp[Mn[1]] "nnz(M)" Mm[1] "Mm" Mn[1] "Mn" length(Mrv) "length(Mrv)"])
		else
			println("Mn = 0")
		end
	end
	gc()		
	rowfilt = Array(Tv,length(comprows)); colfilt = Array(Tv,length(compcols))	
	counter = 0	
	while Mcp[Mn[1]+1]>1		
		if verbose
			println("starting block operation $(counter)")
		end
		counter+=1		
		for i = 1:Mm[1]
			rowfilt[i] = Mrowfiltration[rowlab[i]]
		end
		for j = 1:Mn[1]
			colfilt[j] = Mcolfiltration[collab[j]]
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
					
		if diagnostic
			if numsenpairs[1]==0
				println(["numjunpairs" numjunpairs[1] "numsenpairs" numsenpairs[1]])
				println("please check numsenpairs")
			else
				println(["num zero collab" countnz(collab.==0)
					"num zero spcols" countnz(Spcols[1:numsenpairs[1]].==0)])
				println(["max Spcols[1:numsenpairs]" maximum(Spcols[1:numsenpairs[1]]) 
					"min Spcols[1:numsenpairs]" minimum(Spcols[1:numsenpairs[1]]) 
					"length(lowtranslator0)" length(lowtranslator0)])
			end
			plowsupp[lowlab0[Spcols[1:numsenpairsOLD]]]=true
			for k in Srv[1:(Scp[numsenpairs[1]+1]-1)]
				if !plowsupp[lowlab0[k]]
					println("it may be time to assess trv")
					sleep(3)
					break
				end
			end									
			if Mn[1]>0
				println([Mcp[Mn[1]] "nnz(M)" Mm[1] "Mm" Mn[1] "Mn" length(Mrv) "length(Mrv)"])
			end
		end
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
	if diagnostic
		if length(tlab)!=maxnumpairs
			print("please check tlab")
			sleep(3)
		elseif Mn[1]+numsenpairs[1]!=maxnumpairs
			print("please check Mn[1]")
			sleep(3)
		end
	end			
	return Srv,Scp,Sprows,Spcols,tlab
end

function persistF2_core(farfaces,firstv,prepairs,filtration,maxcard::Integer;record="cyclerep",verbose=false)	
	if record == "all" || record == "cyclerep"
		storetransform = true
	else
		storetransform = false
	end
	m = length(firstv[1])-1
	trv=Array{Any}(maxcard+1); tcp=Array{Any}(maxcard+1); phigs=Array{Any}(maxcard+1); 
	plows=Array{Any}(maxcard+1);lowtranslator=Array{Any}(maxcard+1);higtranslator=Array{Any}(maxcard+1)
	lowlab0 = Array{Any}(maxcard+1);tid = Array{Any}(maxcard+1)
	tcp[1]=fill(1,m);
	tcp[maxcard+1] = [1];
	for i in [1,maxcard+1]
		trv[i]				=Array(Int64,0)
		phigs[i]			=Array(Int64,0)
		plows[i]			=Array(Int64,0)
		lowtranslator[i]	=Array(Int64,0)
		lowlab0[i]			=Array(Int64,0)
		tid[i]				=Array(Int64,0)
	end
	initialnnz = 0	
	for setcard = 2:maxcard
		if setcard>2
			lowbasisnames = phigs[setcard-1]
		else
			lowbasisnames = Array(Int64,0)
		end	
		Mrv,Mcp,lowlab,higlab,Mm,lowtranslator[setcard] = 
		filteredmatrixfromfarfaces(farfaces,firstv,prepairs,filtration,setcard,lowbasisnames;verbose=verbose)
		lowlab0[setcard] = copy(lowlab)
 		lowlabtemp = convert(Array{Int64,1},1:length(lowlab))
 		higlabtemp = convert(Array{Int64,1},1:length(higlab))
 		higfilttemp = filtration[setcard][higlab]
 		lowfilttemp = filtration[setcard-1][lowlab]							
		if verbose
			println("Constructed Morse boundary operator, columns indexed by cells of dimension $(setcard-1)")
		end		
		pplow = convert(Array,length(prepairs[setcard]):-1:1)
		pphig = convert(Array,length(prepairs[setcard]):-1:1)
		Srv,Scp,Sphigs,Splows,tlab = 
		morsereducechaincomplex(
			Mrv,
			higfilttemp,
			Mcp,
			lowfilttemp,
			lowlabtemp,
			higlabtemp,
			pplow,
			pphig,
			Mm,
			setcard,
			lowlab0[setcard],
			storetransform = storetransform,
			verbose = verbose)
		trv[setcard] = lowlab0[setcard][Srv]
		tcp[setcard] = Scp
		tid[setcard] = lowlab[tlab]
		plows[setcard] = lowlab[Splows]
		phigs[setcard] = higlab[Sphigs]
	end
	return trv,tcp,plows,phigs,lowtranslator,lowlab0,tid
end

function persistF2!(
	D::Dict;maxcard=0,
	dictionaryoutput::Bool = true,
	verbose::Bool = false,
	record = "cyclerep")

	farfaces = D["farfaces"]
	firstv = D["firstv"]
	prepairs = D["prepairs"]
	filtration = D["filtration"]		
	if maxcard == 0
		maxcard = length(farfaces)-1
	end
	
	trv,tcp,plows,phigs,lowtranslator,lowlab0,tid = 
	persistF2_core(farfaces,firstv,prepairs,filtration,maxcard::Integer;record=record,verbose = verbose)
	if dictionaryoutput == true
		D["trv"] = trv
		D["tcp"] = tcp
		D["tid"] = tid		
		D["plows"] = plows
		D["phigs"] = phigs
		D["lowtranslator"] = lowtranslator
		D["lowlab0"] = lowlab0
		return D
	else
		return trv,tcp,plows,phigs,lowtranslator,lowlab0,tid
	end
end

function persistF2(
	s,
	maxcard; 
	rowsare = "distances",
	filetype = "textfile",
	minval=-Inf,
	maxval=Inf,
	fillud="up",
	numlevels=Inf,
	conestop=true,
	pointlabels = [],
	metric = Euclidean(),
	verbose = false,
	record = "cyclerep")	
	
	#### Start timer
	tic()
	
	#### Extract data as necessary
	inputisfile = false
	if typeof(s) == String
		inputisfile = true
		filename = copy(s)		
		if filetype == "textfile"
			s = readdlm(s)
		elseif filetype == "perseus_distmat"
			s,minval,maxval,maxcard = parse_perseusdistmat(filename)
			rowsare = "distances"
			fillud = "up"
			numlevels = Inf
		elseif filetype == "perseus_brips"
			s,minval,maxval = parse_perseusbrips(filename)			
			rowsare = "points"
			fillud = "up"
			numlevels = Inf
		end				
	end	
	
	#### Store the input
	inputdata = Dict(
		"inputarray" => copy(s),		
		"maxcard" => maxcard,
		"rowsare" => rowsare,
		"minval" => minval,
		"maxval" => maxval,
		"fillud" => fillud,
		"numlevels" => numlevels,
		"conestop" => conestop,
		"metric" => Euclidean(),
		"eireneversion" => "0.3.5"
	)		
	if inputisfile
		inputdata["filename_inputarray"] = filename
	end

	#### Determine the number of points
	if rowsare == "dimensions"
		numpoints = size(s,2)
	elseif rowsare == "points"
		numpoints = size(s,1)
	elseif rowsare == "distances"
		numpoints = min(size(s,1),size(s,2))
		if !issymmetric(s)
			print("It appears the input matrix is not symmetric.  Only symmetric distance matrices are accepted when the <rowsare> keyword argument has value \"distances\".")
			return
		end
	end
	
	#### Extract labels
	if pointlabels == "none" || pointlabels == []
		pointlabels = 1:numpoints
	elseif pointlabels in ["left","right","top","bottom"]
		(s,pointlabels) = separatelabels(s,pointlabels)
	elseif typeof(pointlabels) == String
		inputdata["filename_pointlabels"] = pointlabels
		pointlabels = ezread(pointlabels)
	end	
	if length(pointlabels) != numpoints
		warn("It appears the number of vertex labels does not match the number of vertices.")
	end
	pointlabels = ezlabel(pointlabels)
	for i = 1:length(pointlabels)
		pointlabels[i] = "$(pointlabels[i])"
	end
	inputdata["pointlabels"] = pointlabels

	#### Rowsare the matrix 					
	if rowsare == "points" || rowsare == "dimensions"
		if rowsare == "points"
			s = transpose(s)
		end			
		s = convert(Array{Float64},s)
		inputdata["pcloud"] = s		
		d = Distances.pairwise(metric,s)			
	elseif rowsare == "distances"	
		d = convert(Array{Float64,2},s)				
	end	
	inputdata["distmat"] = d	
	
	(t,filtrationtranslator) = ordercanonicalform(
		d;
		minval=minval,
		maxval=maxval,
		fillud=fillud,
		numlevels=numlevels,
		conestop=conestop,
		verbose = verbose)
	inputdata["ocf"] = t
	
	#### Build the complex 
	D = buildcomplex3(t,maxcard;verbose = verbose)
	D["filtrationtranslator"]=filtrationtranslator				

	#### Compute persistence 
	persistF2!(D;verbose = verbose,record = record)
	
	#### Store input data 
	D["inputdata"] = inputdata
	
	#### Store generators
	gc()
	if record == "all" || record == "cyclerep"
		unpack!(D)
	end
	gc()
	if record == "cyclerep"
		delete!(D,"trv")
		delete!(D,"tcp")
		delete!(D,"Lrv")
		delete!(D,"Lcp")
		delete!(D,"Rrv")
		delete!(D,"Rcp")
		delete!(D,"Lirv")
		delete!(D,"Licp")
	end
	
	#### Record time
	D["inputdata"]["computationtime"] = toc()
			
	return D
end

function parse_perseusdistmat(filename)
	A = readdlm(filename)
	s = convert(Array{Float64},A[3:end,:])
	minval = convert(Float64,A[2,1])
	maxval = minval+convert(Float64,A[2,2])*convert(Float64,A[2,3])
	maxcard = convert(Int64,A[2,4])+1
	return s,minval,maxval,maxcard
end

function parse_perseusbrips(filename)
	A = readdlm(filename)
	s = convert(Array{Float64},A[3:end,1:end-1])
	minval = 0
	maxval = convert(Array{Float64},A[2,1:3])
	maxval = prod(maxval)
	return s,minval,maxval
end

function boundarylimit_Lionly(farfaces,firstv,trv,tcp,plows,phigs,lowtranslator,lowlab0,tid,setcard;verbose=false)
	numl = length(farfaces[setcard-1])
	nump = length(phigs[setcard])
	numnlpl = numl-length(plows[setcard-1])
	Lrv = copy(trv[setcard])
	tid = tid[setcard]
	plows = plows[setcard]
	phigs = phigs[setcard]	
	lowlab0 = lowlab0[setcard]	
	lowtranslator1 = zeros(Int64,numl)
	lowtranslator1[tid] = 1:numnlpl
	yafterx!(lowtranslator1,Lrv)
	lowtranslator1 = []
	gc()		
	Lirv,Licp   = morseInverseF2orderedColsUnsortedRowsSilentInSilentOut(Lrv,tcp[setcard])	
	Lirv,Licp   = transposeLighter(Lirv,Licp,numnlpl)	
	Lrv = [];		
	gc()
	return Lirv,Licp		
end

function boundarylimit(farfaces,firstv,trv,tcp,plows,phigs,lowtranslator,lowlab0,tid,setcard;verbose=false)
	numl = length(farfaces[setcard-1])
	nump = length(phigs[setcard])
	numnlpl = numl-length(plows[setcard-1])
	trv = copy(trv[setcard])
	tcp = copy(tcp[setcard])
	tid = tid[setcard]
	plows = plows[setcard]
	phigs = phigs[setcard]	
	lowlab0 = lowlab0[setcard]	
	lowtranslator1 = zeros(Int64,numl)
	lowtranslator1[tid] = 1:numnlpl
	yafterx!(lowtranslator1,trv)
	
	Lirvt,Licpt = morseInverseF2orderedColsUnsortedRowsSilentInSilentOut(trv,tcp)	
	Lirv,Licp   = transposeLighter(Lirvt,Licpt,numnlpl)			
	trv,tcp 	= transposeLighter(trv,tcp,numnlpl)
	Lrv = trv; Lcp = tcp			
	brv = ff2aflight(farfaces,firstv,setcard,phigs)
	brv = reshape(brv,length(brv))
	bcp = convert(Array{Int64,1},1:setcard:(nump*setcard+1))
	supportedmatrix!(brv,bcp,tid,1:nump,numl)	
	yafterx!(lowtranslator1,brv)
	append!(bcp,fill(bcp[end],numnlpl-nump)) # <-- note there's no +1 b/c bcp is already 1 elt. too long
	brv,bcp = spmmF2silentLeft(trv,tcp,brv,bcp,numnlpl)	
	Rrv,Rcp = morseInverseF2orderedColsUnsortedRowsInSilentOut(brv,bcp)
	
	trv=[];lowtranslator1=[];brv=[];bcp=[];Lirvt=[];Licpt=[]
	gc()
	return Lrv,Lcp,Lirv,Licp,Rrv,Rcp		
end

function unpack!(D::Dict)
	####  This subroutine is currently specific to simplicial complexes
	farfaces = D["farfaces"];firstv = D["firstv"];trv = D["trv"];tcp=D["tcp"];
	plows=D["plows"];phigs=D["phigs"];lowtranslator=D["lowtranslator"];lowlab0=D["lowlab0"];
	filtration = D["filtration"];lowlab0 = D["lowlab0"];tid = D["tid"];
	l = length(farfaces)
	maxcard = getmaxdim(farfaces)-1
	
	Lirv = Array{Any}(l); Licp = Array{Any}(l)	
	Lrv = Array{Any}(l);  Lcp = Array{Any}(l)	
	Rrv = Array{Any}(l);  Rcp = Array{Any}(l)			
	for i = 2:(maxcard-1)
		Lrv[i],Lcp[i],Lirv[i],Licp[i],Rrv[i],Rcp[i] = 
		boundarylimit(farfaces,firstv,trv,tcp,plows,phigs,lowtranslator,lowlab0,tid,i)
	end
	
	Lirv[maxcard],Licp[maxcard] =
		boundarylimit_Lionly(farfaces,firstv,trv,tcp,plows,phigs,lowtranslator,lowlab0,tid,maxcard)
	Lrv[maxcard]=[];Lcp[maxcard]=[];Rrv[maxcard]=[];Rcp[maxcard]=[]	
	Lirv[1]=Array(Int64,0);Licp[1]=Array(Int64,0);Lrv[1]=Array(Int64,0);
	Lcp[1]=ones(Int64,1);Rrv[1]=Array(Int64,0);Rcp[1]=ones(Int64,1)
	D["Lrv"] = Lrv
	D["Lcp"] = Lcp
	D["Lirv"]= Lirv
	D["Licp"]= Licp
	D["Rrv"] = Rrv
	D["Rcp"] = Rcp	
	
	D["cyclerep"] = fill([],maxcard)
	
	for i = 2:maxcard
		dim = i-2
		m = numbars(D,dim=dim)
		cyclenames = barname2cyclename(D,1:m,dim=dim)
		D["cyclerep"][i] = getcycle(D,cyclenames,dim=dim)
	end
	return	
end

##########################################################################################

####	INVERSION

##########################################################################################

function morseInverseF2orderedColsUnsortedRowsInSilentOut{Tv<:Integer}(Arowval::Array{Tv,1},Acolptr::Array{Tv,1})
	# note q should be a vector sending a row to the column it pairs with; columns will be
	# permuted such that the resulting matrix is upper-triangular
	
	mA = length(Acolptr)-1
	const colptrA = Acolptr
	const rowvalA = Arowval
    const preallocationIncrement = colptrA[end]    
    
    colptrC = Array(Tv,mA+1); colptrC[1]=1
	rowSupp = zeros(Tv, mA)
	rowList = Array(Tv, mA)
	rowvalCj = Array(Bool, mA)		
	rowvalC = Array(Tv, mA)
    totalrowscounter = 0    
    onepast = 0		
	for i in 1:mA
		if colptrC[i]+mA > length(rowvalC)+1
			append!(rowvalC, Array(Int64,preallocationIncrement))
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
	# note q should be a vector sending a row to the column it pairs with; columns will be
	# permuted such that the resulting matrix is upper-triangular
	
	mA = length(Acolptr)-1
	const colptrA = Acolptr
	const rowvalA = Arowval
    const preallocationIncrement = colptrA[end]    
    
    colptrC = Array(Tv,mA+1); colptrC[1]=1
	rowSupp = zeros(Tv, mA)
	rowList = Array(Tv, mA)
	rowvalCj = Array(Bool, mA)		
	rowvalC = Array(Tv, mA)
    totalrowscounter = 0    
    onepast = 0 		
	for i in 1:mA
		if colptrC[i]+mA > length(rowvalC)+1
			append!(rowvalC, Array(Int64,preallocationIncrement))
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

function spmmF2{Tv<:Integer}(Arowval::Array{Tv,1},Acolptr::Array{Tv,1},Browval::Array{Tv,1},Bcolptr::Array{Tv,1},Am)
    const mA = Am
    const nB = length(Bcolptr)-1
    const rowvalA = Arowval; colptrA = Acolptr
    const rowvalB = Browval; colptrB = Bcolptr
    const preallocationIncrement = colptrA[end]+colptrB[end]
    
	colptrC = Array(Tv, nB+1)
    colptrC[1] = 1    
	rowSupp = zeros(Tv, mA)
	rowList = Array(Tv, mA)
	rowvalCj = BitArray(mA)		
	rowvalC = Array(Tv, preallocationIncrement)        	
	for i in 1:nB
		newrowscounter = [0]
		spmmF2_subr2(
			rowvalCj::BitArray{1},rowSupp::Array{Int64,1},i::Int64,
			newrowscounter::Array{Int64,1},colptrA::Array{Int64,1},rowvalA::Array{Int64,1},
			rowList::Array{Int64,1},colptrB::Array{Int64,1},rowvalB::Array{Int64,1},i::Int64)		
		nzRows = find(rowvalCj[rowList[1:newrowscounter[1]]])

		colptrC[i+1] = colptrC[i]+length(nzRows)
		if colptrC[i+1] > length(rowvalC)+1
			append!(rowvalC, Array(Int,preallocationIncrement))
		end
				
		rowvalC[colptrC[i]:(colptrC[i+1]-1)] = sort(rowList[nzRows])	
	end
	deleteat!(rowvalC,colptrC[end]:length(rowvalC))		
	return rowvalC, colptrC
end

function spmmF2_subr2(
	oddfloods::BitArray{1},shoreline::Array{Int64,1},watermark::Int64,
	peakcounter::Array{Int64,1},Acp::Array{Int64,1},Arv::Array{Int64,1},
	flippedlist::Array{Int64,1},Bcp::Array{Int64,1},Brv::Array{Int64,1},col)
	for jp in Bcp[col]:(Bcp[col+1] - 1)
		j = Brv[jp]				
		addcol!(
			oddfloods::BitArray{1},shoreline::Array{Int64,1},watermark::Int64,
			peakcounter::Array{Int64,1},Acp::Array{Int64,1},Arv::Array{Int64,1},
			flippedlist::Array{Int64,1},j::Int64)					
	end	
end

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

function spmmF2silentLeft{Tv<:Integer}(Arowval::Array{Tv,1},Acolptr::Array{Tv,1},Browval::Array{Tv,1},Bcolptr::Array{Tv,1},Am)
    const mA = Am
    const nB = length(Bcolptr)-1
    const rowvalA = Arowval; colptrA = Acolptr
    const rowvalB = Browval; colptrB = Bcolptr
    const preallocationIncrement = colptrA[end]+colptrB[end]
    
	colptrC = Array(Tv, nB+1)
    colptrC[1] = 1    
	rowSupp = zeros(Tv, mA)
	rowList = Array(Tv, mA)
	rowvalCj = Array(Bool, mA)		
	rowvalC = Array(Tv, preallocationIncrement)        
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
			append!(rowvalC, Array(Int,preallocationIncrement))
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
	rowList = Array(Int64, Dm)
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
	rowList = Array(Int64, Dm)
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

function cran(colptr::Array,j)
	return colptr[j]:(colptr[j+1]-1)
end

function cran(colptr::UnitRange,j)
	return colptr[j]
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
		append!(x,Array(Tv,n-length(x)))
	end
end

function copycolumnsubmatrix{Tv<:Integer}(Arv::Array{Tv,1},Acp::Array{Tv,1},columnindices::Array{Tv,1})
	allocationspace = 0
	for j in columnindices
		allocationspace+= Acp[j+1]-Acp[j]
	end	
	Brv = Array(Tv,allocationspace)
	Bcp = Array(Tv,length(columnindices)+1)	
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
	Brv = Array(Tv,allocationspace)
	Bcp = Array(Tv,length(columnindices)+1)	
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
		append!(rowvalB,Array(Tv,numnewrows+growthIncrement))
	end
	if length(colptrB)<startingDestination+length(columnindices)
		append!(colptrB,Array(Tv,startingDestination+length(columnindices)))
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
		append!(r,Array(Int64,max(growthincrement,length(v))))
	end
	r[startpoint:(c[k+1]-1)]=v
end

function supportedmatrix!{Tv<:Integer}(Mrowval::Array{Tv},Mcolptr::Array{Tv,1},rows1,colsinorder,Mm::Tv)    
	### note cols must be in order
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
	Mrowval::Array{Tv,1},
	Mcolptr::Array{Tv,1},
	rows1::Array{Tv,1},
	rows2::Array{Tv,1},
	cols::Array{Tv,1},
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
	rv1 = Array(Tv,nz1)
	rv2 = Array(Tv,nz2)
	cp1 = Array(Tv,n+1); cp1[1]=1
	cp2 = Array(Tv,n+1); cp2[1]=1
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
    Ccolptr = Array(Ti,Am+1)
    Crowval = Array(Ti,Annz)
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
    for Aj in 1:An
        for Ak in Acolptr[Aj]:(Acolptr[Aj+1]-1)
            Ai = Arowval[Ak]
            Ck = Ccolptr[Ai+1]
            Crowval[Ck] = Aj
            Ccolptr[Ai+1] += 1
        end
    end
    # Tracking write positions in Ccolptr as in the last block fixes the colptr shift,
    # but the first colptr remains incorrect
    Ccolptr[1] = 1

	return Crowval, Ccolptr
end

function transposeLighter{Tv,Ti}(Arowval::Array{Ti},Acolptr::Array{Ti},Anzval::Array{Tv},Am::Integer)
    Annz = Acolptr[end]-1
    An = length(Acolptr)-1
    Cm = An
    Cn = Am
    Ccolptr = Array(Ti,Am+1)
    Crowval = Array(Ti,Annz)
    Cnzval = Array(Tv,Annz)
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

function transposeLighter_submatrix{Ti}(Arowval::Array{Ti},Acolptr::Array{Ti},Am;rows = 1:Am,cols = 1:length(Acolptr)-1)	
	if rows == 1:Am && cols == 1:(length(Acolptr)-1)
		Crowval, Ccolptr = transposeLighter(Arowval::Array{Ti},Acolptr::Array{Ti},Am)
		return Crowval, Ccolptr
	end	
    # Attach destination matrix
    Cm = length(cols)
    Cn = length(rows)
    Ccolptr = Array(Ti,Cn+1)
    # Compute the column counts of C and store them shifted forward by one in Ccolptr
    Ccolptr[1:end] = 0    
	rs = rowsupportsum(Arowval,Acolptr,Am,cols)	
	for i = 1:Cn
	    Ccolptr[i+1] = rs[rows[i]]
	end
	Cnnz = sum(Ccolptr)
    Crowval = Array(Ti,Cnnz)   
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

function presparsefull2unsortedsparsetranspose{Tv<:Integer}(
	M::Array{Tv,2},
	row02row1translator,
	col02col1translator;
	verbose::Bool=false)
	Mm,Mn = size(M)
	if Mn == 0
		rowval1 = Array(Int64,0)
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
	colptr1 = Array(Int64,m0+1)
	colptr1[1]=1
	for i = 1:m0
		colptr1[i+1]=colptr1[i]+rowcounter[i]
	end
	rowval1 = Array(Int64,colptr1[end]-1)
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
	gc()
	return rowval1,colptr1,Mn
end

##########################################################################################

####	SHAPE GENERATORS

##########################################################################################

function noisycircle()
	theta = 1:100
	theta = theta*2*pi/100
	x = cos(theta)
	y = sin(theta)
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
	
	x = Array(Float64,0)
	y = Array(Float64,0)
	z = Array(Float64,0)
	
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

function latlon2euc(A;rowsare = "dimensions")
	if rowsare == "dimensions"
		theta1 = A[1,:]
		theta2 = A[2,:]
	elseif rowsare == "points"
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
	C = eirene(d[:,supp];rowsare = "dimensions",bettimax=1,upperlim = 0.25,record="cyclerep")
	return a,d,supp,C
end

##########################################################################################

####	MATRIX WEIGHTS AND FORMATTING

##########################################################################################

function ordercanonicalform{Tv}(
	S::Array{Tv};
	minval=-Inf,
	maxval=Inf,
	fillud="up",
	numlevels=Inf,
	conestop::Bool=true,
	verbose::Bool=false)
       
    symmat_float = convert(Array{Float64},copy(S))
    symmat = copy(S)
	m = size(symmat,1);
	convert(Tv,minval)
	convert(Tv,maxval)	

	if fillud=="up"
		effectivemin = -maxval
		effectivemax = -minval		
		symmat = -symmat;
	else
		effectivemin = minval
		effectivemax = maxval
	end	
	for i = 1:m
		symmat[i,i]=Inf
	end		
	if conestop
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
	if numlevels == 1
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
	end
	if numlevels == Inf
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
					# filtration to appear to
					# live forever
					symmat[i,j] = -Inf
					symmat[j,i] = -Inf					
				end
			end
		end	
		for i = 1:m
			symmat[i,i] = -Inf
		end			
		filtrationtranslator = zeros(Float64,binom(m,2))
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
				filtrationtranslator[ordervalue]=symmat_float[ii]				
				symmat[ii]=ordervalue
			end
		end
		deleteat!(filtrationtranslator,(ordervalue+1):binom(m,2))			
	end
	return round(Int64,symmat),filtrationtranslator   
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
	l = Array(Int64,m)
	lDown = Array(Int64,m)
	val = Array(Any,m)
	supp = Array(Any,m)
	suppDown = Array(Any,m)	
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
		return Array(Int64,0)
	else	
		boundMarker = 1
		upperBound = v[boundMarker]
		complement = Array(Int64,n-L)
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
		return Array(Int64,0)
	else			
		complementsupport = trues(n)
		complementsupport[v]=false
		complement = Array(Int64,n-L)
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

function integersinsameorder!(v::Array{Int64,1},maxvalue::Int64)
	m = length(v)
	x = zeros(Int64,maxvalue)
	for i = 1:m
		x[v[i]]+=1
	end
	y = Array(Int64,maxvalue+1)
	y[1] = 1
	for i = 1:maxvalue	
		y[i+1]=y[i]+x[i]
	end
	for i = 1:length(v)
		u = v[i]
		v[i] = y[u]
		y[u]+=1
	end
	return v
end

function integersinsameorder(v::Array{Int64,1},maxvalue::Int64)
	if isempty(v)
		z = Array(Int64,0)
		return z
	else
		m = length(v)
		x = zeros(Int64,maxvalue)
		y = Array(Int64,maxvalue+1)
		z = Array(Int64,length(v))	
		for i = 1:m
			x[v[i]]+=1
		end
		y[1] = 1
		for i = 1:maxvalue	
			y[i+1]=y[i]+x[i]
		end
		for i = 1:length(v)
			u = v[i]
			z[i] = y[u]
			y[u]+=1
		end
		return z
	end
end

function integersinsameorderbycolumn(v::Array{Int64,1},maxvalue::Int64,colptr)
	numcols = length(colptr)-1
	m = length(v)
	x = Array(Int64,maxvalue)
	y = Array(Int64,maxvalue+1)
	z = Array(Int64,length(v))	
	for j = 1:numcols
		x[:] = 0
		for i = colptr[j]:(colptr[j+1]-1)
			x[v[i]]+=1
		end
		y[1] = colptr[j]
		for i = 1:maxvalue	
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

function integersortperm(v::Array{Int64,1},maxvalue::Int64)
	l = length(v)
	u = integersinsameorder(v,maxvalue)
	w = Array(Int64,l)
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
	maximumcolsum = 0
	for i = 1:n
		colwisesum[i]+=1  # we do this so our column-sorting trick doesn't try to access index 0
		if colwisesum[i]>maximumcolsum
			maximumcolsum = colwisesum[i]
		end
	end
	colwisesumlinearized = integersinsameorderbycolumn(colwisesum,maximumcolsum,colfiltptr)
	colnamesinorder = Array(Tv,n)
	colnamesinorder[colwisesumlinearized]=1:n	
	ncoveredsupp = trues(m)
	pairmarker = 0	
	for jp = 1:n
		j = colnamesinorder[jp]
		if col2firstplace[j]>0			
			firstplace = col2firstplace[j]			
			firstrow = rowval[firstplace]
			if firstplace==colptr[j+1]-1
				if ncoveredsupp[firstrow]			
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
				if ncoveredsupp[firstrow]
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
		return Array(Int64,0)
	end
	n = length(prows)
	rowtranslator = Array(Int64,Mm)
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
	downstreamelements = Array(Int64,counter)
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


function numbars(D::Dict;dim = 1)
	setcard = dim+2
	plows = D["plows"][setcard]
	phigs = D["phigs"][setcard]
	lowfilt = D["filtration"][setcard-1]
	higfilt = D["filtration"][setcard]
	counter::Int64 = 0
	for i = 1:length(plows)
		if higfilt[phigs[i]]!=lowfilt[plows[i]]
			counter+=1
		end
	end
	counter += length(D["tid"][setcard])-length(plows)
	return counter
end

function barname2cyclename(D::Dict,barnumber = [1];dim = 1)
	if typeof(barnumber) <: Array
		setcard = dim+2
		tid = D["tid"][setcard]	
		plows = D["plows"][setcard]
		phigs = D["phigs"][setcard]	
		nummortals = length(plows)
		nzcycles = find(D["filtration"][setcard][phigs].!= D["filtration"][setcard-1][plows])
		append!(nzcycles,Array((nummortals+1):length(tid)))
		return nzcycles[barnumber]		
	elseif typeof(barnumber) <: UnitRange
		setcard = dim+2
		tid = D["tid"][setcard]	
		plows = D["plows"][setcard]
		phigs = D["phigs"][setcard]	
		nummortals = length(plows)
		nzcycles = find(D["filtration"][setcard][phigs].!= D["filtration"][setcard-1][plows])
		append!(nzcycles,nummortals+1:length(tid))	
		return nzcycles[barnumber]
	elseif typeof(barnumber)<:Number
		setcard = dim+2
		tid = D["tid"][setcard]	
		plows = D["plows"][setcard]
		phigs = D["phigs"][setcard]	
		numclasses = length(tid)
		nummortals = length(plows)
		counter = 0
		cyclename = 0
		for i = 1:nummortals
			if D["filtration"][setcard][phigs[i]] != D["filtration"][setcard-1][plows[i]]
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

function getbetticurve(D::Dict,setcard;ocf = false)	
	if length(D["farfaces"][setcard])==0
		return Array(Float64,0,2)
	end

	maxval = Int(maximum(D["filtration"][2]))
	v = zeros(Int64,maxval)
	
	bco = barcode(D;dim = setcard-2,ocf=true)
	bco[bco.==Inf] = maxval
	bco = convert(Array{Int64},bco)	
	
	for i = 1:size(bco,1)
		v[bco[i,1]:bco[i,2]]+=1
	end
	
	if ocf == false
		u = sort(D["filtrationtranslator"])
	else
		u = Array(1:maxval)
	end	
	return hcat(u,v)
end

function getcycle(farfaces,firstv,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plows,phigs,tid,setcard,cyclenumber::Number)	
	if setcard == 2
		numlowlows = 0
	else
		numlowlows = length(farfaces[setcard-2])
	end
	numnlpl = length(farfaces[setcard-1])-length(plows[setcard-1])	
	
	summands = tid[setcard][crows(Licp[setcard],Lirv[setcard],cyclenumber)]
	append!(summands,[tid[setcard][cyclenumber]])	
	brv = ff2aflight(farfaces,firstv,setcard-1,summands)
	supp = falses(length(farfaces[setcard-2]))
	for k in brv
		supp[k] = !supp[k]
	end	
	
	brv = find(supp[tid[setcard-1]])
	bcp = [1,length(brv)+1]
	brv,bcp = spmmF2silentLeft(Lrv[setcard-1],Lcp[setcard-1],brv,bcp,numnlpl)
	brv,bcp = spmmF2silentLeft(Rrv[setcard-1],Rcp[setcard-1],brv,bcp,numnlpl)		
	plow2phigtranslator = Array(Int64,numlowlows)
	plow2phigtranslator[plows[setcard-1]]=phigs[setcard-1]
	
	brv = append!(plow2phigtranslator[tid[setcard-1][brv]],summands)	
	return brv	
end

function getcycle(farfaces,firstv,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plows,phigs,tid,setcard,cyclenumber)
	if setcard == 2
		numlowlows = 0
	else
		numlowlows = length(farfaces[setcard-2])
	end
	numlows    = length(farfaces[setcard-1])	
	numnlpl = length(farfaces[setcard-1])-length(plows[setcard-1])	
	
	numclasses = length(cyclenumber)
	summands = Array(Any,numclasses)
	rep 	 = Array(Any,numclasses)
	summandsupp = falses(numlows)
	for i = 1:numclasses
		summands[i] = tid[setcard][crows(Licp[setcard],Lirv[setcard],cyclenumber[i])]
		append!(summands[i],[tid[setcard][cyclenumber[i]]])	
		summandsupp[summands[i]]=true
	end
	
	lowgenerators = find(summandsupp)
	numlowgenerators = length(lowgenerators)
	translator = zeros(Int64,numlows)
	translator[lowgenerators] = 1:length(lowgenerators)
	
	lowfacemat = ff2aflight(farfaces,firstv,setcard-1,lowgenerators)

	supp = falses(numlowlows)
	m = size(lowfacemat,1)
	plow2phigtranslator = Array(Int64,numlowlows)
	plow2phigtranslator[plows[setcard-1]]=phigs[setcard-1]	
	for i = 1:numclasses

		supp[:] = false
		for j = 1:length(summands[i])
			for k = 1:m
				kk = lowfacemat[k,translator[summands[i][j]]]
				supp[kk] = !supp[kk]
			end
		end
	
		brv = find(supp[tid[setcard-1]])
		bcp = [1,length(brv)+1]	
		brv,bcp = spmmF2silentLeft(Lrv[setcard-1],Lcp[setcard-1],brv,bcp,numlowlows)
		brv,bcp = spmmF2silentLeft(Rrv[setcard-1],Rcp[setcard-1],brv,bcp,numlowlows)		
	
		rep[i] = append!(plow2phigtranslator[tid[setcard-1][brv]],summands[i])
	end	
	return rep	
end

function getcycle(D::Dict,setcard,cyclenumber)	
	if !haskey(D,"Lirv")
		if !haskey(D,"cyclerep")
			println("This object does not store a complete cycle basis.")
		else
			println("This object does not store a complete cycle basis, only those cycles that represent persistent homology classes.")
		end
		return		
	end
	farfaces = D["farfaces"];firstv = D["firstv"];Lirv = D["Lirv"];Licp=D["Licp"];Lrv=D["Lrv"]
	Lcp=D["Lcp"];Rrv=D["Rrv"];Rcp=D["Rcp"];plows=D["plows"];phigs=D["phigs"];tid=D["tid"]
	
	rrv = getcycle(farfaces,firstv,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plows,phigs,tid,setcard,cyclenumber)
	return rrv
end

function getcycle(D::Dict,cyclenumber;dim = 1)	
	if !haskey(D,"Lirv")
		if !haskey(D,"cyclerep")
			println("This object does not store a complete cycle basis.")
		else
			println("This object does not store a complete cycle basis, only those cycles that represent persistent homology classes.")
		end
		return		
	end
	setcard = dim+2	
	farfaces = D["farfaces"];firstv = D["firstv"];Lirv = D["Lirv"];Licp=D["Licp"];Lrv=D["Lrv"]
	Lcp=D["Lcp"];Rrv=D["Rrv"];Rcp=D["Rcp"];plows=D["plows"];phigs=D["phigs"];tid=D["tid"]
	
	rrv = getcycle(farfaces,firstv,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plows,phigs,tid,setcard,cyclenumber)
	return rrv
end

function getcyclesize(farfaces,firstv,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plows,phigs,tid,setcard,cyclenumber)	
	if setcard == 2
		numlowlows = 0
	else
		numlowlows = length(farfaces[setcard-2])
	end
	numlows    = length(farfaces[setcard-1])	
	numnlpl = length(farfaces[setcard-1])-length(plows[setcard-1])	
	
	numclasses = length(cyclenumber)
	summands = Array(Any,numclasses)
	rep 	 = Array(Int64,numclasses)
	summandsupp = falses(numlows)
	for i = 1:numclasses
		summands[i] = tid[setcard][crows(Licp[setcard],Lirv[setcard],cyclenumber[i])]
		append!(summands[i],[tid[setcard][cyclenumber[i]]])	
		summandsupp[summands[i]]=true
	end
	
	lowgenerators = find(summandsupp)
	numlowgenerators = length(lowgenerators)
	translator = zeros(Int64,numlows)
	translator[lowgenerators] = 1:length(lowgenerators)
	
	lowfacemat = ff2aflight(farfaces,firstv,setcard-1,lowgenerators)

	supp = falses(numlowlows)
	m = size(lowfacemat,1)
	plow2phigtranslator = Array(Int64,numlowlows)
	plow2phigtranslator[plows[setcard-1]]=phigs[setcard-1]	
	for i = 1:numclasses

		supp[:] = false
		for j = 1:length(summands[i])
			for k = 1:m
				kk = lowfacemat[k,translator[summands[i][j]]]
				supp[kk] = !supp[kk]
			end
		end
	
		brv = find(supp[tid[setcard-1]])
		bcp = [1,length(brv)+1]
		brv,bcp = spmmF2silentLeft(Lrv[setcard-1],Lcp[setcard-1],brv,bcp,numnlpl)
		brv,bcp = spmmF2silentLeft(Rrv[setcard-1],Rcp[setcard-1],brv,bcp,numnlpl)		
			
		rep[i] = length(brv)+length(summands[i])
	end	
	return rep	
end

function getcyclesize(D::Dict,cyclenumber;dim = 1)	
	setcard = dim+2	
	if !haskey(D,"Lirv")
		if !haskey(D,"cyclerep")
			println("This object does not store a complete cycle basis.")
		else
			println("This object does not store a complete cycle basis, only those cycles that represent persistent homology classes.")
		end
		return		
	end
	farfaces = D["farfaces"];firstv = D["firstv"];Lirv = D["Lirv"];Licp=D["Licp"];Lrv=D["Lrv"]
	Lcp=D["Lcp"];Rrv=D["Rrv"];Rcp=D["Rcp"];plows=D["plows"];phigs=D["phigs"];tid=D["tid"]
	
	rrv = getcyclesize(farfaces,firstv,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plows,phigs,tid,setcard,cyclenumber)
	return rrv
end

function getrepsize(D::Dict,classnumber;dim=1)
	setcard = dim+2
	if !haskey(D,"cyclerep")
		println("This object does not contain data about cycle representatives.")
		return
	elseif typeof(classnumber)<: Number
		return length(D["cyclerep"][dim+2][classnumber])
	else
		l = length(classnumber)
		rsize = Array(Int64,l)
		for i = 1:l
			rsize[i] = length(D["cyclerep"][dim+2][classnumber[i]])
		end
		return rsize
	end
end

##########################################################################################

####	USER-FRIENDLY BARCODE UTILITIES

##########################################################################################

function classrep(
	D::Dict;
	dim = 1,
	class = 1,
	outputtype = "faces")
	
	if outputtype == "faces"
		return classrep_faces(D,dim = dim,class = class)
	elseif outputtype == "vertices"
		return(classrep_vertices(D,dim = dim, class = class))
	end
end	

function classrep_faces(
	D::Dict;
	dim = 1,
	class = 1)
	
	setcard = dim+2
	
	if haskey(D,"cyclerep")
		rep = D["cyclerep"][setcard][class]
	else
		cyclename = barname2cyclename(D,class;dim=dim)
		rep = getcycle(D,setcard,class)		
	end
		
	vrealization = vertexrealization(D::Dict,setcard-1,rep)
	vrealization = D["oldvertsinnewplaces"][vrealization]
	return vrealization
end

function classrep_vertices(
	D::Dict;
	dim = 1,
	class = 1)
	
	setcard = dim+2
	
	if haskey(D,"cyclerep")
		rep = D["cyclerep"][setcard][class]
	else
		cyclename = barname2cyclename(D,class;dim=dim)
		rep = getcycle(D,setcard,class)		
	end
		
	vertices = incidentverts(D::Dict,setcard-1,rep)
	vertices = D["oldvertsinnewplaces"][vertices]
	return vertices
end

function cyclevertices(
	D::Dict;
	dim = 1,
	cycle = 1)
	
	setcard = dim+2
	
	rep = getcycle(D,setcard,cycle)		
	vertices = incidentverts(D::Dict,setcard-1,rep)
	vertices = D["oldvertsinnewplaces"][vertices]
	return vertices
end

function barcode(D::Dict;dim = 1,ocf = false,verbose = false, givenztidindices = false)
	setcard = dim+2
	plows = D["plows"][setcard]
	phigs = D["phigs"][setcard]	
	higfilt = D["filtration"][setcard][phigs]
	lowfilt = D["filtration"][setcard-1][plows]	
	nump = length(plows)
	numnzmortalbars = countnz(higfilt.!=lowfilt)
		
	tid = D["tid"][setcard]	
	numbars = length(tid)
	numnzevergreenbars = numbars-nump
	numnzbars = numnzmortalbars + numnzevergreenbars
	numlows = length(D["farfaces"][setcard-1])
	
	translate2plowindex			= zeros(Int64,numlows)
	translate2plowindex[plows]  = 1:nump
	nzbarcounter = 0
	deathtimes = fill(Inf64,numbars)	
	for i = 1:numbars
		k = translate2plowindex[tid[i]]
		if k>0
			deathtimes[i] = D["filtration"][setcard][phigs[k]]
		end
		if deathtimes[i] != D["filtration"][setcard-1][tid[i]]
			nzbarcounter+=1
		end
	end
	
	if verbose
		println(["typeof(deathtimes)" typeof(deathtimes)])
		println(["length(tid)" length(tid)
		"length(D[plows])" length(plows)
		"numinf(deathtimes)" countnz(deathtimes.==Inf)])
	end	
	summary = Array(Any,nzbarcounter,2)
	nzbarcounter = 0
	for i = 1:numbars
		if deathtimes[i] != D["filtration"][setcard-1][tid[i]]
			nzbarcounter+=1
			summary[nzbarcounter,1] = D["filtration"][setcard-1][tid[i]]
			summary[nzbarcounter,2] = deathtimes[i]
		end
	end		
	if ocf == false
		summary[:,1]=D["filtrationtranslator"][convert(Array{Int64},summary[:,1])]
		summary[1:numnzmortalbars,2] = D["filtrationtranslator"][convert(Array{Int64},summary[1:numnzmortalbars,2])]
	else
		summary = 1+maximum(D["filtration"][2])-summary
		summary = convert(Array{Float64},summary)
		summary[summary.==-Inf] = Inf
	end
	if givenztidindices
		tidindices = vcat(find(higfilt.!=lowfilt),Array(length(plows)+1:length(tid)))
		return (summary,tidindices)
	else
		return summary	
	end
end

function getpersistencediagramprimitives(
	C;
	dim = 1,
	ocf = false,
	flipaxis = false,
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
	setcard = dim+2
	
	if numbrs == 0
		x0=[];y0=[];l0=[];x1=[];y1=[];l1=[];x2=[];y2=[]
		return x0,y0,l0,x1,y1,l1,x2,y2
	end
	
	if showsize
		barsizes = Array(Int64,numbrs)
		for i = 1:numbrs
			barsizes[i] = length(C["cyclerep"][setcard][i])
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
	if C["inputdata"]["fillud"] == "down"
		flipaxis = !flipaxis
	end
	if !flipaxis
		x0 = rows[infrows,1]
		y0 = x0
		l0 = labels[infrows]
		x1 = rows[finrows,1]
		y1 = rows[finrows,2]
		l1 = labels[finrows]
		x2 = [minimum(rows[:,1]),maximum(rows[:,1])]
		y2 = [topheight,topheight]
	else
		y0 = rows[infrows,1]
		x0 = x0
		l0 = labels[infrows]
		y1 = rows[finrows,1]
		x1 = rows[finrows,2]
		l1 = labels[finrows]
		y2 = [minimum(rows[:,1]),maximum(rows[:,1])]
		x2 = [topheight,topheight]
	end
	return x0,y0,l0,x1,y1,l1,x2,y2
end

function plotpersistencediagram_pjs(C;dim = 1,showlabels = false,ocf = false,flipaxis = false)
	if showlabels == true
		mmode = "markers+text"
	else
		mmode = "markers"
	end
	if C["inputdata"]["fillud"] == "down"
		flipaxis = !flipaxis
	end
	
	x0,y0,l0,x1,y1,l1,x2,y2 = getpersistencediagramprimitives(
		C;
		dim = dim,
		ocf = ocf,
		flipaxis = flipaxis)
	
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
	rowsare= D["inputdata"]["rowsare"],		
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
	setcard 	= dim+2
	facecard 	= dim+1

	##
	if showcloud == []
		showcloud = true
	end	
	if embeddingobj == []
		embeddingobj = "distmat"
	end

	###
	cyclename = barname2cyclename(D,class;dim = dim)
	if haskey(D,"cyclerep")
		rep = D["cyclerep"][setcard][class]
	else
		rep = getcycle(D,setcard,cyclename)		
	end
	classvinnewspace = incidentverts(D,setcard-1,rep)
	classvinoldspace = D["oldvertsinnewplaces"][classvinnewspace]
	if showcloud == true
		vsupp = trues(length(D["farfaces"][1]))
		vsupp[classvinnewspace] = false
		compvinnewspace = find(vsupp)
		compvinoldspace = D["oldvertsinnewplaces"][compvinnewspace]
	else
		compvinoldspace = []
	end		
	
	###
	if haskey(D["inputdata"],"pointlabels")
		textlabels = D["inputdata"]["pointlabels"]
	else
		m = length(D["farfaces"][1])
		textlabels = Array(Any,m)
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
 		vrealization = D["oldvertsinnewplaces"][vrealization]
		vertexinverter = Array(Int64,maximum(classvinoldspace))
		vertexinverter[classvinoldspace]=1:length(classvinoldspace) 		
		classedges = d1faces(vrealization)
		edges_orderverts = vertexinverter[classedges]
		L = graphlaplacian_normal(edges_orderverts)		
		efact_class = eigfact(L,1:4)
	end		

	if showcloud && embeddingobj == "hop"
		hopedges = find(hoprange[1]<=D["filtration"][2] & D["filtration"][2].<=hoprange[2])
		cloudedges = vertexrealization(D,dim=1,hopedges)
		cloudedges_orderverts = vetexinverter[cloudedges]
	end
		
	###
	if coords == []
		if !haskey(D["inputdata"],"pcloud")
			print("No point cloud is available.  Please consider using the mds keyword argument to generate a Euclidean embedding from the distance matrix (see documentation).")
			return "nopointcloud","nopointcloud"
		else
			coords = D["inputdata"]["pcloud"]
			if rowsare == "dimensions" && size(coords,1)>3
				print("The input coordinates have dimension greater than 3.  The generated plot will use the first three to represent each point. For more options re: high-dimensional representations, please see the documentation for multidimensional scaling.")
				coords = coords[1:3,:]
			elseif rowsare == "points" && size(coords,2)>3
				print("The input coordinates have dimension greater than 3.  The generated plot will use the first three to represent each point. For more options re: high-dimensional representations, please see the documentation for multidimensional scaling.")
				coords = coords[:,1:3]'
				rowsare = "dimensions"
			end
			if !showcloud
				coords = coords[:,classvinoldspace]
			end
		end
	elseif coords == "mds"	
		if showcloud
			if embeddingobj == "distmat"
				metricmatrix = D["inputdata"]["distmat"]
				metricmatrix = metricmatrix - minimum(metricmatrix)
				for i = 1:size(metricmatrix,1)
					metricmatrix[i,i]=0
				end
			elseif embeddingobj == "hop"
				metricmatrix = hopdistance(cloudedges_orderverts,inputis = "edges")		
			end	
		else
			if embeddingobj == "distmat"
				metricmatrix = D["inputdata"]["distmat"][classvinoldspace,classvinoldspace]
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
		rowsare = "dimensions"
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
		rowsare = rowsare,
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
			rowsare = rowsare,
			subset = compvinoldspace,
			color = "rgb(31,119,180)",
			opacity = 0.5,
			textlabels = textlabels,
			showlabels = showlabtemp)	
		append!(data,[T2])
	end	
	
	if rowsare == "dimensions"
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
	rowsare= "dimensions",
	embeddingdim = 3,
	embeddingobj = [],
	specrange = [-Inf,Inf],
	classcolor = "spectral",
	cloudcolor = [],
	textlabels = [],
	showlabels = "cycle")
	
	if !haskey(D["inputdata"],"pcloud") && coords == []
		print("No point cloud is available.  Coordinates may be supplied by the user with the <coords> keyword argument, or generated automatically via the mds keyword argument.  Please see documentation.")
		return	"nopointcloud","nopointcloud"
	end

	data,layout = classrep_pjs(	
		D;
		dim = dim,
		class=class,
		showcloud = showcloud,
		coords = coords,	
		rowsare="dimensions",		
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
	setcard = dim+2
	return getbetticurve(D,setcard,ocf = ocf)
end

function plotbetticurve_pjs(D::Dict;dim=1,ocf = false)
	bcu = betticurve(D;dim = dim, ocf = ocf)
	T = PlotlyJS.scatter(x = bcu[:,1],y=bcu[:,2],mode = "line",line_width=1)
	L = makelayout_pjs(T)
	PlotlyJS.plot(T,L)
end

function computationtime(D::Dict)
	return D["inputdata"]["computationtime"]
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
	rowsare = "dimensions",
	subset = 1:size(pcloud,2),
	color = "rgb[250,250,250]",
	colorscale = "Jet",
	threshold = Inf,
	opacity = 1,
	outlineonly = false,
	textlabels = [],
	markersize = 5,
	showlabels = false)
	
	if rowsare == "dimensions"
		dim = size(pcloud,1)
		x = pcloud[1,subset]
		y = pcloud[2,subset]
		if dim >= 3
			z = pcloud[3,subset]
			if dim > 3
				print("It appears the dimension of the input point cloud exceeds 3. Using the first three coordinates only.\n")
			end
		end
	elseif rowsare == "points"
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
		T["text"] = Array(Any,length(x))
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
			if rowsare == "points"
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
	
	if data[1]["type"] == "scatter3d"
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
	rowsare = "dimensions",
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
		rowsare = rowsare,
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

function edgetrace_pjs(coordinates,edges;rowsare="dimensions")
	if rowsare != "dimensions"
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
	translator = Array(Int64,M)
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
	edges = Array(Int64,2,numedges)
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

function ezread(s)
	if s[end-2:end] == "csv"
		return readdlm(s,',','\r')
	elseif s[end-2:end] == "txt"
		return readdlm(s)
	else
		println("Please ensure the input file is either comma separated (.csv) or
		space delimited (.prn)")
	end	
end

function ezlabel(y)
	l = length(y)
	x = Array(Any,l)
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
	labels = Array(Any,numlabels)
	
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

function yafterx!{Tv<:Integer}(y::Array{Tv,1},x::Array{Tv,1})
	for i = 1:length(x)
		x[i] = y[x[i]]
	end
end

function yafterx{Tv}(y::AbstractVector{Tv},x)
	z = Array(Tv,length(x))
	for i = 1:length(x)
		z[i] = y[x[i]]
	end
	return z
end

function showfull(rowval,colptr,m)
	n = length(colptr)-1
	M = zeros(Int8,m,n)
	for j = 1:n
		M[rowval[cran(colptr,j)],j]=1
	end
	return M
end

function showfull(rowval,colptr,m,n)
	M = zeros(Int8,m,n)
	for j = 1:n
		M[rowval[cran(colptr,j)],j]=1
	end
	return M
end

function colsupportsum(colptr,n::Integer)
	x = Array(Int64,n)
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
	colptr = Array(Int64,length(v)+1)
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
		rv = Array(Int64,0)
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

	w = getstarweights(symmat)
	vperm = sortperm(-w)
	symmat = symmat[vperm,vperm]		

	filtration = Array{Any}(maxcard+1)
	farfaces = Array{Any}(maxcard+1)
	prepairs = Array{Any}(maxcard+1)
	firstv = Array{Any}(maxcard+1)

	farfaces[maxcard+1] = Array(Int64,0)
	firstv[maxcard+1] = ones(Int64,1)
	filtration[maxcard+1] = Array(Int64,0)
	prepairs[maxcard+1] = Array(Int64,0)

	farfaces[1] = convert(Array,1:m)
	firstv[1] = convert(Array,1:(m+1))
	filtration[1] = maximum(symmat)*ones(Int64,m)
	prepairs[1] = Array(Int64,0)

	r,c,z = generate2faces(symmat)
	farfaces[2] = r
	firstv[2] = c
	filtration[2] = z
	prepairs[2] = Array(Int64,0)
	numFiltrationLevels = length(filtration[2])

	if maxcard == 3
		generate3faces!(farfaces,firstv,filtration,prepairs,m,symmat;verbose = verbose)
		if dictionaryoutput == true
			D = Dict(
				"farfaces" => farfaces, 
				"firstv" => firstv, 
				"filtration" => filtration, 
				"prepairs" => prepairs,
				"symmat" => symmat,
				"oldvertsinnewplaces"=>vperm)
			return D
		else
			return farfaces,firstv,filtration,prepairs,symmat,vperm
		end
	end

	fpi = Array(Int64,0)
	ff2pv = Array(Int64,0)
	pmhist = zeros(Int64,m,m)

	for setcard = 3:maxcard			
		if verbose
			print(["set cardinality = " setcard])
			println(["num setcard-1 cells" length(farfaces[setcard-1])])
		end
	
		nl = length(farfaces[setcard-1])
		nll = length(farfaces[setcard-2])	

		startlength = nl
		stepsize = min(10^7,Int(ceil(nl/4)))
	
		npsupp = trues(nl)
		pflist = Array(Int64,nl)
		jrv = farfaces[setcard-1]
		jcp = firstv[setcard-1]
		jz = filtration[setcard-1]	
		zll= filtration[setcard-2]
		izfull = Array(Int,nll)
		r = Array(Int64,startlength)
		z = Array(Int64,startlength)
		c = Array(Int64,m+1)	
		c[1]=1
		numpairs = [0]	
		facecount = [0]
		if setcard == maxcard-1
			ff2pv = Array(Int64,nl)
			ff2pv[:] = m+1
		end
		if setcard == maxcard			
			#### sort j-matrix by filtration		
			alterweight = Array(Int64,length(zll));
			maxweight = maximum(zll);
			for i = 1:length(alterweight)
				alterweight[i] = 1+maxweight-zll[i]
			end
			lowfilt = yafterx(alterweight,jrv)
			invertiblevec = integersinsameorderbycolumn(lowfilt,numFiltrationLevels+1,jcp)
			inversevec0 = Array(Int64,nl)
			inversevec0[invertiblevec]=1:nl
			jrv = yafterx(jrv,inversevec0)
			jz = yafterx(jz,inversevec0)								
		
			lowfilt = yafterx(ff2pv,jrv)
			invertiblevec = integersinsameorderbycolumn(lowfilt,m+1,jcp)
			inversevec1 = Array(Int64,nl)
			inversevec1[invertiblevec]=1:nl
			jrv = yafterx(jrv,inversevec1)
			jz = yafterx(jz,inversevec1)	
			translatorvecb = yafterx(inversevec0,inversevec1)
			inversevec0 = [];inversevec1 = [];lowfilt = [];invertiblevec = [];
			gc()				
			(rt,ct,zt) = transposeLighter(jrv,jcp,jz,nll) 
			colsum = ct-1		

			pmhist = zeros(Int64,m+1,m) #for sloth (apologies) we'll leave some unsed stuff in row m+1
			fpi = zeros(Int64,m+1,m)
			processfpi!(pmhist,fpi,jcp,jrv,ff2pv,m)					
		
			#### reset ff2pv for next round
			ff2pvold = copy(ff2pv)
			ff2pv = Array(Int64,nl)
			ff2pv[:] = m+1			

			oldclaw = Array(Int64,m)
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
				if setcard < maxcard-1	
					process_setcard_lt_maxcard!(
						i::Int64,j::Int64,dij::Int64,stepsize::Int64,
						facecount::Array{Int64,1},numpairs::Array{Int64,1},
						jrv::Array{Int64,1},jcp::Array{Int64,1},jz::Array{Int64,1},
						r::Array{Int64,1},z::Array{Int64,1},pflist::Array{Int64,1},					
						izfull::Array{Int64,1},ff2pv::Array{Int64,1},pmhist::Array{Int64,2},
						npsupp::BitArray{1})				
				elseif setcard == maxcard-1	
					process_setcard_onelt_maxcard_1!(
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
			#### update the column pattern and the total number of nonzeros encountered per codim2 face
			c[i+1] = facecount[1]+1
			if setcard == maxcard
				colsum[jrv[cran(jcp,i)]]+=1
			end	
		end  
		delrange = c[end]:length(r)
		deleteat!(r,delrange)
		deleteat!(z,delrange)
		deleteat!(pflist,(numpairs[1]+1):nl)
		if setcard == maxcard
			r = translatorvecb[r]	
		end
		firstv[setcard] = c
		farfaces[setcard] = r
		prepairs[setcard] = pflist
		filtration[setcard] = z
		if isempty(farfaces[setcard])
			for nextcard = (setcard+1):maxcard
				firstv[nextcard] = [1;1]
				farfaces[nextcard] = Array(Int64,0)
				prepairs[nextcard] = Array(Int64,0)
				filtration[nextcard] = Array(Int64,0)
			end
			if verbose
				println("no simplices of cardinality $(setcard) or higher")
			end
			break		
		end
	end
	if verbose
		println("collecting garbage")
		println(["number of edges" length(farfaces[2])])
	end
	gc()
	if dictionaryoutput == true
		D = Dict(
			"farfaces" => farfaces, 
			"firstv" => firstv, 
			"filtration" => filtration, 
			"prepairs" => prepairs,
			"symmat" => symmat,
			"oldvertsinnewplaces"=> vperm)
		return D
	else
		return farfaces,firstv,filtration,prepairs,symmat,vperm
	end
end

function process_setcard_lt_maxcard!(
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

function process_setcard_onelt_maxcard_1!(
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
		append!(r,Array(Int64,stepsize))
		append!(z,Array(Int64,stepsize))
	end
	r[facecount]= k
	z[facecount]= farfilt
end

function faceupdatedeluxe!(facecount::Array{Int64,1},r::Array{Int64,1},z::Array{Int64,1},k::Int64,farfilt::Int64,stepsize::Int64,s::Array{Int64,1},i::Int64)									
	facecount[1]+=1
	if facecount[1]>length(r)
		append!(r,Array(Int64,stepsize))
		append!(z,Array(Int64,stepsize))
		append!(s,Array(Int64,stepsize))
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
		rowval = Array(Int64,L)
		nzval = Array(Int64,L)		
		colptr = Array(Int64,m+1)
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
	filtration_cell,
	prepairs_cell,
	m,
	symmat;
	verbose = false)
	
	filtration::Array{Int64,1} = filtration_cell[2]
	farfaces::Array{Int64,1} = farfaces_cell[2]
	firstv::Array{Int64,1} = firstv_cell[2]
	
	numverts = length(firstv)-1
	numedges = length(farfaces)
	stepsize = 10^7
	facecount= [0]	
	numpairs = 0
	
	closefaces = Array(Int64,numedges)
	for i = 1:m
		closefaces[cran(firstv,i)]=i
	end
	iso = integersinsameorder(farfaces,m)
	closefaces_higsorted = 	Array(Int64,numedges)
	filtration_higsorted = 	Array(Int64,numedges)	
	closefaces_higsorted[iso] = closefaces
	filtration_higsorted[iso] = filtration
	
	firstv_hs = zeros(Int64,m+1)
	for face in farfaces
		firstv_hs[face+1]+=1
	end
	firstv_hs[1] = 1
	for i = 2:(m+1)
		firstv_hs[i] = firstv_hs[i-1]+firstv_hs[i]
	end
			
	adist = Array(Int64,m)
	idist = Array(Int64,m)
	r = Array(Int64,numedges)
	z = Array(Int64,numedges)
	s = Array(Int64,numedges)
	
	clawvec = Array(Int64,m)	
	ncheckedges = trues(numedges)	
	
	for a = 1:m
		adist[:]=0
		adist[crows(firstv,farfaces,a)] = crows(firstv,filtration,a)
		for ip = cran(firstv,a)
			i = farfaces[ip]
			dai = filtration[ip]	
			idist[:]=0
			idist[crows(firstv,farfaces,i)]	= crows(firstv,filtration,i)			
			idist[crows(firstv_hs,closefaces_higsorted,i)] = crows(firstv_hs,filtration_higsorted,i)
			for jp = cran(firstv,i)
				if ncheckedges[jp]
					j = farfaces[jp]
					dij = filtration[jp]
					if dij <= dai && dij <= adist[j] # note this condition bakes in the req. that j be adjacent to a 
						numpairs+=1
						ncheckedges[jp] = false
						clawvec[1:i] = 0
						for lp = cran(firstv_hs,j)
							l = closefaces_higsorted[lp]
							if l >= i
								break
							elseif idist[l]!=0
								clawvec[l] = min(idist[l],filtration_higsorted[lp])
							end															
						end
						for kp = cran(firstv,j)
							k = farfaces[kp]
							djk = filtration[kp]
							dak = adist[k]
							dik = idist[k]						
							if dak < dij && dak<djk && dak<dik	# this bakes in req. that dik>0			
								dijk = min(dij,dik,djk)
								keepface = true									
								for bp = cran(firstv_hs,k)
									b = closefaces_higsorted[bp]
									if b >= i
										break
									elseif min(clawvec[b],filtration_higsorted[bp]) >= dijk
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
		dij = filtration[edge]
		if i != holdi
			idist[:]=0
			idist[crows(firstv,farfaces,i)]	= crows(firstv,filtration,i)			
			idist[crows(firstv_hs,closefaces_higsorted,i)] = crows(firstv_hs,filtration_higsorted,i)		
			holdi = i			
		end
		clawvec[1:i] = 0
		for lp = cran(firstv_hs,j)
			l = closefaces_higsorted[lp]
			if l >= i
				break
			elseif idist[l]!=0
				clawvec[l] = min(idist[l],filtration_higsorted[lp])
			end															
		end			
		#### a facsimile of above
		for kp = cran(firstv,j)
			k = farfaces[kp]
			dik = idist[k]
			if dik==0
				continue
			end			
			djk = filtration[kp]
			dijk = min(dij,dik,djk)							
			keepface = true
			for bp = cran(firstv_hs,k)
				b = closefaces_higsorted[bp]
				if b >= i
					break
				elseif min(clawvec[b],filtration_higsorted[bp]) >= dijk		
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
	
	iso = integersinsameorder(s,m)
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
	prepairs = Array(Int64,numedges)
	for i = 1:num3faces
		edge = r[i]
		if npes[edge] && z[i] == filtration[edge]
			npes[edge]=false
			pairmarker+=1
			prepairs[pairmarker]=i
		end
	end
	deleteat!(prepairs,(pairmarker+1):numedges)	

	farfaces_cell[3]=r
	firstv_cell[3]=fv3
	filtration_cell[3]=z
	prepairs_cell[3]=prepairs
	
	buffer1 = Array{Any}(1)
	buffer2 = Array{Any}(1)
	buffer1[1] = Array(Int64,0)
	buffer2[1] = ones(Int64,numverts+1)
	append!(farfaces_cell,buffer1)
	append!(filtration_cell,buffer1)
	append!(prepairs_cell,buffer1)
	append!(firstv_cell,buffer1)

	return r,fv3,z,prepairs,numpairs
end

##########################################################################################

####	TESTING AND DIAGNOSTICS

##########################################################################################

function persistencestats(x)
	L = length(x)
	A = Array(Float64,8,L)
	for ip = 1:L
		i = x[ip]
		println(i)								
		println(i)		
		pcloud = rand(20,i)
		tic()
		D = persistF2(pcloud,5,rowsare = "dimensions",conestop=false)
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
	trials = Array{Any}(numtrials)
	for i = 1:numtrials
		trials[i] = persistencestats(x)
	end
	return trials
end

function construction_sanitycheck(;numtrials = 10,samplesize = 50,setcard=4)
	for i = 1:numtrials
		pcloud = rand(20,samplesize)
		d = Distances.pairwise(Euclidean(),pcloud)
		(t,filtrationtranslator) = ordercanonicalform(d;conestop=false)
		construction_sanitycheck_subroutine(t,setcard,samplesize)
		gc()
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

function perseuscrosscheck(;numtrials = 1,samplesize = 50,setcard=3)
	for i = 1:numtrials
		pcloud = rand(20,samplesize)
		d = Distances.pairwise(Euclidean(),pcloud)
		(t,filtrationtranslator) = ordercanonicalform(d;conestop=false)
		
		#### the following directories must be adjusted by the user
		writedlm("/Users/greghenselman/Julia/testtext.txt",""," ")
		outfile = open("/Users/greghenselman/Julia/testtext.txt","a+")		
		str = "$(samplesize)\n0.1 1 $(maximum(t)) $(setcard-1)"
		write(outfile,str)
		for i = 1:samplesize
			str = string(t[i,:])
			str = "\n"*str[2:end-1]
			write(outfile,str)			
		end
		close(outfile)
		
		tic()
		run(`/Users/greghenselman/Downloads/perseusMac distmat /Users/greghenselman/Julia/testtext.txt perseusoutput_caliber`)
		toc()
		
		p1 = readdlm("/Users/greghenselman/perseusoutput_caliber_$(setcard-2).txt")
		p1 = sortrows(p1)
		p3 = readdlm("/Users/greghenselman/perseusoutput_caliber_$(setcard-1).txt")
		p3 = sortrows(p3)		
		
		tic()
		D = persistF2(pcloud,setcard,rowsare = "dimensions",conestop=false,fillud = "down")
		E = persistF2(pcloud,setcard,rowsare = "dimensions",fillud = "down")		
		toc()
		
		unpack!(D)
		unpack!(E)		
		bcs = getbarcodesummary(D,setcard,ocf=true)
		p2 = sortrows(maximum(D["filtration"][2])+1-convert(Array{Float64},bcs[:,2:3]));
		bcs = getbarcodesummary(E,setcard,ocf=true)		
		p4 = sortrows(maximum(E["filtration"][2])+1-convert(Array{Float64},bcs[:,2:3]));		
		
		if p1 != p2 || p1 !=p4
			println("p1 != p2")
			println(size(p1))
			println(size(p2))
			println(size(p3))
			println(size(p4))						
			println("$(samplesize)\n0.1 1 $(maximum(t)) $(setcard-1)")
			println("/Users/greghenselman/Julia/perseusoutput_caliber_$(setcard-2).txt")
			break
		else
			println("p1 == p2")
		end
	end
end

##########################################################################################

####	MAIN

##########################################################################################

"""

    eirene(X[, keyword arugemts])
    
Computes the persistent homology of a Vietoris-Rips complex.

"""
function eirene(
	s;
	bettimax = 1,
	rowsare="distances",
	filetype="textfile",
	lowerlim=-Inf,
	upperlim=Inf,
	fillud="up",
	numlevels=Inf,
	conestop=true,
	record="cyclerep",
	pointlabels=[],
	verbose=false)

	maxcard = bettimax+2
	D = persistF2(
		s,
		maxcard; 
		rowsare = rowsare,
		filetype = filetype,
		minval = lowerlim,
		maxval = upperlim,
		fillud = fillud,
		pointlabels = pointlabels,
		numlevels = numlevels,
		conestop = conestop,
		record = record,
		verbose = verbose)

	return D
end
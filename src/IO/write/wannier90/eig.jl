
function Base.write(io::IO, band::Eigen, ::Type{wannier90_eig}; comment = "",
	bandindex::AbstractVector{<:Integer} = eachindex(band.values))

	now = Dates.format(Dates.now(), "yyyy/m/d H:M:S")
	comment = comment * " " * now


	for (ii, ei) in enumerate(bandindex)
		@printf(io, "%12u %12u %22.12f \n", ii, i, band.values[ei])
	end


	return nothing
end
function Base.write(io::IO, band::AbstractVector{<:Eigen}, ::Type{wannier90_eig}; comment = "",
	bandindex::AbstractVector{<:Integer} = eachindex(band[1].values))

	now = Dates.format(Dates.now(), "yyyy/m/d H:M:S")
	comment = comment * " " * now

	for (i, E) in enumerate(band)
		for (ii, ei) in enumerate(bandindex)
			@printf(io, "%12u %12u %22.12f \n", ii, i, E.values[ei])
		end
	end

	return nothing
end
function Base.write(io::IO, band::AbstractVecOrMat{<:Real}, ::Type{wannier90_eig};
	bandindex::AbstractVector{<:Integer} = axes(band, 1), comment = "")

	now = Dates.format(Dates.now(), "yyyy/m/d H:M:S")
	comment = comment * " " * now

	for ki in axes(band, 2)
		for ei in bandindex
			@printf(io, "%12u %12u %22.12f \n", ei, ki, band[ei, ki])
		end
	end

	return nothing
end

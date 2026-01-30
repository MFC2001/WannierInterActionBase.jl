
function Base.write(io::IO, band::Eigen, ::Type{BAND_dat}; kwargs...)
	bandmatrix = _eigen2val(band)
	write(io, bandmatrix, BAND_dat; kwargs...)
end
function Base.write(io::IO, band::AbstractVector{<:Eigen}, ::Type{BAND_dat}; kwargs...)
	bandmatrix = _eigen2val(band)
	write(io, bandmatrix, BAND_dat; kwargs...)
end
function Base.write(io::IO, band::AbstractVecOrMat{<:Real}, ::Type{BAND_dat}; comment = "",
	kline = Kline(; line = collect(1:size(band, 2)), name = ["Γ"], index = [1]))

	size(band, 2) == length(kline.line) || begin
		@info "Mismatched kline and band. Reset to default."
		kline = Kline(; line = collect(1:size(band, 2)), name = ["Γ"], index = [1])
	end

	Nband = size(band, 1)
	Nk = length(kline.line)

	now = Dates.format(Dates.now(), "yyyy/m/d H:M:S")
	comment = "#" * comment * " " * now

	println(io, comment * " #K-Path(1/A) Energy-Level(eV)")
	println(io, "# NKPTS & NBANDS: ", Nk, " ", Nband)

	print(io, "# kindex: ")
	for index in kline.index
		@printf(io, "%8u", index)
	end
	print(io, "\n# kpoint: ")
	for index in kline.index
		@printf(io, "%8.4f", kline.line[index])
	end
	print(io, "\n# kname:  ")
	for name in kline.name
		@printf(io, "%8s", name)
	end
	println(io, "")

	if Nk == 1
		println(io, "# Only Γ.")
		output = Vector{String}(undef, Nband)
		for i in 1:Nband
			output[i] = @sprintf("%8u%14.6f", i, band[i])
		end
		writedlm(io, output)
	else
		line = Vector{String}(undef, Nk)
		for i in 1:Nk
			line[i] = @sprintf("%10.5f", kline.line[i])
		end
		bandi = Vector{String}(undef, Nk)
		for i in 1:Nband
			for j in 1:Nk
				bandi[j] = @sprintf("%14.6f", band[i, j])
			end
			println(io, "# Band-Index  $(i)")
			writedlm(io, line .* bandi)
			println(io, "")
		end
	end

	return nothing
end

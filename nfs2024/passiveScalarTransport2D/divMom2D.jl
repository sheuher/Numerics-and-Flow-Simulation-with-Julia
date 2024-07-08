function divMom2D!(divPred, rhouP, delta)
	# calc resulting divergence field for a given momentum prediction field using CDS

	nG = 1
	ImaAll, JmaAll = size(divPred)
	Ifi = nG +1; Ifim = Ifi-1; Ifip = Ifi +1
	Ila = ImaAll - nG; Ilap = Ila+1; Ilam = Ila-1
	Jfi = nG +1; Jfim = Jfi-1; Jfip = Jfi +1
	Jla = JmaAll - nG; Jlap = Jla+1; Jlam = Jla-1

	rhouPE = @view rhouP[Ifip:Ilap, Jfi:Jla, 1]
	rhouPW = @view rhouP[Ifim:Ilam, Jfi:Jla, 1]
	rhouPN = @view rhouP[Ifi:Ila, Jfip:Jlap, 2]
	rhouPS = @view rhouP[Ifi:Ila, Jfim:Jlam, 2]

	divPred[Ifi:Ila, Jfi:Jla] .= 0.5 /delta *( rhouPE .- rhouPW .+ rhouPN .- rhouPS ) 

	return nothing
end
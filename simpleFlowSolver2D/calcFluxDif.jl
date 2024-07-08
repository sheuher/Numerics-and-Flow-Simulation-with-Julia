function calcFluxDif!(fluxDifX, fluxDifY, Phi, D, deltaX)

	#use central differencing
	
	nG = 1

	ImaAll, JmaAll = size(fluxDifX)
	Ifi = nG+1; Ifim = Ifi-1; Ifip = Ifi+1
	Ila = ImaAll -nG; Ilap = Ila+1; Ilam = Ila-1
	Jfi = nG+1; Jfim = Jfi-1; Jfip = Jfi+1
	Jla = JmaAll -nG; Jlap = Jla+1; Jlam = Jla-1

	PhiC = @view Phi[Ifi:Ila,Jfi:Jla]
	PhiE = @view Phi[Ifip:Ilap, Jfi:Jla]; PhiW = @view Phi[Ifim:Ilam, Jfi:Jla]
	PhiN = @view Phi[Ifi:Ila, Jfip:Jlap]; PhiS = @view Phi[Ifi:Ila, Jfim:Jlam]

	fluxDifX[Ifi:Ila, Jfi:Jla] .= D *deltaX * (PhiE .- 2*PhiC .+ PhiW)
	fluxDifY[Ifi:Ila, Jfi:Jla] .= D *deltaX * (PhiN .- 2*PhiC .+ PhiS)

	return nothing
end
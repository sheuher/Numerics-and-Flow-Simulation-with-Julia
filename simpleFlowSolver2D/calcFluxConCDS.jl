function calcFluxConCDS!(fluxConX, fluxConY, Phi, U, V, Dx)

	# using central differencing
	
	nG = 1

	ImaAll, JmaAll = size(fluxConX)
	Ifi = nG+1; Ifim = Ifi-1; Ifip = Ifi+1
	Ila = ImaAll -nG; Ilap = Ila+1; Ilam = Ila-1
	Jfi = nG+1; Jfim = Jfi-1; Jfip = Jfi+1
	Jla = JmaAll -nG; Jlap = Jla+1; Jlam = Jla-1

	PhiC = @view Phi[Ifi:Ila, Jfi:Jla]
	PhiE = @view Phi[Ifip:Ilap, Jfi:Jla]; PhiW = @view Phi[Ifim:Ilam, Jfi:Jla]
	PhiN = @view Phi[Ifi:Ila, Jfip:Jlap]; PhiS = @view Phi[Ifi:Ila, Jfim:Jlam]

	Ue = @view U[Ifi:Ila, Jfi:Jla]; Uw = @view U[Ifim:Ilam, Jfi:Jla]
	Vn = @view V[Ifi:Ila, Jfi:Jla]; Vs = @view V[Ifi:Ila, Jfim:Jlam]

	fluxConX[Ifi:Ila, Jfi:Jla] .= 0.5 *Dx^2 *((PhiE +PhiC) .*Ue .- (PhiC +PhiW) .*Uw)

	fluxConY[Ifi:Ila, Jfi:Jla] .= 0.5 *Dx^2 *((PhiN +PhiC) .*Vn .- (PhiC + PhiS) .*Vs)
	
	return nothing
end
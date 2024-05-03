function poissonSolver2D!(p, divPred, dt, delta, nItMax, epsMax)
	# solve poisson eq. with iterative jacobi method for given divergence field, apply zero gradient BC at the inlet and fix pressure = 0 at outlet

	nG = 1
	ImaAll, JmaAll = size(divPred)
	Ifi = nG +1; Ifim = Ifi-1; Ifip = Ifi +1
	Ila = ImaAll - nG; Ilap = Ila+1; Ilam = Ila-1
	Jfi = nG +1; Jfim = Jfi-1; Jfip = Jfi +1
	Jla = JmaAll - nG; Jlap = Jla+1; Jlam = Jla-1
	
	RHS = divPred /dt
	nIt = 0; eps = 10.0
	### while ( (nIt <nItMax) && (eps>epsMax) )
	for _ in 1:nItMax
		pold = copy(p)
		pE = @view pold[Ifip:Ilap, Jfi:Jla]
		pW = @view pold[Ifim:Ilam, Jfi:Jla]
		pN = @view pold[Ifi:Ila, Jfip:Jlap]
		pS = @view pold[Ifi:Ila, Jfim:Jlam]
		p[Ifi:Ila, Jfi:Jla] = 0.25 *(  pE + pW + pN + pS -RHS[Ifi:Ila, Jfi:Jla]*delta^2 )
		p[Ifim, Jfi:Jla] .= p[Ifi, Jfi:Jla]
		eps = maximum( abs.(p[Ifi:Ila] -pold[Ifi:Ila]) )
		#if eps < epsMax
		#	break
		#end
		nIt = nIt + 1
	end
	return (nIt=nIt, eps=eps)
end
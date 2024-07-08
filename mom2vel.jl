function mom2vel!(u, v, rhou, rhov)

	"""
	calculate velocity fields u[:,:] and v[:,:] on the CELL SURFACES
	from the momentum fields rhou[:,:] and rhov[:,:] stored on cell centre
	"""
	nG = 1

	ImaAll, JmaAll = size(u)
	Ifi = nG+1; Ifim = Ifi-1; Ifip = Ifi+1
	Ila = ImaAll -nG; Ilap = Ila+1; Ilam = Ila-1
	Jfi = nG+1; Jfim = Jfi-1; Jfip = Jfi+1
	Jla = JmaAll -nG; Jlap = Jla+1; Jlam = Jla-1

	u[Ifim:Ila, Jfi:Jla] .= 0.5 * (
		rhou[Ifi:Ilap, Jfi:Jla] + rhou[Ifim:Ila, Jfi:Jla]
	)

	v[Ifi:Ila, Jfim:Jla] .= 0.5 * (
		rhov[Ifi:Ila, Jfi:Jlap] + rhov[Ifi:Ila, Jfim:Jla]
	)

	return nothing
end
function applyBCV!(Phi)
    """
        apply ... boundary condition to pressure
    """
    nG = 1
    ImaAll, JmaAll = size(Phi)
	Ifi = nG+1; Ifip = Ifi+1; Ifim = Ifi-1
	Ila = ImaAll - nG; Ilap = Ila+1; Ilam = Ila-1
	Jfi = nG+1; Jfip = Jfi+1; Jfim = Jfi-1
	Jla = JmaAll - nG; Jlap = Jla+1; Jlam = Jla-1
    
    Phi[Ifim, Jfi:Jla] .= 0
    Phi[Ilap, Jfi:Jla] .= 0
    Phi[Ifi:Ila, Jfim] .= 0
    Phi[Ifi:Ila, Jlap] .= 0
end